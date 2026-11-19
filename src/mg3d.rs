//! Geometric multigrid preconditioner for the 3D voxel solver (`stress3d`).
//!
//! The fine level is the voxel system the solver already has (one shared brick matrix times a
//! stiffness per voxel; void voxels have stiffness 0). Coarse levels halve each axis, rounding up
//! (so odd grids work: the coarse grid then reaches half a cell past the fine one). Their cell
//! matrices are the Galerkin product P^T K P with trilinear interpolation, computed per coarse
//! cell: interpolation inside a coarse cell only uses that cell's eight corners, so summing
//! `P_f^T K_f P_f` over the (up to eight) fine voxels grouped under a coarse cell is exactly the
//! global product. That keeps solid/void and skin/infill contrasts on the coarse grids, which
//! averaging stiffness would blur.
//!
//! The V-cycle uses damped Jacobi sweeps (the same before and after, so the preconditioner is
//! symmetric) and a dense Cholesky solve on the coarsest level. Fixed DOFs are masks: operators
//! return zero there and transfers skip them.
//!
//! Native-only code (the plugin never runs it): indexing is by grid dimensions throughout and the
//! tests check every operator against a slow reference, so the crate's no-indexing lint is off.
#![allow(clippy::indexing_slicing)]

use crate::stress3d::{CORNERS, Mat24};

/// What the hierarchy needs to know about the fine system.
pub struct FineSystem<'a> {
    /// Voxels per axis.
    pub n: [usize; 3],
    /// The brick matrix of a unit voxel; a voxel's stiffness is `ke * cell_k[c]`.
    pub ke: &'a Mat24,
    /// Per voxel (grid order `(i * ny + j) * nz + k`): stiffness factor, 0 = no material.
    pub cell_k: &'a [f64],
    /// Per DOF: held (supports, and nodes no material touches).
    pub fixed: &'a [bool],
    /// Per DOF: held by a support (not just outside the material).
    pub support: &'a [bool],
}

/// Settings of the V-cycle.
#[derive(Clone, Copy, Debug)]
pub struct Settings {
    /// Jacobi sweeps before and after the coarse correction.
    pub sweeps: usize,
    /// Jacobi damping as a fraction of 1 / (largest eigenvalue of D^-1 K), estimated per level:
    /// a fixed damping such as 0.6 nearly diverges on thin, uneven parts.
    pub damping: f64,
    /// Stop coarsening when at most this many DOFs are free (dense Cholesky on them).
    pub coarsest_dofs: usize,
    /// Refuse to build when the first coarse level would need more than this many bytes.
    pub max_bytes: usize,
}

impl Default for Settings {
    fn default() -> Self {
        Self { sweeps: 2, damping: 1.4, coarsest_dofs: 1500, max_bytes: 4_000_000_000 }
    }
}

fn node_id(n: &[usize; 3], i: usize, j: usize, k: usize) -> usize {
    (i * (n[1] + 1) + j) * (n[2] + 1) + k
}

fn cell_id(n: &[usize; 3], i: usize, j: usize, k: usize) -> usize {
    (i * n[1] + j) * n[2] + k
}

fn dofs_of(n: &[usize; 3], i: usize, j: usize, k: usize) -> [usize; 24] {
    let mut out = [0usize; 24];
    for (q, c) in CORNERS.iter().enumerate() {
        let node = node_id(n, i + c[0], j + c[1], k + c[2]);
        for d in 0..3 {
            out[3 * q + d] = 3 * node + d;
        }
    }
    out
}

/// A coarse level: a matrix per coarse cell that holds any stiffness.
#[derive(Debug)]
pub struct Level {
    pub n: [usize; 3],
    /// Cell coordinates of each stored matrix.
    cells: Vec<[usize; 3]>,
    /// Position of each cell's matrix in `mats` (`u32::MAX` = none).
    index: Vec<u32>,
    mats: Vec<Mat24>,
    pub fixed: Vec<bool>,
    inv: Vec<f64>,
}

impl Level {
    pub fn ndof(&self) -> usize {
        3 * (self.n[0] + 1) * (self.n[1] + 1) * (self.n[2] + 1)
    }

    fn free_dofs(&self) -> usize {
        self.fixed.iter().filter(|f| !**f).count()
    }

    /// y = K u (zero on fixed DOFs).
    pub fn apply(&self, u: &[f64], y: &mut [f64]) {
        y.iter_mut().for_each(|v| *v = 0.0);
        for ([i, j, k], k24) in self.cells.iter().zip(&self.mats) {
            let d = dofs_of(&self.n, *i, *j, *k);
            let ue = d.map(|x| u[x]);
            let mut ye = [0.0f64; 24];
            for (row, v) in k24.iter().zip(ue) {
                for (out, kv) in ye.iter_mut().zip(row) {
                    *out += kv * v;
                }
            }
            for (x, v) in d.iter().zip(ye) {
                y[*x] += v;
            }
        }
        for (v, f) in y.iter_mut().zip(&self.fixed) {
            if *f {
                *v = 0.0;
            }
        }
    }
}

/// The matrices of the level being coarsened: the fine system or an earlier coarse level.
enum Source<'a> {
    Fine(&'a FineSystem<'a>),
    Coarse(&'a Level),
}

impl Source<'_> {
    fn n(&self) -> [usize; 3] {
        match self {
            Source::Fine(f) => f.n,
            Source::Coarse(l) => l.n,
        }
    }
    fn fixed(&self, dof: usize) -> bool {
        match self {
            Source::Fine(f) => f.fixed[dof],
            Source::Coarse(l) => l.fixed[dof],
        }
    }
    /// A coarse DOF is held where the coincident fine DOF is held by a support (or, from a coarse
    /// level, held at all); nodes outside the material do not hold their coarse DOF.
    fn holds_coarse(&self, dof: usize) -> bool {
        match self {
            Source::Fine(f) => f.support[dof],
            Source::Coarse(l) => l.fixed[dof],
        }
    }
    fn cell(&self, i: usize, j: usize, k: usize) -> Option<(&Mat24, f64)> {
        match self {
            Source::Fine(f) => {
                let s = f.cell_k[cell_id(&f.n, i, j, k)];
                (s > 0.0).then_some((f.ke, s))
            }
            Source::Coarse(l) => {
                let at = l.index[cell_id(&l.n, i, j, k)];
                (at != u32::MAX).then(|| (&l.mats[at as usize], 1.0))
            }
        }
    }
}

type Entries = Vec<(usize, usize, f64)>;

/// The non-zeros of the interpolation of one fine voxel (at offset (a, b, d) inside its coarse
/// cell) from the coarse cell's corners: (fine DOF, coarse DOF, weight).
fn interpolation_entries(a: usize, b: usize, d: usize) -> Entries {
    let mut list = Vec::new();
    for (q, fc) in CORNERS.iter().enumerate() {
        let pos = [a + fc[0], b + fc[1], d + fc[2]];
        for (r, cc) in CORNERS.iter().enumerate() {
            let w: f64 = (0..3).map(|ax| (1.0 - (pos[ax] as f64 - 2.0 * cc[ax] as f64).abs() / 2.0).max(0.0)).product();
            if w > 0.0 {
                for c in 0..3 {
                    list.push((3 * q + c, 3 * r + c, w));
                }
            }
        }
    }
    list
}

/// kc += scale * P^T K P for the given interpolation entries.
fn accumulate(kc: &mut Mat24, k: &Mat24, scale: f64, entries: &[(usize, usize, f64)]) {
    let mut kp = [[0.0f64; 24]; 24]; // K P
    for &(fi, ci, w) in entries {
        for (kprow, krow) in kp.iter_mut().zip(k) {
            kprow[ci] += krow[fi] * w;
        }
    }
    for &(fi, ci, w) in entries {
        let f = w * scale;
        for (out, v) in kc[ci].iter_mut().zip(&kp[fi]) {
            *out += f * v;
        }
    }
}

fn coarsen(src: &Source, settings: &Settings) -> Level {
    let n = src.n();
    let cn = n.map(|v| v.div_ceil(2));
    let nodes = (cn[0] + 1) * (cn[1] + 1) * (cn[2] + 1);

    // A coarse DOF is held where its node coincides with a held fine node.
    let mut fixed = vec![false; 3 * nodes];
    for i in 0..=cn[0] {
        for j in 0..=cn[1] {
            for k in 0..=cn[2] {
                let (fi, fj, fk) = (2 * i, 2 * j, 2 * k);
                if fi <= n[0] && fj <= n[1] && fk <= n[2] {
                    let fine_node = node_id(&n, fi, fj, fk);
                    for c in 0..3 {
                        if src.holds_coarse(3 * fine_node + c) {
                            fixed[3 * node_id(&cn, i, j, k) + c] = true;
                        }
                    }
                }
            }
        }
    }

    let table: Vec<Entries> = (0..8).map(|w| interpolation_entries((w >> 2) & 1, (w >> 1) & 1, w & 1)).collect();
    // From the shared brick matrix, the eight products P^T Ke P are the same everywhere.
    let shared: Option<Vec<Mat24>> = match src {
        Source::Fine(f) => Some(
            table
                .iter()
                .map(|list| {
                    let mut g = [[0.0f64; 24]; 24];
                    accumulate(&mut g, f.ke, 1.0, list);
                    g
                })
                .collect(),
        ),
        Source::Coarse(_) => None,
    };

    let mut cells = Vec::new();
    let mut mats: Vec<Mat24> = Vec::new();
    let mut index = vec![u32::MAX; cn[0] * cn[1] * cn[2]];
    let mut filtered: Entries = Vec::with_capacity(81);
    for ci in 0..cn[0] {
        for cj in 0..cn[1] {
            for ck in 0..cn[2] {
                // The fine voxels grouped under this coarse cell that hold material.
                let mut present: Vec<(usize, [usize; 3], &Mat24, f64)> = Vec::with_capacity(8);
                for w in 0..8usize {
                    let (a, b, d) = ((w >> 2) & 1, (w >> 1) & 1, w & 1);
                    let (fi, fj, fk) = (2 * ci + a, 2 * cj + b, 2 * ck + d);
                    if fi < n[0] && fj < n[1] && fk < n[2] {
                        if let Some((k, s)) = src.cell(fi, fj, fk) {
                            present.push((w, [fi, fj, fk], k, s));
                        }
                    }
                }
                if present.is_empty() {
                    continue;
                }
                let corner_dofs = dofs_of(&cn, ci, cj, ck);
                let touches = corner_dofs.iter().any(|d| fixed[*d]) || present.iter().any(|(_, [fi, fj, fk], _, _)| dofs_of(&n, *fi, *fj, *fk).iter().any(|d| src.fixed(*d)));
                let mut kc = [[0.0f64; 24]; 24];
                if let (false, Some(g)) = (touches, &shared) {
                    for (w, _, _, s) in &present {
                        for (out_row, row) in kc.iter_mut().zip(&g[*w]) {
                            for (out, v) in out_row.iter_mut().zip(row) {
                                *out += s * v;
                            }
                        }
                    }
                } else {
                    for (w, [fi, fj, fk], k, s) in &present {
                        let fine = dofs_of(&n, *fi, *fj, *fk);
                        filtered.clear();
                        filtered.extend(table[*w].iter().filter(|(f, c, _)| !src.fixed(fine[*f]) && !fixed[corner_dofs[*c]]));
                        accumulate(&mut kc, k, *s, &filtered);
                    }
                }
                index[cell_id(&cn, ci, cj, ck)] = mats.len() as u32;
                cells.push([ci, cj, ck]);
                mats.push(kc);
            }
        }
    }

    // Diagonal; a DOF with no stiffness at all (nothing interpolates into it) is held.
    let mut diag = vec![0.0f64; 3 * nodes];
    for ([i, j, k], m) in cells.iter().zip(&mats) {
        for (x, d) in dofs_of(&cn, *i, *j, *k).iter().zip(0..24) {
            diag[*x] += m[d][d];
        }
    }
    for (f, d) in fixed.iter_mut().zip(&diag) {
        if !(d.is_finite() && *d > 0.0) {
            *f = true;
        }
    }
    let inv = diag.iter().zip(&fixed).map(|(d, f)| if *f { 1.0 } else { 1.0 / d }).collect();
    let _ = settings;
    Level { n: cn, cells, index, mats, fixed, inv }
}

/// The coarse nodes (indices, weights) a fine node index interpolates from.
fn parents(i: usize) -> ([usize; 2], [f64; 2], usize) {
    if i % 2 == 0 { ([i / 2, 0], [1.0, 0.0], 1) } else { ([i / 2, i / 2 + 1], [0.5, 0.5], 2) }
}

/// Trilinear interpolation of a coarse vector, added to `fine` (not on fixed fine DOFs).
fn prolong_add(nf: &[usize; 3], fixed_fine: &[bool], nc: &[usize; 3], coarse: &[f64], fine: &mut [f64]) {
    for i in 0..=nf[0] {
        let (pi, wi, ci) = parents(i);
        for j in 0..=nf[1] {
            let (pj, wj, cj) = parents(j);
            for k in 0..=nf[2] {
                let (pk, wk, ck) = parents(k);
                let node = node_id(nf, i, j, k);
                let mut v = [0.0f64; 3];
                for a in 0..ci {
                    for b in 0..cj {
                        for c in 0..ck {
                            let w = wi[a] * wj[b] * wk[c];
                            let cnode = node_id(nc, pi[a], pj[b], pk[c]);
                            for (d, out) in v.iter_mut().enumerate() {
                                *out += w * coarse[3 * cnode + d];
                            }
                        }
                    }
                }
                for (d, x) in v.iter().enumerate() {
                    if !fixed_fine[3 * node + d] {
                        fine[3 * node + d] += x;
                    }
                }
            }
        }
    }
}

/// The transpose of [`prolong_add`]: a fine residual to the coarse grid (zero on fixed DOFs).
fn restrict(nf: &[usize; 3], fixed_fine: &[bool], nc: &[usize; 3], fixed_coarse: &[bool], fine: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0f64; 3 * (nc[0] + 1) * (nc[1] + 1) * (nc[2] + 1)];
    for i in 0..=nf[0] {
        let (pi, wi, ci) = parents(i);
        for j in 0..=nf[1] {
            let (pj, wj, cj) = parents(j);
            for k in 0..=nf[2] {
                let (pk, wk, ck) = parents(k);
                let node = node_id(nf, i, j, k);
                for d in 0..3 {
                    if fixed_fine[3 * node + d] {
                        continue;
                    }
                    let v = fine[3 * node + d];
                    for a in 0..ci {
                        for b in 0..cj {
                            for c in 0..ck {
                                out[3 * node_id(nc, pi[a], pj[b], pk[c]) + d] += wi[a] * wj[b] * wk[c] * v;
                            }
                        }
                    }
                }
            }
        }
    }
    for (v, f) in out.iter_mut().zip(fixed_coarse) {
        if *f {
            *v = 0.0;
        }
    }
    out
}

/// Dense Cholesky factor of the coarsest level's free DOFs.
#[derive(Debug)]
struct Dense {
    /// Free DOF number for each DOF (`usize::MAX` = held).
    map: Vec<usize>,
    m: usize,
    l: Vec<f64>,
}

impl Dense {
    fn factor(level: &Level) -> Dense {
        let mut map = vec![usize::MAX; level.ndof()];
        let mut m = 0;
        for (d, f) in level.fixed.iter().enumerate() {
            if !*f {
                map[d] = m;
                m += 1;
            }
        }
        let mut a = vec![0.0f64; m * m];
        for ([i, j, k], mat) in level.cells.iter().zip(&level.mats) {
            let d = dofs_of(&level.n, *i, *j, *k);
            for (r, dr) in d.iter().enumerate() {
                let Some(row) = Some(map[*dr]).filter(|v| *v != usize::MAX) else { continue };
                for (c, dc) in d.iter().enumerate() {
                    let col = map[*dc];
                    if col != usize::MAX {
                        a[row * m + col] += mat[r][c];
                    }
                }
            }
        }
        // In-place lower Cholesky, row by row.
        for i in 0..m {
            for j in 0..=i {
                let s: f64 = a[i * m + j] - a[i * m..i * m + j].iter().zip(&a[j * m..j * m + j]).map(|(x, y)| x * y).sum::<f64>();
                a[i * m + j] = if i == j {
                    if s.is_finite() && s > 0.0 { s.sqrt() } else { 1.0 }
                } else {
                    let dj = a[j * m + j];
                    if dj != 0.0 { s / dj } else { 0.0 }
                };
            }
        }
        Dense { map, m, l: a }
    }

    fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        let m = self.m;
        let mut y = vec![0.0f64; m];
        for (d, v) in rhs.iter().enumerate() {
            if self.map[d] != usize::MAX {
                y[self.map[d]] = *v;
            }
        }
        for i in 0..m {
            let s: f64 = y[i] - self.l[i * m..i * m + i].iter().zip(&y[..i]).map(|(a, b)| a * b).sum::<f64>();
            y[i] = s / self.l[i * m + i];
        }
        for i in (0..m).rev() {
            y[i] /= self.l[i * m + i];
            let yi = y[i];
            for (v, lv) in y[..i].iter_mut().zip(&self.l[i * m..i * m + i]) {
                *v -= lv * yi;
            }
        }
        let mut z = vec![0.0f64; rhs.len()];
        for (d, slot) in z.iter_mut().enumerate() {
            if self.map[d] != usize::MAX {
                *slot = y[self.map[d]];
            }
        }
        z
    }
}

/// The fine level's operator and inverse diagonal (the solver's own, so the fine level costs no
/// extra memory).
pub struct FineOps<'a> {
    pub apply: &'a dyn Fn(&[f64], &mut [f64]),
    pub inv: &'a [f64],
}

/// The multigrid preconditioner: coarse levels, the coarsest factor, and the V-cycle.
#[derive(Debug)]
pub struct Hierarchy3 {
    n: [usize; 3],
    fixed: Vec<bool>,
    levels: Vec<Level>,
    coarsest: Dense,
    settings: Settings,
    /// Jacobi damping of level 0 (fine), 1, ...
    omega: Vec<f64>,
}

/// The largest eigenvalue of D^-1 K on the free DOFs, by power iteration (with a small margin).
fn largest_eigenvalue(apply: &dyn Fn(&[f64], &mut [f64]), inv: &[f64], fixed: &[bool]) -> f64 {
    let n = fixed.len();
    let mut seed = 0x2545F4914F6CDD1Du64;
    let mut v: Vec<f64> = (0..n)
        .map(|i| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            if fixed[i] { 0.0 } else { (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5 }
        })
        .collect();
    let mut w = vec![0.0f64; n];
    let mut lambda = 1.0;
    for _ in 0..20 {
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm == 0.0 {
            return 1.0;
        }
        v.iter_mut().for_each(|x| *x /= norm);
        apply(&v, &mut w);
        w.iter_mut().zip(inv).for_each(|(x, d)| *x *= d);
        lambda = w.iter().map(|x| x * x).sum::<f64>().sqrt();
        std::mem::swap(&mut v, &mut w);
    }
    (lambda * 1.05).max(f64::MIN_POSITIVE)
}

impl Hierarchy3 {
    /// Builds the coarse levels. `Err` when the system is too small to need multigrid or the first
    /// coarse level would not fit in `settings.max_bytes`.
    pub fn build(fine: &FineSystem, settings: Settings, ops: &FineOps) -> Result<Hierarchy3, String> {
        let free = fine.fixed.iter().filter(|f| !**f).count();
        if free <= settings.coarsest_dofs {
            return Err("the system is small enough for the plain solver".into());
        }
        // Memory of the first coarse level: a 24 x 24 matrix per coarse cell with material.
        let n = fine.n;
        let cn = n.map(|v| v.div_ceil(2));
        let mut used = vec![false; cn[0] * cn[1] * cn[2]];
        for i in 0..n[0] {
            for j in 0..n[1] {
                for k in 0..n[2] {
                    if fine.cell_k[cell_id(&n, i, j, k)] > 0.0 {
                        used[cell_id(&cn, i / 2, j / 2, k / 2)] = true;
                    }
                }
            }
        }
        let bytes = used.iter().filter(|u| **u).count() * std::mem::size_of::<Mat24>();
        if bytes > settings.max_bytes {
            return Err(format!("the first coarse level would need {:.1} GB", bytes as f64 / 1e9));
        }
        let mut levels = vec![coarsen(&Source::Fine(fine), &settings)];
        while let Some(last) = levels.last() {
            if last.free_dofs() <= settings.coarsest_dofs || last.n == [1, 1, 1] {
                break;
            }
            let next = coarsen(&Source::Coarse(last), &settings);
            levels.push(next);
        }
        let coarsest = Dense::factor(levels.last().ok_or("no coarse level")?);
        let mut omega = vec![settings.damping / largest_eigenvalue(ops.apply, ops.inv, fine.fixed)];
        for lv in &levels {
            omega.push(settings.damping / largest_eigenvalue(&|u, y| lv.apply(u, y), &lv.inv, &lv.fixed));
        }
        Ok(Hierarchy3 { n, fixed: fine.fixed.to_vec(), levels, coarsest, settings, omega })
    }

    pub fn depth(&self) -> usize {
        self.levels.len() + 1
    }

    /// Bytes of the stored coarse matrices.
    pub fn bytes(&self) -> usize {
        self.levels.iter().map(|l| l.mats.len() * std::mem::size_of::<Mat24>()).sum()
    }

    /// z = M^-1 r: one V-cycle from zero.
    pub fn precondition(&self, ops: &FineOps, r: &[f64]) -> Vec<f64> {
        self.cycle(0, ops, r)
    }

    fn level_n(&self, l: usize) -> &[usize; 3] {
        if l == 0 { &self.n } else { &self.levels[l - 1].n }
    }

    fn level_fixed(&self, l: usize) -> &[bool] {
        if l == 0 { &self.fixed } else { &self.levels[l - 1].fixed }
    }

    fn smooth(&self, l: usize, ops: &FineOps, r: &[f64], z: &mut [f64], scratch: &mut [f64], from_zero: bool) {
        let inv: &[f64] = if l == 0 { ops.inv } else { &self.levels[l - 1].inv };
        let omega = self.omega.get(l).copied().unwrap_or(0.5);
        let mut left = self.settings.sweeps;
        if from_zero && left > 0 {
            for ((zi, ri), di) in z.iter_mut().zip(r).zip(inv) {
                *zi = omega * di * ri;
            }
            left -= 1;
        }
        for _ in 0..left {
            if l == 0 {
                (ops.apply)(z, scratch);
            } else {
                self.levels[l - 1].apply(z, scratch);
            }
            for ((zi, (ri, ki)), di) in z.iter_mut().zip(r.iter().zip(scratch.iter())).zip(inv) {
                *zi += omega * di * (ri - ki);
            }
        }
    }

    fn cycle(&self, l: usize, ops: &FineOps, r: &[f64]) -> Vec<f64> {
        if l == self.levels.len() {
            return self.coarsest.solve(r);
        }
        let n = r.len();
        let mut z = vec![0.0f64; n];
        let mut scratch = vec![0.0f64; n];
        self.smooth(l, ops, r, &mut z, &mut scratch, true);
        if l == 0 {
            (ops.apply)(&z, &mut scratch);
        } else {
            self.levels[l - 1].apply(&z, &mut scratch);
        }
        let residual: Vec<f64> = r.iter().zip(&scratch).map(|(a, b)| a - b).collect();
        let coarse = &self.levels[l];
        let rc = restrict(self.level_n(l), self.level_fixed(l), &coarse.n, &coarse.fixed, &residual);
        let zc = self.cycle(l + 1, ops, &rc);
        prolong_add(self.level_n(l), self.level_fixed(l), &coarse.n, &zc, &mut z);
        self.smooth(l, ops, r, &mut z, &mut scratch, false);
        z
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
    use super::*;
    use crate::stress3d::brick_stiffness;

    struct Block {
        n: [usize; 3],
        ke: Box<Mat24>,
        cell_k: Vec<f64>,
        fixed: Vec<bool>,
        support: Vec<bool>,
    }

    /// A box of voxels with a hole (void voxels) and soft voxels, held on its x = 0 face.
    fn block(n: [usize; 3]) -> Block {
        let ke = brick_stiffness(0.3);
        let mut cell_k = vec![1.0; n[0] * n[1] * n[2]];
        for i in 0..n[0] {
            for j in 0..n[1] {
                for k in 0..n[2] {
                    let c = cell_id(&n, i, j, k);
                    if i == n[0] / 2 && j == n[1] / 2 {
                        cell_k[c] = 0.0; // a void channel
                    } else if (i + j + k) % 7 == 0 {
                        cell_k[c] = 1e-3; // soft infill
                    }
                }
            }
        }
        let nodes = (n[0] + 1) * (n[1] + 1) * (n[2] + 1);
        let mut active = vec![false; nodes];
        for i in 0..n[0] {
            for j in 0..n[1] {
                for k in 0..n[2] {
                    if cell_k[cell_id(&n, i, j, k)] > 0.0 {
                        for c in CORNERS {
                            active[node_id(&n, i + c[0], j + c[1], k + c[2])] = true;
                        }
                    }
                }
            }
        }
        let mut fixed = vec![false; 3 * nodes];
        let mut support = vec![false; 3 * nodes];
        for i in 0..=n[0] {
            for j in 0..=n[1] {
                for k in 0..=n[2] {
                    let node = node_id(&n, i, j, k);
                    for d in 0..3 {
                        if i == 0 && active[node] {
                            support[3 * node + d] = true;
                        }
                        fixed[3 * node + d] = !active[node] || (i == 0 && active[node]);
                    }
                }
            }
        }
        Block { n, ke, cell_k, fixed, support }
    }

    impl Block {
        fn system(&self) -> FineSystem<'_> {
            FineSystem { n: self.n, ke: &self.ke, cell_k: &self.cell_k, fixed: &self.fixed, support: &self.support }
        }
        fn apply(&self, u: &[f64], y: &mut [f64]) {
            y.iter_mut().for_each(|v| *v = 0.0);
            for i in 0..self.n[0] {
                for j in 0..self.n[1] {
                    for k in 0..self.n[2] {
                        let s = self.cell_k[cell_id(&self.n, i, j, k)];
                        if s == 0.0 {
                            continue;
                        }
                        let d = dofs_of(&self.n, i, j, k);
                        let ue = d.map(|x| u[x]);
                        for (r, dr) in d.iter().enumerate() {
                            y[*dr] += s * self.ke[r].iter().zip(&ue).map(|(a, b)| a * b).sum::<f64>();
                        }
                    }
                }
            }
            for (v, f) in y.iter_mut().zip(&self.fixed) {
                if *f {
                    *v = 0.0;
                }
            }
        }
        fn inv_diag(&self) -> Vec<f64> {
            let mut diag = vec![0.0; self.fixed.len()];
            for i in 0..self.n[0] {
                for j in 0..self.n[1] {
                    for k in 0..self.n[2] {
                        let s = self.cell_k[cell_id(&self.n, i, j, k)];
                        for (d, x) in dofs_of(&self.n, i, j, k).iter().enumerate() {
                            diag[*x] += s * self.ke[d][d];
                        }
                    }
                }
            }
            diag.iter().zip(&self.fixed).map(|(d, f)| if *f || *d <= 0.0 { 1.0 } else { 1.0 / d }).collect()
        }
    }

    fn pseudo(n: usize, seed: u64, fixed: &[bool]) -> Vec<f64> {
        let mut s = seed;
        (0..n)
            .map(|i| {
                s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                if fixed[i] { 0.0 } else { ((s >> 33) as f64 / (1u64 << 31) as f64) - 0.5 }
            })
            .collect()
    }

    /// The coarse cell matrices must equal P^T K P done the slow way, on odd and even grids.
    #[test]
    fn coarse_operator_is_the_galerkin_product() {
        for n in [[8, 6, 4], [7, 5, 5], [5, 4, 3]] {
            let b = block(n);
            let level = coarsen(&Source::Fine(&b.system()), &Settings::default());
            let v = pseudo(level.ndof(), 7, &level.fixed);
            let mut pv = vec![0.0; b.fixed.len()];
            prolong_add(&n, &b.fixed, &level.n, &v, &mut pv);
            let mut kpv = vec![0.0; b.fixed.len()];
            b.apply(&pv, &mut kpv);
            let slow = restrict(&n, &b.fixed, &level.n, &level.fixed, &kpv);
            let mut fast = vec![0.0; level.ndof()];
            level.apply(&v, &mut fast);
            for (a, c) in slow.iter().zip(&fast) {
                assert!((a - c).abs() < 1e-9 * (1.0 + a.abs()), "{n:?}: {a} vs {c}");
            }
            // A second coarsening is also the Galerkin product of the first.
            let level2 = coarsen(&Source::Coarse(&level), &Settings::default());
            let w = pseudo(level2.ndof(), 9, &level2.fixed);
            let mut pw = vec![0.0; level.ndof()];
            prolong_add(&level.n, &level.fixed, &level2.n, &w, &mut pw);
            let mut kpw = vec![0.0; level.ndof()];
            level.apply(&pw, &mut kpw);
            let slow2 = restrict(&level.n, &level.fixed, &level2.n, &level2.fixed, &kpw);
            let mut fast2 = vec![0.0; level2.ndof()];
            level2.apply(&w, &mut fast2);
            for (a, c) in slow2.iter().zip(&fast2) {
                assert!((a - c).abs() < 1e-9 * (1.0 + a.abs()), "{n:?} level 2: {a} vs {c}");
            }
        }
    }

    #[test]
    fn coarse_matrices_are_symmetric() {
        let b = block([6, 5, 4]);
        let level = coarsen(&Source::Fine(&b.system()), &Settings::default());
        for m in &level.mats {
            for (i, row) in m.iter().enumerate() {
                for (j, v) in row.iter().enumerate() {
                    assert!((v - m[j][i]).abs() < 1e-9 * (1.0 + v.abs()));
                }
            }
        }
    }

    /// The V-cycle is a symmetric positive-definite preconditioner: <Mx, y> = <x, My> and <Mx, x> > 0.
    #[test]
    fn the_v_cycle_is_symmetric_and_positive() {
        let b = block([12, 10, 9]);
        let settings = Settings { coarsest_dofs: 200, ..Settings::default() };
        let inv = b.inv_diag();
        let ap = |u: &[f64], y: &mut [f64]| b.apply(u, y);
        let ops = FineOps { apply: &ap, inv: &inv };
        let h = Hierarchy3::build(&b.system(), settings, &ops).unwrap();
        assert!(h.depth() >= 3, "depth {}", h.depth());
        let x = pseudo(b.fixed.len(), 1, &b.fixed);
        let y = pseudo(b.fixed.len(), 2, &b.fixed);
        let mx = h.precondition(&ops, &x);
        let my = h.precondition(&ops, &y);
        let dot = |a: &[f64], c: &[f64]| a.iter().zip(c).map(|(p, q)| p * q).sum::<f64>();
        let (l, r) = (dot(&mx, &y), dot(&x, &my));
        assert!((l - r).abs() < 1e-9 * (l.abs() + r.abs() + 1e-30), "{l} vs {r}");
        assert!(dot(&mx, &x) > 0.0);
    }

    #[test]
    fn small_systems_are_left_to_the_plain_solver() {
        let b = block([3, 3, 3]);
        let (inv, ap) = (b.inv_diag(), |u: &[f64], y: &mut [f64]| b.apply(u, y));
        assert!(Hierarchy3::build(&b.system(), Settings::default(), &FineOps { apply: &ap, inv: &inv }).is_err());
    }

    #[test]
    fn a_memory_limit_refuses_to_build() {
        let b = block([12, 10, 9]);
        let settings = Settings { coarsest_dofs: 200, max_bytes: 1000, ..Settings::default() };
        let (inv, ap) = (b.inv_diag(), |u: &[f64], y: &mut [f64]| b.apply(u, y));
        assert!(Hierarchy3::build(&b.system(), settings, &FineOps { apply: &ap, inv: &inv }).unwrap_err().contains("GB"));
    }

    fn cg_steps(b: &Block, h: Option<&Hierarchy3>, f: &[f64]) -> usize {
        let inv = b.inv_diag();
        let ap = |u: &[f64], y: &mut [f64]| b.apply(u, y);
        let dot = |a: &[f64], c: &[f64]| a.iter().zip(c).map(|(p, q)| p * q).sum::<f64>();
        let pre = |r: &[f64]| -> Vec<f64> { match h { Some(h) => h.precondition(&FineOps { apply: &ap, inv: &inv }, r), None => r.iter().zip(&inv).map(|(a, d)| a * d).collect() } };
        let n = f.len();
        let mut u = vec![0.0; n];
        let mut r = f.to_vec();
        let mut z = pre(&r);
        let mut p = z.clone();
        let mut rz = dot(&r, &z);
        let fnorm = dot(f, f).sqrt();
        let mut q = vec![0.0; n];
        for it in 0..5000 {
            if dot(&r, &r).sqrt() <= 1e-8 * fnorm {
                return it;
            }
            ap(&p, &mut q);
            let alpha = rz / dot(&p, &q);
            for i in 0..n {
                u[i] += alpha * p[i];
                r[i] -= alpha * q[i];
            }
            z = pre(&r);
            let rz_new = dot(&r, &z);
            let beta = rz_new / rz;
            rz = rz_new;
            for i in 0..n {
                p[i] = z[i] + beta * p[i];
            }
        }
        5000
    }

    /// Multigrid must cut the CG steps well below Jacobi's.
    #[test]
    fn multigrid_beats_jacobi() {
        let b = block([16, 12, 10]);
        let f = pseudo(b.fixed.len(), 5, &b.fixed);
        let settings = Settings { coarsest_dofs: 300, ..Settings::default() };
        let (inv, ap) = (b.inv_diag(), |u: &[f64], y: &mut [f64]| b.apply(u, y));
        let h = Hierarchy3::build(&b.system(), settings, &FineOps { apply: &ap, inv: &inv }).unwrap();
        let (plain, mg) = (cg_steps(&b, None, &f), cg_steps(&b, Some(&h), &f));
        eprintln!("levels {}, jacobi {plain} steps, multigrid {mg} steps", h.depth());
        assert!(mg * 3 < plain, "jacobi {plain}, multigrid {mg}");
    }

}
