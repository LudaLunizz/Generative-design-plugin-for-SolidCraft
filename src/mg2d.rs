//! Stiffness operators on a regular 2D grid and a geometric multigrid preconditioner for them.
//!
//! A level is a grid of `nx` x `ny` square cells, each with its own 8x8 stiffness matrix (the fine
//! level stores one shared matrix and a modulus per cell). DOFs are numbered node by node,
//! column-major (node (i, j) has DOFs 2(i(ny+1)+j) and +1). Fixed DOFs are a mask: an operator
//! returns zero there, so vectors stay zero on supports.
//!
//! Coarse levels halve each direction. Their cell matrices are the Galerkin product P^T K P
//! computed per coarse cell: on nested bilinear grids, interpolation inside a coarse cell only
//! uses that cell's four corners, so assembling `P_e^T K_patch P_e` over coarse cells is exactly
//! the global `P^T K P`. That keeps strong/void contrasts (SIMP's 1e9) on the coarse grids, which
//! averaging densities would blur. The V-cycle uses damped Jacobi smoothing (the same sweeps
//! before and after, so the preconditioner is symmetric) and a banded Cholesky solve on the
//! coarsest grid.

/// The DOFs of cell (ex, ey) on a grid `ny` cells high.
#[inline]
pub fn cell_dofs(ny: usize, ex: usize, ey: usize) -> [usize; 8] {
    let n = ny + 1;
    let n0 = ex * n + ey;
    let n1 = n0 + n;
    [2 * n0, 2 * n0 + 1, 2 * n1, 2 * n1 + 1, 2 * n1 + 2, 2 * n1 + 3, 2 * n0 + 2, 2 * n0 + 3]
}

pub type Mat8 = [[f64; 8]; 8];

pub enum Cells {
    /// Fine level: one matrix times a modulus per cell.
    Scaled { ke: Box<Mat8>, modulus: Vec<f64> },
    /// Coarse levels: a matrix per cell.
    Full(Vec<Mat8>),
}

pub struct Level {
    pub nx: usize,
    pub ny: usize,
    pub cells: Cells,
    pub fixed: Vec<bool>,
}

impl Level {
    pub fn ndof(&self) -> usize {
        2 * (self.nx + 1) * (self.ny + 1)
    }

    fn cell_matrix(&self, e: usize) -> Mat8 {
        match &self.cells {
            Cells::Scaled { ke, modulus } => {
                let m = modulus.get(e).copied().unwrap_or(0.0);
                ke.map(|row| row.map(|v| v * m))
            }
            Cells::Full(all) => all.get(e).copied().unwrap_or([[0.0; 8]; 8]),
        }
    }

    /// y = K u (zero on fixed DOFs). Walks the grid one column of cells at a time, reading the
    /// two node columns each cell column touches as slices, so the inner loop does no lookups.
    pub fn apply(&self, u: &[f64], y: &mut [f64]) {
        y.iter_mut().for_each(|v| *v = 0.0);
        let col = 2 * (self.ny + 1);
        for ex in 0..self.nx {
            let (left, right) = y.split_at_mut((ex + 1) * col);
            let (Some(y0), Some(y1), Some(u0), Some(u1)) =
                (left.get_mut(ex * col..), right.get_mut(..col), u.get(ex * col..(ex + 1) * col), u.get((ex + 1) * col..(ex + 2) * col))
            else {
                return;
            };
            let cells = ex * self.ny..(ex + 1) * self.ny;
            match &self.cells {
                Cells::Scaled { ke, modulus } => {
                    let Some(m) = modulus.get(cells) else { return };
                    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
                    {
                        let packed = simd::pack(ke);
                        column_product(m.iter().map(|m| (&packed, *m)), simd::cell_product, u0, u1, y0, y1);
                    }
                    #[cfg(not(all(target_arch = "wasm32", target_feature = "simd128")))]
                    column_product(m.iter().map(|m| (&**ke, *m)), cell_product, u0, u1, y0, y1);
                }
                Cells::Full(all) => {
                    let Some(k) = all.get(cells) else { return };
                    column_product(k.iter().map(|k| (k, 1.0)), cell_product, u0, u1, y0, y1);
                }
            }
        }
        for (out, fixed) in y.iter_mut().zip(&self.fixed) {
            if *fixed {
                *out = 0.0;
            }
        }
    }

    /// The diagonal, with 1 on fixed DOFs and wherever it is not positive.
    pub fn diagonal(&self) -> Vec<f64> {
        let mut diag = vec![0.0; self.ndof()];
        let mut e = 0;
        for ex in 0..self.nx {
            for ey in 0..self.ny {
                let k = self.cell_matrix(e);
                for (i, d) in cell_dofs(self.ny, ex, ey).iter().enumerate() {
                    if let (Some(out), Some(row)) = (diag.get_mut(*d), k.get(i)) {
                        *out += row.get(i).copied().unwrap_or(0.0);
                    }
                }
                e += 1;
            }
        }
        for (v, fixed) in diag.iter_mut().zip(&self.fixed) {
            if *fixed || !(v.is_finite() && *v > 0.0) {
                *v = 1.0;
            }
        }
        diag
    }

    /// The next coarser level, or None when this grid cannot be halved.
    pub fn coarsen(&self) -> Option<Level> {
        if self.nx % 2 != 0 || self.ny % 2 != 0 || self.nx < 4 || self.ny < 4 {
            return None;
        }
        let (cx, cy) = (self.nx / 2, self.ny / 2);
        // A coarse DOF is fixed where its node coincides with a fixed fine node.
        let mut fixed = vec![false; 2 * (cx + 1) * (cy + 1)];
        for i in 0..=cx {
            for j in 0..=cy {
                for c in 0..2 {
                    let f = 2 * ((2 * i) * (self.ny + 1) + 2 * j) + c;
                    if self.fixed.get(f).copied().unwrap_or(false)
                        && let Some(slot) = fixed.get_mut(2 * (i * (cy + 1) + j) + c)
                    {
                        *slot = true;
                    }
                }
            }
        }
        // Per coarse cell: kc = sum over its 4 fine cells of P_f^T K_f P_f, with P_f the 8 x 8
        // interpolation from the coarse cell's DOFs to that fine cell's DOFs. Each fine node takes
        // its value from 1, 2 or 4 coarse corners, so P_f has at most 18 non-zeros; rows of fixed
        // fine DOFs and columns of fixed coarse DOFs are left out, which drops the supports.
        // Coarse corner order matches cell_dofs: (0,0), (2,0), (2,2), (0,2) in fine-node units.
        let corners = [(0usize, 0usize), (2, 0), (2, 2), (0, 2)];
        let mut cells = Vec::with_capacity(cx * cy);
        let mut entries: Vec<(usize, usize, f64)> = Vec::with_capacity(18);
        for ex in 0..cx {
            for ey in 0..cy {
                let coarse = cell_dofs(cy, ex, ey);
                let mut kc = [[0.0f64; 8]; 8];
                for a in 0..2 {
                    for b in 0..2 {
                        let (fx, fy) = (2 * ex + a, 2 * ey + b);
                        let k = self.cell_matrix(fx * self.ny + fy);
                        let fine = cell_dofs(self.ny, fx, fy);
                        // The fine cell's nodes in the coarse cell's 3x3 node patch.
                        let nodes = [(a, b), (a + 1, b), (a + 1, b + 1), (a, b + 1)];
                        entries.clear();
                        for (node, (p, q)) in nodes.iter().enumerate() {
                            for (corner, (px, qy)) in corners.iter().enumerate() {
                                let w = (1.0 - (*p as f64 - *px as f64).abs() / 2.0).max(0.0) * (1.0 - (*q as f64 - *qy as f64).abs() / 2.0).max(0.0);
                                if w == 0.0 {
                                    continue;
                                }
                                for c in 0..2 {
                                    let (fi, ci) = (2 * node + c, 2 * corner + c);
                                    let fine_fixed = fine.get(fi).is_none_or(|d| self.fixed.get(*d).copied().unwrap_or(true));
                                    let coarse_fixed = coarse.get(ci).is_none_or(|d| fixed.get(*d).copied().unwrap_or(true));
                                    if !fine_fixed && !coarse_fixed {
                                        entries.push((fi, ci, w));
                                    }
                                }
                            }
                        }
                        // kp = K_f P_f (8 x 8), then kc += P_f^T kp.
                        let mut kp = [[0.0f64; 8]; 8];
                        for &(fi, ci, w) in &entries {
                            for (kprow, krow) in kp.iter_mut().zip(&k) {
                                if let (Some(out), Some(kv)) = (kprow.get_mut(ci), krow.get(fi)) {
                                    *out += kv * w;
                                }
                            }
                        }
                        for &(fi, ci, w) in &entries {
                            if let (Some(kcrow), Some(kprow)) = (kc.get_mut(ci), kp.get(fi)) {
                                for (out, v) in kcrow.iter_mut().zip(kprow) {
                                    *out += w * v;
                                }
                            }
                        }
                    }
                }
                cells.push(kc);
            }
        }
        Some(Level { nx: cx, ny: cy, cells: Cells::Full(cells), fixed })
    }

    /// Bilinear interpolation of a coarse vector, added to `fine` (this level is the fine one).
    fn prolong_add(&self, coarse: &[f64], fine: &mut [f64]) {
        let cy = self.ny / 2;
        for i in 0..=self.nx {
            for j in 0..=self.ny {
                for c in 0..2 {
                    let f = 2 * (i * (self.ny + 1) + j) + c;
                    if self.fixed.get(f).copied().unwrap_or(true) {
                        continue;
                    }
                    let mut v = 0.0;
                    for (ci, wi) in parents(i) {
                        for (cj, wj) in parents(j) {
                            v += wi * wj * coarse.get(2 * (ci * (cy + 1) + cj) + c).copied().unwrap_or(0.0);
                        }
                    }
                    if let Some(slot) = fine.get_mut(f) {
                        *slot += v;
                    }
                }
            }
        }
    }

    /// The transpose of `prolong_add`: fine residual to the coarse grid.
    fn restrict(&self, fine: &[f64], coarse_fixed: &[bool]) -> Vec<f64> {
        let cy = self.ny / 2;
        let mut out = vec![0.0; coarse_fixed.len()];
        for i in 0..=self.nx {
            for j in 0..=self.ny {
                for c in 0..2 {
                    let f = 2 * (i * (self.ny + 1) + j) + c;
                    if self.fixed.get(f).copied().unwrap_or(true) {
                        continue;
                    }
                    let v = fine.get(f).copied().unwrap_or(0.0);
                    for (ci, wi) in parents(i) {
                        for (cj, wj) in parents(j) {
                            if let Some(slot) = out.get_mut(2 * (ci * (cy + 1) + cj) + c) {
                                *slot += wi * wj * v;
                            }
                        }
                    }
                }
            }
        }
        for (v, fixed) in out.iter_mut().zip(coarse_fixed) {
            if *fixed {
                *v = 0.0;
            }
        }
        out
    }
}

/// m * K u for one cell. K is symmetric, so K u = sum_j u_j * (row j): eight independent running
/// sums instead of eight dependent dot-product chains.
#[inline]
fn cell_product(k: &Mat8, ue: [f64; 8], m: f64) -> [f64; 8] {
    let mut ye = [0.0f64; 8];
    for (row, v) in k.iter().zip(ue) {
        for (y, kv) in ye.iter_mut().zip(row) {
            *y += kv * v;
        }
    }
    ye.map(|v| v * m)
}

/// The same product with WebAssembly SIMD: two f64 per operation. The fuel SolidCraft meters
/// counts operations, so this roughly halves the cost of the hottest loop.
#[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
mod simd {
    use core::arch::wasm32::{f64x2, f64x2_add, f64x2_extract_lane, f64x2_mul, f64x2_splat, v128};

    pub type Packed = [[v128; 4]; 8];

    pub fn pack(k: &super::Mat8) -> Packed {
        k.map(|[a, b, c, d, e, f, g, h]| [f64x2(a, b), f64x2(c, d), f64x2(e, f), f64x2(g, h)])
    }

    #[inline]
    pub fn cell_product(k: &Packed, ue: [f64; 8], m: f64) -> [f64; 8] {
        let mut acc = [f64x2_splat(0.0); 4];
        for (row, v) in k.iter().zip(ue) {
            let s = f64x2_splat(v * m);
            for (a, kv) in acc.iter_mut().zip(row) {
                *a = f64x2_add(*a, f64x2_mul(*kv, s));
            }
        }
        let [a, b, c, d] = acc;
        [
            f64x2_extract_lane::<0>(a),
            f64x2_extract_lane::<1>(a),
            f64x2_extract_lane::<0>(b),
            f64x2_extract_lane::<1>(b),
            f64x2_extract_lane::<0>(c),
            f64x2_extract_lane::<1>(c),
            f64x2_extract_lane::<0>(d),
            f64x2_extract_lane::<1>(d),
        ]
    }
}

/// y += K u for one column of cells: `u0`/`y0` are the left node column, `u1`/`y1` the right one
/// (x, y per node, bottom to top). Cell `ey` couples nodes ey and ey+1 of both columns; the
/// upper node's share is carried to the next cell so every output is written once, in order.
#[inline]
fn column_product<'a, T: 'a>(
    cells: impl Iterator<Item = (&'a T, f64)>,
    product: impl Fn(&T, [f64; 8], f64) -> [f64; 8],
    u0: &[f64],
    u1: &[f64],
    y0: &mut [f64],
    y1: &mut [f64],
) {
    let mut out0 = y0.chunks_exact_mut(2);
    let mut out1 = y1.chunks_exact_mut(2);
    let mut carry0 = [0.0f64; 2];
    let mut carry1 = [0.0f64; 2];
    for ((w0, w1), (k, m)) in u0.windows(4).step_by(2).zip(u1.windows(4).step_by(2)).zip(cells) {
        let (Ok([a0, a1, a2, a3]), Ok([b0, b1, b2, b3])) = (<[f64; 4]>::try_from(w0), <[f64; 4]>::try_from(w1)) else { break };
        // Cell DOF order: n0 = left lower, n1 = right lower, n2 = right upper, n3 = left upper.
        let ue = [a0, a1, b0, b1, b2, b3, a2, a3];
        let [e0, e1, e2, e3, e4, e5, e6, e7] = product(k, ue, m);
        if let Some([p, q]) = out0.next() {
            *p += e0 + carry0[0];
            *q += e1 + carry0[1];
        }
        if let Some([p, q]) = out1.next() {
            *p += e2 + carry1[0];
            *q += e3 + carry1[1];
        }
        carry0 = [e6, e7];
        carry1 = [e4, e5];
    }
    if let Some([p, q]) = out0.next() {
        *p += carry0[0];
        *q += carry0[1];
    }
    if let Some([p, q]) = out1.next() {
        *p += carry1[0];
        *q += carry1[1];
    }
}

/// The coarse nodes (index, weight) a fine node index interpolates from.
fn parents(i: usize) -> impl Iterator<Item = (usize, f64)> {
    let (a, b) = if i % 2 == 0 { ((i / 2, 1.0), None) } else { ((i / 2, 0.5), Some((i / 2 + 1, 0.5))) };
    std::iter::once(a).chain(b)
}

/// A banded Cholesky factor of a level's assembled matrix (fixed DOFs as identity rows).
pub struct Banded {
    n: usize,
    b: usize,
    l: Vec<f64>, // row i, offset k = i - j, at i * (b + 1) + k
}

impl Banded {
    pub fn factor(level: &Level) -> Banded {
        let n = level.ndof();
        let b = 2 * (level.ny + 1) + 3;
        let mut a = vec![0.0; n * (b + 1)];
        let mut e = 0;
        for ex in 0..level.nx {
            for ey in 0..level.ny {
                let k = level.cell_matrix(e);
                let dofs = cell_dofs(level.ny, ex, ey);
                for (row, i) in k.iter().zip(dofs) {
                    for (v, j) in row.iter().zip(dofs) {
                        let fixed = level.fixed.get(i).copied().unwrap_or(true) || level.fixed.get(j).copied().unwrap_or(true);
                        if j <= i && !fixed && i - j <= b
                            && let Some(slot) = a.get_mut(i * (b + 1) + (i - j))
                        {
                            *slot += v;
                        }
                    }
                }
                e += 1;
            }
        }
        for (i, fixed) in level.fixed.iter().enumerate() {
            if *fixed && let Some(slot) = a.get_mut(i * (b + 1)) {
                *slot = 1.0;
            }
        }
        // In-place LL^T within the band.
        let at = |v: &Vec<f64>, i: usize, j: usize| if i >= j && i - j <= b { v.get(i * (b + 1) + (i - j)).copied().unwrap_or(0.0) } else { 0.0 };
        for i in 0..n {
            let lo = i.saturating_sub(b);
            for j in lo..=i {
                let mut s = at(&a, i, j);
                for m in lo.max(j.saturating_sub(b))..j {
                    s -= at(&a, i, m) * at(&a, j, m);
                }
                let value = if i == j {
                    if s.is_finite() && s > 0.0 { s.sqrt() } else { 1.0 }
                } else {
                    let d = at(&a, j, j);
                    if d != 0.0 { s / d } else { 0.0 }
                };
                if let Some(slot) = a.get_mut(i * (b + 1) + (i - j)) {
                    *slot = value;
                }
            }
        }
        Banded { n, b, l: a }
    }

    fn get(&self, i: usize, j: usize) -> f64 {
        if i >= j && i - j <= self.b { self.l.get(i * (self.b + 1) + (i - j)).copied().unwrap_or(0.0) } else { 0.0 }
    }

    fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        let mut y = rhs.to_vec();
        for i in 0..self.n {
            let mut s = y.get(i).copied().unwrap_or(0.0);
            for j in i.saturating_sub(self.b)..i {
                s -= self.get(i, j) * y.get(j).copied().unwrap_or(0.0);
            }
            let d = self.get(i, i);
            if let Some(slot) = y.get_mut(i) {
                *slot = if d != 0.0 { s / d } else { 0.0 };
            }
        }
        for i in (0..self.n).rev() {
            let mut s = y.get(i).copied().unwrap_or(0.0);
            for j in i + 1..(i + self.b + 1).min(self.n) {
                s -= self.get(j, i) * y.get(j).copied().unwrap_or(0.0);
            }
            let d = self.get(i, i);
            if let Some(slot) = y.get_mut(i) {
                *slot = if d != 0.0 { s / d } else { 0.0 };
            }
        }
        y
    }
}

/// Settings of the V-cycle.
#[derive(Clone, Copy, Debug)]
pub struct Settings {
    /// Jacobi sweeps before and after the coarse correction.
    pub sweeps: usize,
    /// Jacobi damping.
    pub omega: f64,
    /// Coarsen until the banded Cholesky of the coarsest grid costs at most this many
    /// multiply-adds (n * b^2).
    pub coarse_budget: usize,
    pub max_levels: usize,
}

impl Default for Settings {
    fn default() -> Self {
        Self { sweeps: 2, omega: 0.6, coarse_budget: 2_000_000, max_levels: 10 }
    }
}

enum Coarsest {
    Direct(Box<Banded>),
    /// The grid could not be halved far enough for a cheap direct solve: smooth instead.
    Sweeps(usize),
}

/// The multigrid preconditioner: the levels, their inverse diagonals and the coarsest solver.
pub struct Hierarchy {
    pub levels: Vec<Level>,
    inv_diag: Vec<Vec<f64>>,
    coarsest: Coarsest,
    settings: Settings,
}

impl Hierarchy {
    pub fn new(fine: Level, settings: Settings) -> Hierarchy {
        let mut levels = vec![fine];
        while let Some(last) = levels.last() {
            let cost = last.ndof() * (2 * (last.ny + 1) + 3).pow(2);
            if cost <= settings.coarse_budget || levels.len() >= settings.max_levels {
                break;
            }
            match last.coarsen() {
                Some(c) => levels.push(c),
                None => break,
            }
        }
        let inv_diag = levels.iter().map(|l| l.diagonal().iter().map(|d| 1.0 / d).collect()).collect();
        let coarsest = match levels.last() {
            Some(l) if l.ndof() * (2 * (l.ny + 1) + 3).pow(2) <= settings.coarse_budget.saturating_mul(4) => Coarsest::Direct(Box::new(Banded::factor(l))),
            _ => Coarsest::Sweeps(30),
        };
        Hierarchy { levels, inv_diag, coarsest, settings }
    }

    pub fn depth(&self) -> usize {
        self.levels.len()
    }

    pub fn coarsest_is_direct(&self) -> bool {
        matches!(self.coarsest, Coarsest::Direct(_))
    }

    pub fn fine(&self) -> Option<&Level> {
        self.levels.first()
    }

    /// z = M^-1 r: one V-cycle from zero.
    pub fn precondition(&self, r: &[f64]) -> Vec<f64> {
        self.vcycle(0, r)
    }

    fn smooth(&self, l: usize, r: &[f64], z: &mut [f64], sweeps: usize, scratch: &mut [f64]) {
        let (Some(level), Some(inv)) = (self.levels.get(l), self.inv_diag.get(l)) else { return };
        for _ in 0..sweeps {
            level.apply(z, scratch);
            for ((zi, (ri, ki)), di) in z.iter_mut().zip(r.iter().zip(scratch.iter())).zip(inv) {
                *zi += self.settings.omega * di * (ri - ki);
            }
        }
    }

    fn vcycle(&self, l: usize, r: &[f64]) -> Vec<f64> {
        let n = r.len();
        let mut z = vec![0.0; n];
        let mut scratch = vec![0.0; n];
        if l + 1 >= self.levels.len() {
            match &self.coarsest {
                Coarsest::Direct(chol) => return chol.solve(r),
                Coarsest::Sweeps(s) => {
                    self.smooth(l, r, &mut z, *s, &mut scratch);
                    return z;
                }
            }
        }
        let (Some(level), Some(coarse)) = (self.levels.get(l), self.levels.get(l + 1)) else { return z };
        self.smooth(l, r, &mut z, self.settings.sweeps, &mut scratch);
        level.apply(&z, &mut scratch);
        let residual: Vec<f64> = r.iter().zip(&scratch).map(|(a, b)| a - b).collect();
        let rc = level.restrict(&residual, &coarse.fixed);
        let zc = self.vcycle(l + 1, &rc);
        level.prolong_add(&zc, &mut z);
        self.smooth(l, r, &mut z, self.settings.sweeps, &mut scratch);
        z
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
    use super::*;
    use crate::topo2d::element_stiffness;

    fn fine(nx: usize, ny: usize, modulus: impl Fn(usize, usize) -> f64) -> Level {
        let mut fixed = vec![false; 2 * (nx + 1) * (ny + 1)];
        for j in 0..=ny {
            fixed[2 * j] = true;
            fixed[2 * j + 1] = true;
        }
        let mut m = Vec::new();
        for ex in 0..nx {
            for ey in 0..ny {
                m.push(modulus(ex, ey));
            }
        }
        Level { nx, ny, cells: Cells::Scaled { ke: Box::new(element_stiffness(0.3)), modulus: m }, fixed }
    }

    fn pseudo_random(n: usize, seed: u64) -> Vec<f64> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                ((s >> 33) as f64 / (1u64 << 31) as f64) - 0.5
            })
            .collect()
    }

    /// The coarse cell matrices must equal P^T K P computed the slow way: K applied to the
    /// prolongation of a coarse vector, restricted back.
    #[test]
    fn coarse_operator_is_the_galerkin_product() {
        let f = fine(8, 6, |ex, ey| if (ex + ey) % 3 == 0 { 1e-9 } else { 1.0 + ex as f64 * 0.1 });
        let c = f.coarsen().unwrap();
        let mut v = pseudo_random(c.ndof(), 7);
        for (x, fixed) in v.iter_mut().zip(&c.fixed) {
            if *fixed {
                *x = 0.0;
            }
        }
        let mut pv = vec![0.0; f.ndof()];
        f.prolong_add(&v, &mut pv);
        let mut kpv = vec![0.0; f.ndof()];
        f.apply(&pv, &mut kpv);
        let slow = f.restrict(&kpv, &c.fixed);
        let mut fast = vec![0.0; c.ndof()];
        c.apply(&v, &mut fast);
        for (a, b) in slow.iter().zip(&fast) {
            assert!((a - b).abs() < 1e-10 * (1.0 + a.abs()), "{a} vs {b}");
        }
    }

    #[test]
    fn coarse_matrices_are_symmetric() {
        let f = fine(8, 8, |_, _| 1.0);
        let c = f.coarsen().unwrap();
        if let Cells::Full(all) = &c.cells {
            for k in all {
                for (i, row) in k.iter().enumerate() {
                    for (j, v) in row.iter().enumerate() {
                        assert!((v - k[j][i]).abs() < 1e-12);
                    }
                }
            }
        } else {
            panic!("coarse level should hold full matrices");
        }
    }

    #[test]
    fn banded_cholesky_solves_the_system() {
        let l = fine(6, 4, |ex, _| 1.0 + ex as f64);
        let chol = Banded::factor(&l);
        let mut x = pseudo_random(l.ndof(), 3);
        for (v, fixed) in x.iter_mut().zip(&l.fixed) {
            if *fixed {
                *v = 0.0;
            }
        }
        let mut b = vec![0.0; l.ndof()];
        l.apply(&x, &mut b);
        let y = chol.solve(&b);
        for (a, e) in y.iter().zip(&x) {
            assert!((a - e).abs() < 1e-9, "{a} vs {e}");
        }
    }

    #[test]
    fn odd_grids_stop_coarsening() {
        let h = Hierarchy::new(fine(30, 15, |_, _| 1.0), Settings { coarse_budget: 0, ..Default::default() });
        assert_eq!(h.depth(), 1);
        let h = Hierarchy::new(fine(64, 32, |_, _| 1.0), Settings { coarse_budget: 0, ..Default::default() });
        assert!(h.depth() >= 4);
    }
}
