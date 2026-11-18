//! A long beam, fixed at one end and loaded at the other, designed by 2D topology optimisation of
//! its side view and then extruded to a plate of constant thickness.
//!
//! Physical units throughout: mm, N, MPa, kg. The optimiser works on a unit problem (`topo2d`);
//! plane stress makes the displacement of a cell grid independent of the cell size, so the unit
//! result scales to physical units as `F / (E t)`. The lightest design that meets the stated
//! limits is found by bisection on the volume fraction, for each candidate plate thickness; the
//! lateral-torsional stability of the plate is estimated from the design's cross-sections.
//!
//! Limits, stated plainly: the loads are in the plane of the plate only; self-weight is checked
//! on the finished design, not designed for; stress at stair-stepped corners is read as a high
//! percentile (the bilinear grid has artificial peaks there); member buckling is not checked.
#![allow(clippy::indexing_slicing)] // grid indices are bounded by nx, ny by construction

use crate::mg2d::cell_dofs;
use crate::topo2d::{self, LoadCase, Problem, Solver, System};

/// Multigrid coarsest-level budget for beam grids: the slender grid halves down to about 125 x 20
/// cells, where a banded direct solve (about 11 million multiply-adds) beats 30 Jacobi sweeps.
const COARSE_BUDGET: usize = 12_000_000;

#[derive(Clone, Debug)]
pub struct BeamSpec {
    pub length: f64,
    /// Depth at the fixed end (the grid's full height), mm.
    pub depth: f64,
    /// Depth at the free end, mm; the design area tapers in a straight line (>= depth: no taper).
    pub tip_depth: f64,
    /// Thinnest member allowed, mm: thinner features are removed by a morphological opening.
    pub min_member_mm: f64,
    pub elem: f64,
    pub e_mpa: f64,
    pub nu: f64,
    pub density_g_cm3: f64,
    pub yield_mpa: f64,
    pub load_n: f64,
    pub pad_mm: f64,
    pub safety: f64,
    pub buckling_safety: f64,
    pub max_defl_mm: f64,
    pub rmin_cells: f64,
    pub iterations: usize,
    pub self_weight: bool,
    /// Assumed unbraced length of a compression member for the Euler check, mm.
    pub member_len_mm: f64,
    /// Solve on the GPU (builds with the `gpu` feature; otherwise the CPU multigrid).
    pub gpu: bool,
}

impl BeamSpec {
    pub fn nx(&self) -> usize {
        (self.length / self.elem).round() as usize
    }
    pub fn ny(&self) -> usize {
        (self.depth / self.elem).round() as usize
    }
    /// The minimum member width in cells (the opening window), at least 1.
    pub fn min_member_cells(&self) -> usize {
        ((self.min_member_mm / self.elem).round() as usize).max(1)
    }
    /// Rows open at the free end (0 = no taper), as `Problem::tip_rows`.
    pub fn tip_rows(&self) -> usize {
        if self.tip_depth >= self.depth { 0 } else { ((self.tip_depth / self.elem).round() as usize).max(3) }
    }
    /// The top row of material at the free end plus one: the row the load pad sits under.
    pub fn top_row(&self) -> usize {
        match self.tip_rows() {
            0 => self.ny(),
            t => t,
        }
    }
    fn pad_cells(&self) -> usize {
        ((self.pad_mm / self.elem).round() as usize).max(1)
    }
    pub fn check(&self) -> Result<(), String> {
        for (v, n) in [(self.length, "length"), (self.depth, "depth"), (self.tip_depth, "tip_depth"), (self.elem, "elem"), (self.e_mpa, "e_mpa"), (self.yield_mpa, "yield_mpa"), (self.load_n, "load_n"), (self.max_defl_mm, "max_defl_mm")] {
            if !(v.is_finite() && v > 0.0) {
                return Err(format!("{n} must be a positive number"));
            }
        }
        self.problem(0.3, 0.0).map(|_| ())
    }
    fn problem(&self, volfrac: f64, change_tol: f64) -> Result<Problem, String> {
        let p = Problem {
            nx: self.nx(),
            ny: self.ny(),
            volfrac,
            penal: 3.0,
            rmin: self.rmin_cells,
            nu: self.nu,
            load: LoadCase::TipPad(self.pad_cells()),
            max_iterations: self.iterations,
            change_tol,
            cg_tol: 1e-4,
            cg_max: 5000,
            solver: if self.gpu { Solver::Gpu } else { Solver::Multigrid },
            beta_final: 16.0,
            beta_every: (self.iterations / 8).max(5),
            tip_rows: self.tip_rows(),
            coarse_budget: COARSE_BUDGET,
            pad_rows: self.min_member_cells().max(2),
        };
        p.check()?;
        Ok(p)
    }
}

/// A finished, cleaned binary design (column-major like the optimiser: index = ex * ny + ey,
/// ey = 0 at the bottom) and what the finite element check says about it.
#[derive(Clone, Debug)]
pub struct Design {
    pub nx: usize,
    pub ny: usize,
    pub solid: Vec<bool>,
    pub volfrac: f64,
    pub thickness: f64,
    pub mass_kg: f64,
    pub tip_defl_mm: f64,
    pub tip_defl_no_weight_mm: f64,
    pub vm_max: f64,
    pub vm_p99: f64,
    pub weight_n: f64,
    pub ltb_pcr_n: f64,
    pub connected: bool,
    /// Width of the thinnest members (10th percentile over the design's medial ridge), mm.
    pub member_width_mm: f64,
    /// Euler buckling stress of such a member over the p99 stress (>= safety needed).
    pub member_margin: f64,
}

impl Design {
    pub fn ok(&self, s: &BeamSpec) -> bool {
        self.connected && self.tip_defl_mm <= s.max_defl_mm && self.vm_p99 * s.safety <= s.yield_mpa && self.member_margin >= s.safety
    }
}

/// Squared distance in cells from each solid cell to the nearest empty cell (exact Euclidean
/// transform, Felzenszwalb-Huttenlocher); 0 on empty cells. Cells outside the grid are ignored.
fn edt_sq(solid: &[bool], nx: usize, ny: usize) -> Vec<f64> {
    const INF: f64 = 1e18;
    let mut f: Vec<f64> = solid.iter().map(|&b| if b { INF } else { 0.0 }).collect();
    let line = |v: &[f64]| -> Vec<f64> {
        let n = v.len();
        let mut d = vec![0.0; n];
        let mut k = 0usize;
        let mut z = vec![0.0; n + 1];
        let mut vtx = vec![0usize; n];
        z[0] = -INF;
        z[1] = INF;
        for q in 1..n {
            loop {
                let p = vtx[k];
                let s = ((v[q] + (q * q) as f64) - (v[p] + (p * p) as f64)) / (2.0 * q as f64 - 2.0 * p as f64);
                if s <= z[k] && k > 0 {
                    k -= 1;
                    continue;
                }
                if s <= z[k] {
                    vtx[0] = q;
                    z[0] = -INF;
                    z[1] = INF;
                } else {
                    k += 1;
                    vtx[k] = q;
                    z[k] = s;
                    z[k + 1] = INF;
                }
                break;
            }
        }
        k = 0;
        for (q, out) in d.iter_mut().enumerate() {
            while z[k + 1] < q as f64 {
                k += 1;
            }
            let p = vtx[k];
            *out = ((q as f64 - p as f64).powi(2)) + v[p];
        }
        d
    };
    for i in 0..nx {
        let col = line(&f[i * ny..(i + 1) * ny]);
        f[i * ny..(i + 1) * ny].copy_from_slice(&col);
    }
    for j in 0..ny {
        let row: Vec<f64> = (0..nx).map(|i| f[i * ny + j]).collect();
        let r = line(&row);
        for (i, v) in r.iter().enumerate() {
            f[i * ny + j] = *v;
        }
    }
    f
}

/// Distance in cells from each cell to the nearest empty cell; cells outside the top, bottom and
/// right edges count as empty. Returned as the half width: distance minus half a cell.
fn distance_to_void(solid: &[bool], nx: usize, ny: usize) -> Vec<f64> {
    const INF: f64 = 1e18;
    let f = edt_sq(solid, nx, ny);
    let mut out = vec![0.0; nx * ny];
    for i in 0..nx {
        for j in 0..ny {
            let d = f[i * ny + j].sqrt().min(INF.sqrt());
            let edge = (j + 1).min(ny - j).min(nx - i) as f64 - 0.5; // half width to the outer edges
            out[i * ny + j] = if solid[i * ny + j] { (d - 0.5).min(edge).max(0.0) } else { 0.0 };
        }
    }
    out
}

/// Width of the thinnest members in mm: 2 x the half width on the design's medial ridge,
/// 10th percentile.
pub fn member_width(solid: &[bool], nx: usize, ny: usize, h: f64) -> f64 {
    let half = distance_to_void(solid, nx, ny);
    let mut widths = Vec::new();
    for i in 1..nx.saturating_sub(1) {
        for j in 1..ny.saturating_sub(1) {
            let c = half[i * ny + j];
            if c <= 0.0 {
                continue;
            }
            let ridge = (-1i64..=1).all(|di| (-1i64..=1).all(|dj| half[((i as i64 + di) as usize) * ny + (j as i64 + dj) as usize] <= c + 1e-9));
            if ridge {
                widths.push(2.0 * c * h);
            }
        }
    }
    widths.sort_by(|a, b| a.total_cmp(b));
    widths.get(widths.len() / 10).copied().unwrap_or(0.0)
}

/// Keeps the cells connected to the fixed edge, closes diagonal-only contacts and small holes so
/// the outline is made of simple loops. Returns the cleaned design and whether the free end's top
/// is still attached.
pub fn clean(density: &[f64], nx: usize, ny: usize, pad: usize, top: usize) -> (Vec<bool>, bool) {
    let mut solid: Vec<bool> = density.iter().map(|d| *d >= 0.5).collect();
    // Diagonal-only contacts (checkerboard 2x2): add one of the two empty cells.
    loop {
        let mut changed = false;
        for i in 0..nx.saturating_sub(1) {
            for j in 0..ny.saturating_sub(1) {
                let (a, b, c, d) = (solid[i * ny + j], solid[(i + 1) * ny + j], solid[i * ny + j + 1], solid[(i + 1) * ny + j + 1]);
                if a && d && !b && !c {
                    solid[(i + 1) * ny + j] = true;
                    changed = true;
                } else if b && c && !a && !d {
                    solid[i * ny + j] = true;
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    // Only material connected (4-neighbours) to the fixed edge stays.
    let mut keep = vec![false; nx * ny];
    let mut stack: Vec<(usize, usize)> = (0..ny).filter(|&j| solid[j]).map(|j| (0, j)).collect();
    for &(i, j) in &stack {
        keep[i * ny + j] = true;
    }
    while let Some((i, j)) = stack.pop() {
        for (di, dj) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
            let (a, b) = (i as i64 + di, j as i64 + dj);
            if a < 0 || b < 0 || a >= nx as i64 || b >= ny as i64 {
                continue;
            }
            let k = a as usize * ny + b as usize;
            if solid[k] && !keep[k] {
                keep[k] = true;
                stack.push((a as usize, b as usize));
            }
        }
    }
    // Holes: empty cells not connected to the outside. Fill those under 12 cells.
    let mut outside = vec![false; nx * ny];
    let mut st: Vec<(usize, usize)> = Vec::new();
    for i in 0..nx {
        for j in [0, ny - 1] {
            if !keep[i * ny + j] && !outside[i * ny + j] {
                outside[i * ny + j] = true;
                st.push((i, j));
            }
        }
    }
    for j in 0..ny {
        let i = nx - 1;
        if !keep[i * ny + j] && !outside[i * ny + j] {
            outside[i * ny + j] = true;
            st.push((i, j));
        }
    }
    while let Some((i, j)) = st.pop() {
        for (di, dj) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
            let (a, b) = (i as i64 + di, j as i64 + dj);
            if a < 0 || b < 0 || a >= nx as i64 || b >= ny as i64 {
                continue;
            }
            let k = a as usize * ny + b as usize;
            if !keep[k] && !outside[k] {
                outside[k] = true;
                st.push((a as usize, b as usize));
            }
        }
    }
    let mut seen = vec![false; nx * ny];
    for s in 0..nx * ny {
        if keep[s] || outside[s] || seen[s] {
            continue;
        }
        let mut comp = vec![s];
        seen[s] = true;
        let mut at = 0;
        while at < comp.len() {
            let (i, j) = (comp[at] / ny, comp[at] % ny);
            at += 1;
            for (di, dj) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
                let (a, b) = (i as i64 + di, j as i64 + dj);
                if a < 0 || b < 0 || a >= nx as i64 || b >= ny as i64 {
                    continue;
                }
                let k = a as usize * ny + b as usize;
                if !keep[k] && !outside[k] && !seen[k] {
                    seen[k] = true;
                    comp.push(k);
                }
            }
        }
        if comp.len() < 12 {
            for k in comp {
                keep[k] = true;
            }
        }
    }
    // The load pad: the top row under the pad must be material for the load to land on the part.
    let top_ok = (nx - pad..nx).all(|i| keep[i * ny + top - 1]);
    (keep, top_ok)
}

/// Finite element check of a cleaned design for the thickness `t`: tip deflection with and
/// without self-weight, stress per cell, mass, lateral-torsional estimate.
pub fn evaluate(s: &BeamSpec, solid: &[bool], top_ok: bool, t: f64) -> Result<Design, String> {
    let (nx, ny, h) = (s.nx(), s.ny(), s.elem);
    let p = s.problem(0.3, 0.0)?;
    let (fixed, unit) = p.supports_and_load();
    let n = ny + 1;
    let modulus: Vec<f64> = solid.iter().map(|&b| if b { s.e_mpa * t } else { s.e_mpa * t * 1e-6 }).collect();
    let rho = s.density_g_cm3 * 1e-6; // kg/mm3
    let cell_weight = rho * h * h * t * 9.81;
    let count = solid.iter().filter(|b| **b).count();
    let run = |load: &[f64]| {
        let sys = System { nx, ny, nu: s.nu, modulus: &modulus, fixed: fixed.clone(), load, solver: if s.gpu { Solver::Gpu } else { Solver::Multigrid }, cg_tol: 1e-8, cg_max: 20_000, coarse_budget: COARSE_BUDGET };
        topo2d::solve(&sys, None, &|| 0.0)
    };
    let tip_load: Vec<f64> = unit.iter().map(|v| v * s.load_n).collect();
    let mut total = tip_load.clone();
    let mut weight = 0.0;
    if s.self_weight {
        for i in 0..nx {
            for j in 0..ny {
                if solid[i * ny + j] {
                    weight += cell_weight;
                    for dof in cell_dofs(ny, i, j).iter().skip(1).step_by(2) {
                        total[*dof] -= cell_weight / 4.0;
                    }
                }
            }
        }
    }
    let tip_dof = |u: &[f64]| -> f64 {
        let pad = s.pad_cells();
        let mut sum = 0.0;
        for i in (nx - pad)..=nx {
            sum += u[2 * (i * n + s.top_row()) + 1];
        }
        -sum / (pad + 1) as f64
    };
    let a0 = run(&tip_load);
    let a1 = if s.self_weight { run(&total) } else { run(&tip_load) };
    if !(a0.converged && a1.converged) {
        return Err("the finite element solve did not converge".into());
    }
    let tip0 = tip_dof(&a0.u);
    let tip1 = tip_dof(&a1.u);
    // Stress at each solid cell's centre (one-point strain), von Mises, plane stress.
    let u = &a1.u;
    let mut vm = Vec::with_capacity(count);
    let c = s.e_mpa / (1.0 - s.nu * s.nu);
    let g = s.e_mpa / (2.0 * (1.0 + s.nu));
    for i in 0..nx {
        for j in 0..ny {
            if !solid[i * ny + j] {
                continue;
            }
            let d = cell_dofs(ny, i, j);
            let (ux, uy): (Vec<f64>, Vec<f64>) = (0..4).map(|k| (u[d[2 * k]], u[d[2 * k + 1]])).unzip();
            let ex = (ux[1] + ux[2] - ux[0] - ux[3]) / (2.0 * h);
            let ey = (uy[2] + uy[3] - uy[0] - uy[1]) / (2.0 * h);
            let gxy = ((ux[2] + ux[3] - ux[0] - ux[1]) + (uy[1] + uy[2] - uy[0] - uy[3])) / (2.0 * h);
            let sx = c * (ex + s.nu * ey);
            let sy = c * (ey + s.nu * ex);
            let tau = g * gxy;
            vm.push((sx * sx - sx * sy + sy * sy + 3.0 * tau * tau).sqrt());
        }
    }
    vm.sort_by(|a, b| a.total_cmp(b));
    let vm_max = vm.last().copied().unwrap_or(0.0);
    let vm_p99 = vm.get(((vm.len() as f64) * 0.99) as usize).copied().unwrap_or(vm_max);
    // Lateral-torsional buckling of a cantilever with the load on the top edge (destabilising):
    // thin-rectangle section properties from the material height at each station, averaged over
    // the root half where the moment is.
    let half = (nx / 2).max(1);
    let mean_height: f64 = (0..half).map(|i| (0..ny).filter(|&j| solid[i * ny + j]).count() as f64 * h).sum::<f64>() / half as f64;
    let (ei, gj) = (s.e_mpa * mean_height * t.powi(3) / 12.0, s.e_mpa / (2.0 * (1.0 + s.nu)) * mean_height * t.powi(3) / 3.0);
    let a = s.tip_depth.min(s.depth) / 2.0;
    let pcr = 4.013 * (ei * gj).sqrt() / (s.length * s.length) * (1.0 - a / s.length * (ei / gj).sqrt());
    let member_width_mm = member_width(solid, nx, ny, h);
    let r = member_width_mm.min(t) / 12f64.sqrt();
    let sigma_cr = std::f64::consts::PI.powi(2) * s.e_mpa * (r / s.member_len_mm).powi(2);
    Ok(Design {
        nx,
        ny,
        solid: solid.to_vec(),
        volfrac: count as f64 / (nx * ny) as f64,
        thickness: t,
        mass_kg: count as f64 * h * h * t * rho,
        tip_defl_mm: tip1,
        tip_defl_no_weight_mm: tip0,
        vm_max,
        vm_p99,
        weight_n: weight,
        ltb_pcr_n: pcr,
        connected: top_ok,
        member_width_mm,
        member_margin: sigma_cr / vm_p99.max(1e-9),
    })
}

/// Morphological opening with a k x k square window: removes every feature narrower than k cells
/// (a cell stays only if some fully solid k x k window covers it). The fixed edge counts as solid
/// beyond the grid, the other edges as empty. Cell-centre distances cannot tell a 1-cell strip
/// from a 2-cell one, so the window is used instead.
pub fn open_design(solid: &[bool], nx: usize, ny: usize, k: usize) -> Vec<bool> {
    if k <= 1 {
        return solid.to_vec();
    }
    let at = |i: i64, j: i64| -> bool {
        if j < 0 || j >= ny as i64 || i >= nx as i64 {
            return false;
        }
        solid[i.max(0) as usize * ny + j as usize]
    };
    let k = k as i64;
    let mut out = vec![false; nx * ny];
    for i in -(k - 1)..nx as i64 {
        for j in -(k - 1)..ny as i64 {
            if (0..k).all(|di| (0..k).all(|dj| at(i + di, j + dj))) {
                for di in 0..k {
                    for dj in 0..k {
                        let (a, b) = (i + di, j + dj);
                        if a >= 0 && b >= 0 && a < nx as i64 && b < ny as i64 {
                            out[a as usize * ny + b as usize] = true;
                        }
                    }
                }
            }
        }
    }
    for (o, s) in out.iter_mut().zip(solid) {
        *o = *o && *s;
    }
    out
}

/// SIMP at one volume fraction, thresholded at the density that keeps the target volume (a fixed
/// 0.5 cuts a grey truss apart) and cleaned. The topology does not depend on the plate thickness:
/// plane stress makes the unit problem scale-free.
pub fn topology(s: &BeamSpec, volfrac: f64) -> Result<(Vec<bool>, bool), String> {
    let p = s.problem(volfrac, 0.005)?;
    let out = topo2d::optimise(&p, |_, _| true);
    let mut sorted = out.density.clone();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let at = (((1.0 - volfrac) * sorted.len() as f64) as usize).min(sorted.len() - 1);
    let thr = sorted[at].clamp(0.05, 0.5);
    let binary: Vec<f64> = out.density.iter().map(|d| if *d >= thr { 1.0 } else { 0.0 }).collect();
    let (solid, _) = clean(&binary, p.nx, p.ny, s.pad_cells(), s.top_row());
    let solid = if s.min_member_mm > 0.0 { open_design(&solid, p.nx, p.ny, s.min_member_cells()) } else { solid };
    let as_density: Vec<f64> = solid.iter().map(|&b| if b { 1.0 } else { 0.0 }).collect();
    Ok(clean(&as_density, p.nx, p.ny, s.pad_cells(), s.top_row()))
}

pub type Loop = Vec<(f64, f64)>;

/// Boundary loops of the design in mm (x along the beam, y up): the outer boundary counter-
/// clockwise first, holes clockwise after it. Edges lie on grid lines; `simplify` then merges
/// staircases into straight runs within `tol` mm.
pub fn outline(d: &Design, h: f64, tol: f64) -> Vec<Loop> {
    let (nx, ny) = (d.nx, d.ny);
    let solid = |i: i64, j: i64| i >= 0 && j >= 0 && i < nx as i64 && j < ny as i64 && d.solid[i as usize * ny + j as usize];
    let mut next: std::collections::HashMap<(i64, i64), (i64, i64)> = std::collections::HashMap::new();
    for i in 0..nx as i64 {
        for j in 0..ny as i64 {
            if !solid(i, j) {
                continue;
            }
            if !solid(i, j - 1) {
                next.insert((i, j), (i + 1, j));
            }
            if !solid(i + 1, j) {
                next.insert((i + 1, j), (i + 1, j + 1));
            }
            if !solid(i, j + 1) {
                next.insert((i + 1, j + 1), (i, j + 1));
            }
            if !solid(i - 1, j) {
                next.insert((i, j + 1), (i, j));
            }
        }
    }
    let mut loops = Vec::new();
    while let Some((&start, _)) = next.iter().next() {
        let mut pts = vec![start];
        let mut at = start;
        while let Some(n) = next.remove(&at) {
            if n == start {
                break;
            }
            pts.push(n);
            at = n;
        }
        loops.push(pts);
    }
    let mut out: Vec<Loop> = loops.into_iter().map(|l| simplify(&l.iter().map(|&(a, b)| (a as f64 * h, b as f64 * h)).collect::<Vec<_>>(), tol)).filter(|l| l.len() >= 3).collect();
    out.sort_by(|a, b| area(b).abs().total_cmp(&area(a).abs()));
    out
}

pub fn area(l: &[(f64, f64)]) -> f64 {
    let n = l.len();
    (0..n).map(|k| l[k].0 * l[(k + 1) % n].1 - l[(k + 1) % n].0 * l[k].1).sum::<f64>() / 2.0
}

fn dist_to_line(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (dx, dy) = (b.0 - a.0, b.1 - a.1);
    let len = dx.hypot(dy);
    if len == 0.0 { (p.0 - a.0).hypot(p.1 - a.1) } else { ((p.0 - a.0) * dy - (p.1 - a.1) * dx).abs() / len }
}

fn dp(pts: &[(f64, f64)], tol: f64, keep: &mut Vec<bool>, lo: usize, hi: usize) {
    if hi <= lo + 1 {
        return;
    }
    let (mut worst, mut at) = (0.0, lo);
    for k in lo + 1..hi {
        let d = dist_to_line(pts[k], pts[lo], pts[hi]);
        if d > worst {
            worst = d;
            at = k;
        }
    }
    if worst > tol {
        keep[at] = true;
        dp(pts, tol, keep, lo, at);
        dp(pts, tol, keep, at, hi);
    }
}

/// Removes vertices that lie exactly on the line between their neighbours.
fn drop_collinear(pts: &[(f64, f64)]) -> Loop {
    let n = pts.len();
    (0..n)
        .filter(|&k| {
            let (a, b, c) = (pts[(k + n - 1) % n], pts[k], pts[(k + 1) % n]);
            ((b.0 - a.0) * (c.1 - b.1) - (b.1 - a.1) * (c.0 - b.0)).abs() > 1e-9
        })
        .map(|k| pts[k])
        .collect()
}

/// Douglas-Peucker on a closed loop, anchored at the first vertex and the one farthest from it.
pub fn simplify(pts: &[(f64, f64)], tol: f64) -> Loop {
    let pts = &drop_collinear(pts);
    let n = pts.len();
    if n < 4 {
        return pts.to_vec();
    }
    let far = (1..n).max_by(|&a, &b| ((pts[a].0 - pts[0].0).hypot(pts[a].1 - pts[0].1)).total_cmp(&(pts[b].0 - pts[0].0).hypot(pts[b].1 - pts[0].1))).unwrap_or(1);
    let mut keep = vec![false; n + 1];
    let mut ring = pts.to_vec();
    ring.push(pts[0]);
    keep[0] = true;
    keep[far] = true;
    keep[n] = true;
    dp(&ring, tol, &mut keep, 0, far);
    dp(&ring, tol, &mut keep, far, n);
    (0..n).filter(|&k| keep[k]).map(|k| pts[k]).collect()
}

/// A picture of the design: the outline, filled, as SVG (y up).
pub fn svg(loops: &[Loop], length: f64, depth: f64) -> String {
    let mut path = String::new();
    for l in loops {
        for (k, (x, y)) in l.iter().enumerate() {
            path.push_str(&format!("{}{:.1} {:.1} ", if k == 0 { "M" } else { "L" }, x, depth - y));
        }
        path.push_str("Z ");
    }
    format!("<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"-20 -20 {:.0} {:.0}\" width=\"{:.0}\" height=\"{:.0}\"><rect x=\"-20\" y=\"-20\" width=\"{:.0}\" height=\"{:.0}\" fill=\"#000\"/><path d=\"{}\" fill=\"#fff\" fill-rule=\"evenodd\"/></svg>\n", length + 40.0, depth + 40.0, (length + 40.0) / 2.0, (depth + 40.0) / 2.0, length + 40.0, depth + 40.0, path)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]
    use super::*;

    fn spec() -> BeamSpec {
        BeamSpec { length: 1000.0, depth: 100.0, tip_depth: 100.0, min_member_mm: 0.0, elem: 5.0, e_mpa: 68_900.0, nu: 0.33, density_g_cm3: 2.7, yield_mpa: 276.0, load_n: 100.0, pad_mm: 20.0, safety: 2.0, buckling_safety: 3.0, max_defl_mm: 50.0, rmin_cells: 2.0, iterations: 10, self_weight: false, member_len_mm: 300.0, gpu: false }
    }

    /// A solid bar: tip deflection P L^3 / (3 E I) plus shear, root stress M c / I, and the
    /// self-weight load equal to the bar's weight.
    #[test]
    fn solid_bar_matches_beam_theory() {
        let mut s = spec();
        let (nx, ny) = (s.nx(), s.ny());
        let t = 10.0;
        let solid = vec![true; nx * ny];
        let d = evaluate(&s, &solid, true, t).unwrap();
        let i = t * s.depth.powi(3) / 12.0;
        let g = s.e_mpa / (2.0 * (1.0 + s.nu));
        let theory = s.load_n * s.length.powi(3) / (3.0 * s.e_mpa * i) + s.load_n * s.length / (5.0 / 6.0 * g * t * s.depth);
        assert!((d.tip_defl_no_weight_mm / theory - 1.0).abs() < 0.04, "tip {} theory {theory}", d.tip_defl_no_weight_mm);
        let root = s.load_n * s.length * (s.depth / 2.0) / i;
        assert!((d.vm_max / root - 1.0).abs() < 0.15, "stress {} theory {root}", d.vm_max);
        s.self_weight = true;
        let w = evaluate(&s, &solid, true, t).unwrap();
        let expected = 2.7e-6 * s.length * s.depth * t * 9.81;
        assert!((w.weight_n / expected - 1.0).abs() < 1e-9);
        assert!(w.tip_defl_mm > d.tip_defl_mm);
        assert!((w.mass_kg - 2.7e-6 * s.length * s.depth * t).abs() < 1e-9);
        // A solid 100 mm deep bar: the thinnest ridge width is the depth.
        assert!((d.member_width_mm - 100.0).abs() < 10.0, "width {}", d.member_width_mm);
    }

    #[test]
    fn opening_removes_thin_bars_and_keeps_wide_ones() {
        let (nx, ny) = (40, 30);
        // Rows 4-5 (2 cells wide) and rows 12-15 (4 cells wide), both attached to the fixed edge.
        let solid: Vec<bool> = (0..nx * ny).map(|k| matches!(k % ny, 4 | 5 | 12..=15)).collect();
        let out = open_design(&solid, nx, ny, 3); // keep members of 3 cells or more
        // A 2-cell window keeps the 2-cell bar and removes a 1-cell one.
        let two = open_design(&solid, nx, ny, 2);
        assert!((0..nx).all(|i| two[i * ny + 4] && two[i * ny + 5]), "2-cell bar lost by a 2-cell window");
        let one_cell: Vec<bool> = (0..nx * ny).map(|k| k % ny == 7).collect();
        assert!(open_design(&one_cell, nx, ny, 2).iter().all(|b| !b), "1-cell bar kept by a 2-cell window");
        for i in 0..nx {
            for j in 0..ny {
                let was = solid[i * ny + j];
                let now = out[i * ny + j];
                if matches!(j, 4 | 5) {
                    assert!(!now, "thin bar kept at ({i},{j})");
                } else {
                    assert_eq!(was, now, "wide bar changed at ({i},{j})");
                }
            }
        }
    }

    #[test]
    fn member_width_of_a_bar_of_known_width() {
        let (nx, ny) = (60, 20);
        let solid: Vec<bool> = (0..nx * ny).map(|k| (5..9).contains(&(k % ny))).collect();
        let w = member_width(&solid, nx, ny, 2.0);
        // Cell centres cannot sit on the ridge of an even-width bar: the result is low by up to
        // one cell (2 mm here), which errs on the safe side for the buckling check.
        assert!((6.0 - 1e-9..=8.0 + 1e-9).contains(&w), "width {w}");
    }

    #[test]
    fn outline_of_a_ring_has_an_outer_loop_and_a_hole() {
        let (nx, ny) = (12, 12);
        let mut solid = vec![false; nx * ny];
        for i in 0..nx {
            for j in 0..ny {
                solid[i * ny + j] = !(4..8).contains(&i) || !(4..8).contains(&j);
            }
        }
        let d = Design { nx, ny, solid, volfrac: 0.0, thickness: 1.0, mass_kg: 0.0, tip_defl_mm: 0.0, tip_defl_no_weight_mm: 0.0, vm_max: 0.0, vm_p99: 0.0, weight_n: 0.0, ltb_pcr_n: 0.0, connected: true, member_width_mm: 0.0, member_margin: 0.0 };
        let loops = outline(&d, 2.0, 0.1);
        assert_eq!(loops.len(), 2);
        assert!((area(&loops[0]) - 576.0).abs() < 1e-9 && (area(&loops[1]) + 64.0).abs() < 1e-9);
        assert_eq!(loops[0].len(), 4);
    }

    #[test]
    fn diagonal_contacts_are_closed_and_islands_dropped() {
        let (nx, ny) = (6, 6);
        let mut dens = vec![0.0; nx * ny];
        dens.iter_mut().take(ny).for_each(|d| *d = 1.0); // fixed edge column
        dens[ny + 5] = 1.0; // touches the column's top cell
        dens[2 * ny + 4] = 1.0; // diagonal to (1,5)
        dens[5 * ny + 1] = 1.0; // island
        let (solid, _) = clean(&dens, nx, ny, 1, ny);
        assert!(!solid[5 * ny + 1], "island kept");
        assert!(solid[ny + 4] || solid[2 * ny + 5], "diagonal contact left open");
    }
}
