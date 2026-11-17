//! 2D topology optimisation: linear elastic finite elements on a regular grid (bilinear plane
//! stress quads, unit thickness, square cells of size 1) and the SIMP method (density per cell,
//! stiffness = density^penal, density filter, optimality-criteria update).
//!
//! Written from the textbook formulation: element stiffness by 2x2 Gauss quadrature, the global
//! system solved matrix-free with Jacobi-preconditioned conjugate gradients (no assembled matrix,
//! so memory stays a few vectors long), warm-started from the previous design's displacements.
//!
//! The solver is matrix-free conjugate gradients, preconditioned with Jacobi or a geometric
//! multigrid V-cycle (`mg2d`).
//!
//! Everything is bounded by `Problem::check`; nothing here panics on any input that passes it.

use crate::mg2d::{Cells, Hierarchy, Level, Settings, cell_dofs};

/// One element's 8x8 stiffness for Young's modulus 1, in DOF order
/// (n0x, n0y, n1x, n1y, n2x, n2y, n3x, n3y) with nodes counter-clockwise from bottom-left.
pub fn element_stiffness(nu: f64) -> [[f64; 8]; 8] {
    let c = 1.0 / (1.0 - nu * nu);
    let d = [[c, c * nu, 0.0], [c * nu, c, 0.0], [0.0, 0.0, c * (1.0 - nu) / 2.0]];
    let corners = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];
    let g = 1.0 / 3f64.sqrt();
    let mut k = [[0.0; 8]; 8];
    for (gx, gy) in [(-g, -g), (g, -g), (g, g), (-g, g)] {
        // Strain-displacement rows (exx, eyy, gxy), one column per DOF. Unit square:
        // x = (xi + 1) / 2, so d/dx = 2 d/dxi and det J = 1/4.
        let mut b = [[0.0; 8]; 3];
        for (n, (xi, eta)) in corners.iter().enumerate() {
            let dndx = 2.0 * 0.25 * xi * (1.0 + eta * gy);
            let dndy = 2.0 * 0.25 * eta * (1.0 + xi * gx);
            let [bx, by, bxy] = &mut b;
            for (row, col, value) in [(&mut *bx, 2 * n, dndx), (&mut *by, 2 * n + 1, dndy)] {
                if let Some(v) = row.get_mut(col) {
                    *v = value;
                }
            }
            for (col, value) in [(2 * n, dndy), (2 * n + 1, dndx)] {
                if let Some(v) = bxy.get_mut(col) {
                    *v = value;
                }
            }
        }
        // k += B^T D B * weight(1) * det J(1/4)
        let db: Vec<[f64; 8]> = d
            .iter()
            .map(|drow| {
                let mut out = [0.0; 8];
                for (j, o) in out.iter_mut().enumerate() {
                    *o = drow.iter().zip(&b).map(|(dv, brow)| dv * brow.get(j).copied().unwrap_or(0.0)).sum();
                }
                out
            })
            .collect();
        for (i, krow) in k.iter_mut().enumerate() {
            for (j, kv) in krow.iter_mut().enumerate() {
                let s: f64 = b.iter().zip(&db).map(|(brow, dbrow)| brow.get(i).copied().unwrap_or(0.0) * dbrow.get(j).copied().unwrap_or(0.0)).sum();
                *kv += s * 0.25;
            }
        }
    }
    k
}

/// Where the part is held and loaded.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LoadCase {
    /// Left edge fixed, unit load straight down at the middle of the right edge.
    Cantilever,
    /// Half of a simply supported beam (MBB): left edge symmetric (x fixed), bottom-right corner
    /// on a roller (y fixed), unit load down at the top-left corner.
    Mbb,
}

#[derive(Clone, Debug)]
pub struct Problem {
    pub nx: usize,
    pub ny: usize,
    pub volfrac: f64,
    pub penal: f64,
    pub rmin: f64,
    pub nu: f64,
    pub load: LoadCase,
    pub max_iterations: usize,
    pub change_tol: f64,
    pub cg_tol: f64,
    pub cg_max: usize,
    pub solver: Solver,
}

pub const MAX_ELEMENTS: usize = 250_000;

impl Problem {
    pub fn check(&self) -> Result<(), String> {
        let range = |v: f64, lo: f64, hi: f64, name: &str| {
            if v.is_finite() && (lo..=hi).contains(&v) { Ok(()) } else { Err(format!("{name}: expected {lo} to {hi}")) }
        };
        if !(2..=2000).contains(&self.nx) || !(2..=2000).contains(&self.ny) {
            return Err("nx and ny: expected 2 to 2000 cells".into());
        }
        if self.nx.saturating_mul(self.ny) > MAX_ELEMENTS {
            return Err(format!("nx * ny: at most {MAX_ELEMENTS} cells"));
        }
        range(self.volfrac, 0.05, 0.95, "volfrac")?;
        range(self.penal, 1.0, 5.0, "penal")?;
        range(self.rmin, 1.0, 10.0, "rmin")?;
        range(self.nu, 0.0, 0.49, "nu")?;
        range(self.change_tol, 0.0, 1.0, "change_tol")?;
        range(self.cg_tol, 1e-12, 1e-1, "cg_tol")?;
        if !(1..=1000).contains(&self.max_iterations) || !(1..=100_000).contains(&self.cg_max) {
            return Err("max_iterations 1-1000, cg_max 1-100000".into());
        }
        Ok(())
    }

    fn nodes_y(&self) -> usize {
        self.ny + 1
    }
    pub fn ndof(&self) -> usize {
        2 * (self.nx + 1) * (self.ny + 1)
    }
    fn element_dofs(&self, ex: usize, ey: usize) -> [usize; 8] {
        cell_dofs(self.ny, ex, ey)
    }
    /// Fixed DOFs as a mask, and the load vector.
    fn supports_and_load(&self) -> (Vec<bool>, Vec<f64>) {
        let n = self.nodes_y();
        let mut fixed = vec![false; self.ndof()];
        let mut f = vec![0.0; self.ndof()];
        let set = |v: &mut Vec<bool>, d: usize| {
            if let Some(slot) = v.get_mut(d) {
                *slot = true;
            }
        };
        match self.load {
            LoadCase::Cantilever => {
                for j in 0..n {
                    set(&mut fixed, 2 * j);
                    set(&mut fixed, 2 * j + 1);
                }
                if let Some(v) = f.get_mut(2 * (self.nx * n + self.ny / 2) + 1) {
                    *v = -1.0;
                }
            }
            LoadCase::Mbb => {
                for j in 0..n {
                    set(&mut fixed, 2 * j);
                }
                set(&mut fixed, 2 * (self.nx * n) + 1);
                if let Some(v) = f.get_mut(2 * self.ny + 1) {
                    *v = -1.0;
                }
            }
        }
        (fixed, f)
    }
}

/// How the stiffness system is solved: conjugate gradients with a preconditioner.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Solver {
    /// Diagonal (Jacobi) preconditioner: cheap per step, many steps on fine grids.
    Jacobi,
    /// One geometric multigrid V-cycle per step: a few dozen steps at any grid size.
    Multigrid,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Preconditioned conjugate gradients from the start value in `u`. Returns the iterations.
fn pcg(op: &Level, precondition: &dyn Fn(&[f64]) -> Vec<f64>, f: &[f64], u: &mut [f64], tol: f64, max: usize) -> (usize, bool) {
    let n = f.len();
    for (v, fixed) in u.iter_mut().zip(&op.fixed) {
        if *fixed {
            *v = 0.0;
        }
    }
    let mut r = vec![0.0; n];
    op.apply(u, &mut r);
    for ((ri, fi), fixed) in r.iter_mut().zip(f).zip(&op.fixed) {
        *ri = if *fixed { 0.0 } else { fi - *ri };
    }
    let fnorm = dot(f, f).sqrt().max(f64::MIN_POSITIVE);
    let mut z = precondition(&r);
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    let mut q = vec![0.0; n];
    for it in 0..max {
        if dot(&r, &r).sqrt() <= tol * fnorm {
            return (it, true);
        }
        op.apply(&p, &mut q);
        let pq = dot(&p, &q);
        if !(pq.is_finite() && pq > 0.0) {
            return (it, false);
        }
        let alpha = rz / pq;
        for ((ui, pi), (ri, qi)) in u.iter_mut().zip(&p).zip(r.iter_mut().zip(&q)) {
            *ui += alpha * pi;
            *ri -= alpha * qi;
        }
        z = precondition(&r);
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz;
        rz = rz_new;
        for (pi, zi) in p.iter_mut().zip(&z) {
            *pi = zi + beta * *pi;
        }
    }
    (max, dot(&r, &r).sqrt() <= tol * fnorm)
}

/// Static analysis of one design: displacements and compliance (f . u).
pub struct Analysis {
    pub u: Vec<f64>,
    pub compliance: f64,
    pub cg_iterations: usize,
    pub converged: bool,
    /// Multigrid levels used (1 for Jacobi).
    pub levels: usize,
    /// Seconds spent building the preconditioner and in CG (0 without a clock).
    pub setup_s: f64,
    pub solve_s: f64,
}

pub fn analyse(p: &Problem, modulus: &[f64], start: Option<Vec<f64>>) -> Analysis {
    analyse_timed(p, modulus, start, &|| 0.0)
}

/// [`analyse`] with a clock (seconds) for the setup/solve split.
pub fn analyse_timed(p: &Problem, modulus: &[f64], start: Option<Vec<f64>>, clock: &dyn Fn() -> f64) -> Analysis {
    let t0 = clock();
    let (fixed, f) = p.supports_and_load();
    let fine = Level { nx: p.nx, ny: p.ny, cells: Cells::Scaled { ke: Box::new(element_stiffness(p.nu)), modulus: modulus.to_vec() }, fixed };
    let mut u = start.filter(|s| s.len() == f.len()).unwrap_or_else(|| vec![0.0; f.len()]);
    let (cg_iterations, converged, levels, t1) = match p.solver {
        Solver::Jacobi => {
            let inv: Vec<f64> = fine.diagonal().iter().map(|d| 1.0 / d).collect();
            let jacobi = |r: &[f64]| r.iter().zip(&inv).map(|(a, d)| a * d).collect();
            let t1 = clock();
            let (it, ok) = pcg(&fine, &jacobi, &f, &mut u, p.cg_tol, p.cg_max);
            (it, ok, 1, t1)
        }
        Solver::Multigrid => {
            let h = Hierarchy::new(fine, Settings::default());
            let t1 = clock();
            match h.fine() {
                Some(op) => {
                    let (it, ok) = pcg(op, &|r: &[f64]| h.precondition(r), &f, &mut u, p.cg_tol, p.cg_max);
                    (it, ok, h.depth(), t1)
                }
                None => (0, false, 0, t1),
            }
        }
    };
    let t2 = clock();
    let compliance = dot(&f, &u);
    Analysis { u, compliance, cg_iterations, converged, levels, setup_s: t1 - t0, solve_s: t2 - t1 }
}

pub struct Iteration {
    pub compliance: f64,
    pub volume: f64,
    pub change: f64,
    pub cg_iterations: usize,
    pub cg_converged: bool,
    pub levels: usize,
}

pub struct Outcome {
    pub density: Vec<f64>, // filtered (physical) density per cell, column-major (ex, ey)
    pub history: Vec<Iteration>,
    pub converged: bool,
    pub cancelled: bool,
    /// Seconds per phase over the whole run: preconditioner setup, CG solve, sensitivities,
    /// filter + optimality-criteria update (all 0 without a clock).
    pub phases: Phases,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Phases {
    pub setup_s: f64,
    pub solve_s: f64,
    pub sensitivity_s: f64,
    pub update_s: f64,
}

/// Density filter: a cone of radius `rmin` around each cell, applied as `out = H in / Hs`.
struct Filter {
    nx: usize,
    ny: usize,
    stencil: Vec<(isize, isize, f64)>,
    hs: Vec<f64>,
}

impl Filter {
    fn new(nx: usize, ny: usize, r: f64) -> Self {
        let reach = (r.ceil() as usize).saturating_sub(1);
        // The cone's weights, once: (dx, dy, w) for every offset with w > 0.
        let reach_i = reach as isize;
        let mut stencil = Vec::new();
        for dx in -reach_i..=reach_i {
            for dy in -reach_i..=reach_i {
                let w = r - ((dx * dx + dy * dy) as f64).sqrt();
                if w > 0.0 {
                    stencil.push((dx, dy, w));
                }
            }
        }
        let mut f = Self { nx, ny, stencil, hs: Vec::new() };
        let ones = vec![1.0; nx * ny];
        f.hs = f.weighted(&ones);
        f
    }
    fn weighted(&self, v: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; self.nx * self.ny];
        let (nx, ny) = (self.nx as isize, self.ny as isize);
        for (i, column) in out.chunks_exact_mut(self.ny.max(1)).enumerate() {
            for (j, o) in column.iter_mut().enumerate() {
                let (i, j) = (i as isize, j as isize);
                let mut s = 0.0;
                for (dx, dy, w) in &self.stencil {
                    let (k, l) = (i + dx, j + dy);
                    if k >= 0 && k < nx && l >= 0 && l < ny {
                        s += w * v.get((k * ny + l) as usize).copied().unwrap_or(0.0);
                    }
                }
                *o = s;
            }
        }
        out
    }
    /// Filtered densities: H x / Hs.
    fn forward(&self, x: &[f64]) -> Vec<f64> {
        self.weighted(x).iter().zip(&self.hs).map(|(a, h)| a / h.max(f64::MIN_POSITIVE)).collect()
    }
    /// The chain rule through the filter: H (g / Hs). H is symmetric.
    fn backward(&self, g: &[f64]) -> Vec<f64> {
        let scaled: Vec<f64> = g.iter().zip(&self.hs).map(|(a, h)| a / h.max(f64::MIN_POSITIVE)).collect();
        self.weighted(&scaled)
    }
}

/// Runs SIMP. `step` is called after every iteration (iteration index, compliance) and returns
/// false to stop (cancel).
pub fn optimise(p: &Problem, step: impl FnMut(usize, &Iteration) -> bool) -> Outcome {
    optimise_timed(p, &|| 0.0, step)
}

/// [`optimise`] with a clock (seconds) for the phase split in `Outcome::phases`.
pub fn optimise_timed(p: &Problem, clock: &dyn Fn() -> f64, mut step: impl FnMut(usize, &Iteration) -> bool) -> Outcome {
    let mut phases = Phases::default();
    const E0: f64 = 1.0;
    const EMIN: f64 = 1e-9;
    let n = p.nx * p.ny;
    let ke = element_stiffness(p.nu);
    let filter = Filter::new(p.nx, p.ny, p.rmin);
    let mut x = vec![p.volfrac; n];
    let mut phys = filter.forward(&x);
    let mut u: Option<Vec<f64>> = None;
    let mut history = Vec::new();
    let mut converged = false;
    let mut cancelled = false;
    let dv = filter.backward(&vec![1.0; n]);
    for it in 0..p.max_iterations {
        let modulus: Vec<f64> = phys.iter().map(|d| EMIN + d.powf(p.penal) * (E0 - EMIN)).collect();
        let a = analyse_timed(p, &modulus, u.take(), clock);
        phases.setup_s += a.setup_s;
        phases.solve_s += a.solve_s;
        let ts = clock();
        // Element strain energy for E = 1, then the compliance sensitivity.
        let mut dc = vec![0.0; n];
        let mut e = 0;
        for ex in 0..p.nx {
            for ey in 0..p.ny {
                let dofs = p.element_dofs(ex, ey);
                let ue: Vec<f64> = dofs.iter().map(|d| a.u.get(*d).copied().unwrap_or(0.0)).collect();
                let mut ce = 0.0;
                for (row, ui) in ke.iter().zip(&ue) {
                    ce += ui * row.iter().zip(&ue).map(|(k, uj)| k * uj).sum::<f64>();
                }
                let d = phys.get(e).copied().unwrap_or(0.0);
                if let Some(slot) = dc.get_mut(e) {
                    *slot = -p.penal * d.powf(p.penal - 1.0) * (E0 - EMIN) * ce;
                }
                e += 1;
            }
        }
        let tu = clock();
        phases.sensitivity_s += tu - ts;
        let dc = filter.backward(&dc);
        // Optimality criteria with bisection on the volume multiplier.
        let mv = 0.2;
        let (mut l1, mut l2) = (0.0_f64, 1e9_f64);
        let mut xnew = x.clone();
        for _ in 0..200 {
            if (l2 - l1) / (l1 + l2).max(f64::MIN_POSITIVE) <= 1e-4 {
                break;
            }
            let lmid = 0.5 * (l1 + l2);
            for i in 0..n {
                let xi = x.get(i).copied().unwrap_or(0.0);
                let ratio = (-dc.get(i).copied().unwrap_or(0.0) / (dv.get(i).copied().unwrap_or(1.0) * lmid)).max(0.0);
                let v = (xi * ratio.sqrt()).clamp((xi - mv).max(0.0), (xi + mv).min(1.0));
                if let Some(slot) = xnew.get_mut(i) {
                    *slot = if v.is_finite() { v } else { xi };
                }
            }
            // The filtered volume is linear in x: sum_i (H x)_i / Hs_i = x . dv, with dv the
            // filter's backward image of ones. No filter pass needed per bisection step.
            let vol = dot(&xnew, &dv);
            if vol > p.volfrac * n as f64 { l1 = lmid } else { l2 = lmid }
        }
        let change = x.iter().zip(&xnew).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        x = xnew;
        phys = filter.forward(&x);
        phases.update_s += clock() - tu;
        let iteration = Iteration {
            compliance: a.compliance,
            volume: phys.iter().sum::<f64>() / n as f64,
            change,
            cg_iterations: a.cg_iterations,
            cg_converged: a.converged,
            levels: a.levels,
        };
        u = Some(a.u);
        let go_on = step(it, &iteration);
        history.push(iteration);
        if change < p.change_tol {
            converged = true;
            break;
        }
        if !go_on {
            cancelled = true;
            break;
        }
    }
    Outcome { density: phys, history, converged, cancelled, phases }
}

/// Profiling aid for the benchmark: runs one part of a design iteration `repeat` times on a
/// mid-way design (density alternating 1 and 0.2), so its fuel can be measured on its own.
pub fn profile_part(p: &Problem, part: &str, repeat: usize) -> Result<serde_json::Value, String> {
    let n = p.nx * p.ny;
    let modulus: Vec<f64> = (0..n).map(|i| if i % 3 == 0 { 0.2f64.powf(p.penal) } else { 1.0 }).collect();
    let (fixed, f) = p.supports_and_load();
    let fine = || Level { nx: p.nx, ny: p.ny, cells: Cells::Scaled { ke: Box::new(element_stiffness(p.nu)), modulus: modulus.clone() }, fixed: fixed.clone() };
    let mut checksum = 0.0;
    match part {
        "apply" => {
            let level = fine();
            let mut y = vec![0.0; f.len()];
            let mut u = f.clone();
            for _ in 0..repeat {
                level.apply(&u, &mut y);
                std::mem::swap(&mut u, &mut y);
            }
            checksum = dot(&u, &u);
        }
        "setup" => {
            for _ in 0..repeat {
                checksum += Hierarchy::new(fine(), Settings::default()).depth() as f64;
            }
        }
        "vcycle" => {
            let h = Hierarchy::new(fine(), Settings::default());
            let mut r = f.clone();
            for _ in 0..repeat {
                r = h.precondition(&r);
            }
            checksum = dot(&r, &r);
        }
        "filter" => {
            let filter = Filter::new(p.nx, p.ny, p.rmin);
            let mut x = vec![p.volfrac; n];
            for _ in 0..repeat {
                x = filter.forward(&x);
            }
            checksum = x.iter().sum();
        }
        _ => return Err("expected apply, setup, vcycle or filter".into()),
    }
    Ok(serde_json::json!({ "part": part, "repeat": repeat, "cells": n, "checksum": checksum }))
}

/// A small picture of a density field: one character per block of cells, darker = more material.
pub fn ascii(density: &[f64], nx: usize, ny: usize, max_cols: usize) -> Vec<String> {
    let step = nx.div_ceil(max_cols.max(1)).max(1);
    let ystep = (step * 2).max(1); // characters are about twice as tall as wide
    let shades = [' ', '.', ':', '+', '#', '@'];
    let mut rows = Vec::new();
    let mut j = ny;
    while j > 0 {
        let j0 = j.saturating_sub(ystep);
        let mut line = String::new();
        let mut i = 0;
        while i < nx {
            let (mut s, mut c) = (0.0, 0usize);
            for a in i..(i + step).min(nx) {
                for b in j0..j {
                    s += density.get(a * ny + b).copied().unwrap_or(0.0);
                    c += 1;
                }
            }
            let v = if c == 0 { 0.0 } else { s / c as f64 };
            let k = ((v * (shades.len() - 1) as f64).round() as usize).min(shades.len() - 1);
            line.push(shades.get(k).copied().unwrap_or(' '));
            i += step;
        }
        rows.push(line);
        j = j0;
    }
    rows
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
    use super::*;

    fn beam(nx: usize, ny: usize) -> Problem {
        Problem {
            nx,
            ny,
            volfrac: 0.5,
            penal: 3.0,
            rmin: 1.5,
            nu: 0.3,
            load: LoadCase::Cantilever,
            max_iterations: 60,
            change_tol: 0.01,
            cg_tol: 1e-8,
            cg_max: 20_000,
            solver: Solver::Multigrid,
        }
    }

    #[test]
    fn element_stiffness_is_symmetric_with_three_rigid_body_modes() {
        let k = element_stiffness(0.3);
        for (i, row) in k.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                assert!((v - k[j][i]).abs() < 1e-12);
            }
        }
        // Translations in x and y and a small rotation give no forces.
        let tx = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0];
        let ty = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0];
        let rot = [0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0]; // (x, y) -> (-y, x)
        for mode in [tx, ty, rot] {
            for row in &k {
                assert!(row.iter().zip(mode).map(|(a, b)| a * b).sum::<f64>().abs() < 1e-12);
            }
        }
    }

    /// A solid cantilever (length 100, depth 10, unit thickness, E = 1, unit tip load) against
    /// beam theory with shear: d = F L^3 / (3 E I) + F L / (k G A). Bilinear quads are a little
    /// too stiff in bending, so they must come out within a few percent, from below.
    #[test]
    fn solid_cantilever_matches_beam_theory() {
        let p = beam(100, 10);
        let modulus = vec![1.0; p.nx * p.ny];
        let a = analyse(&p, &modulus, None);
        assert!(a.converged);
        let (l, h, nu): (f64, f64, f64) = (100.0, 10.0, 0.3);
        let i = h.powi(3) / 12.0;
        let g = 1.0 / (2.0 * (1.0 + nu));
        let expected = l.powi(3) / (3.0 * i) + l / (5.0 / 6.0 * g * h);
        let tip = a.compliance; // unit load: compliance = tip deflection
        let ratio = tip / expected;
        assert!((0.93..=1.01).contains(&ratio), "tip {tip}, theory {expected}, ratio {ratio}");
    }

    #[test]
    fn simp_lowers_compliance_and_keeps_the_volume() {
        let p = Problem { max_iterations: 40, ..beam(40, 20) };
        let out = optimise(&p, |_, _| true);
        let first = out.history.first().unwrap().compliance;
        let last = out.history.last().unwrap().compliance;
        assert!(last < first, "{first} -> {last}");
        let vol = out.history.last().unwrap().volume;
        assert!((vol - 0.5).abs() < 0.01, "volume {vol}");
        assert!(out.history.iter().all(|i| i.cg_converged));
    }

    #[test]
    fn mbb_beam_runs() {
        let p = Problem { load: LoadCase::Mbb, max_iterations: 20, ..beam(60, 20) };
        let out = optimise(&p, |_, _| true);
        assert!(out.history.last().unwrap().compliance.is_finite());
        assert_eq!(ascii(&out.density, 60, 20, 30).len(), 5);
    }

    #[test]
    fn stop_from_the_callback_is_reported() {
        let p = beam(20, 10);
        let out = optimise(&p, |it, _| it < 2);
        assert!(out.cancelled);
        assert_eq!(out.history.len(), 3);
    }

    #[test]
    fn bad_problems_are_refused() {
        for p in [
            Problem { nx: 1, ..beam(10, 10) },
            Problem { nx: 600, ny: 600, ..beam(10, 10) },
            Problem { volfrac: f64::NAN, ..beam(10, 10) },
            Problem { nu: 0.5, ..beam(10, 10) },
            Problem { cg_max: 0, ..beam(10, 10) },
        ] {
            assert!(p.check().is_err());
        }
        assert!(beam(10, 10).check().is_ok());
    }
}
