#![allow(clippy::unwrap_used, clippy::indexing_slicing)]
//! GPU solve against the CPU solvers on a beam-like design: same displacements, and timings.
//! gpu_check [ELEM] [DESIGN: solid|truss]
use generative_design_engine::gpu::{GpuMg, GpuPcg};
use generative_design_engine::mg2d::Settings;
use generative_design_engine::mg2d::{Cells, Level};
use generative_design_engine::topo2d::{self, LoadCase, Problem, Solver, System, element_stiffness};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let elem: f64 = args.first().and_then(|v| v.parse().ok()).unwrap_or(10.0);
    let truss = args.get(1).map(|s| s == "truss").unwrap_or(true);
    let (nx, ny) = ((5000.0 / elem) as usize, (800.0 / elem) as usize);
    let p = Problem { nx, ny, volfrac: 0.2, penal: 3.0, rmin: 2.0, nu: 0.33, load: LoadCase::TipPad((100.0 / elem) as usize), max_iterations: 1, change_tol: 0.0, cg_tol: 1e-8, cg_max: 200_000, solver: Solver::Multigrid, beta_final: 0.0, beta_every: 30, tip_rows: 0, coarse_budget: 12_000_000, pad_rows: 0 };
    let (fixed, f) = p.supports_and_load();
    // A sharp, truss-like design: two chords and diagonals, void 1e-9 elsewhere (the hard case).
    let modulus: Vec<f64> = (0..nx * ny)
        .map(|e| {
            let (i, j) = (e / ny, e % ny);
            let chord = j < ny / 20 + 1 || j >= ny - ny / 20 - 1;
            let bay = ny; // diagonal every `ny` cells along x
            let diag = (i % bay) as i64 - j as i64;
            let brace = diag.abs() <= (ny / 30 + 1) as i64;
            if !truss || chord || brace { 1.0 } else { 1e-9 }
        })
        .collect();
    println!("{nx} x {ny} cells, {} DOFs, {} design", 2 * (nx + 1) * (ny + 1), if truss { "truss" } else { "solid" });
    let tol = 1e-8;
    let t = Instant::now();
    let cpu = topo2d::solve(&System { nx, ny, nu: 0.33, modulus: &modulus, fixed: fixed.clone(), load: &f, solver: Solver::Multigrid, cg_tol: tol, cg_max: 200_000, coarse_budget: 12_000_000 }, None, &|| 0.0);
    let cpu_ms = t.elapsed().as_secs_f64() * 1e3;
    println!("CPU multigrid PCG: {} steps, {:.0} ms, converged {}", cpu.cg_iterations, cpu_ms, cpu.converged);
    let fine = Level { nx, ny, cells: Cells::Scaled { ke: Box::new(element_stiffness(0.33)), modulus: modulus.clone() }, fixed: fixed.clone() };
    let inv: Vec<f64> = fine.diagonal().iter().map(|d| 1.0 / d).collect();
    // The slow Jacobi comparisons (CPU and GPU) only when asked: JACOBI=1.
    if std::env::var_os("JACOBI").is_some() {
        let t = Instant::now();
        let cpu_j = topo2d::solve(&System { nx, ny, nu: 0.33, modulus: &modulus, fixed: fixed.clone(), load: &f, solver: Solver::Jacobi, cg_tol: tol, cg_max: 200_000, coarse_budget: 0 }, None, &|| 0.0);
        println!("CPU Jacobi PCG:    {} steps, {:.0} ms, converged {}", cpu_j.cg_iterations, t.elapsed().as_secs_f64() * 1e3, cpu_j.converged);
        let gpu = GpuPcg::new(nx, ny, &element_stiffness(0.33), &fixed).unwrap();
        let t = Instant::now();
        let g = gpu.solve(&modulus, &inv, &f, None, tol, 200_000).unwrap();
        println!("GPU Jacobi PCG:    {} steps, {:.0} ms, converged {} ({})", g.iterations, t.elapsed().as_secs_f64() * 1e3, g.converged, gpu.adapter_name);
    }
    let mg = GpuMg::new().unwrap();
    let settings = Settings { coarse_budget: 12_000_000, ..Settings::default() };
    for round in 0..2 {
        let t = Instant::now();
        let (m, times) = mg.solve(nx, ny, 0.33, &modulus, &fixed, &f, None, settings, tol, 20_000).unwrap();
        let ms = t.elapsed().as_secs_f64() * 1e3;
        let scale = cpu.u.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let diff = cpu.u.iter().zip(&m.u).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        println!("GPU multigrid PCG (run {round}): {} steps, {ms:.0} ms (hierarchy {:.0}, upload {:.0}, iterate {:.0} of which coarse CPU {:.0}), converged {}, max diff {:.2e}; GPU per step: up {:.2} ms, down {:.2} ms", m.iterations, times.hierarchy * 1e3, times.upload * 1e3, times.iterate * 1e3, times.coarse_cpu * 1e3, m.converged, diff / scale, times.gpu_up * 1e3, times.gpu_down * 1e3);
    }
}
