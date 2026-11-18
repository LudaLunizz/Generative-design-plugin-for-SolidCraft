#![allow(clippy::unwrap_used, clippy::indexing_slicing)]
//! Where does a beam design iteration spend its time? beam_prof [ELEM] [ITS] [TIP_ROWS_OR_0] [RMIN]
use generative_design_engine::topo2d::{self, LoadCase, Problem, Solver};
use std::time::Instant;
fn main() {
    let a: Vec<f64> = std::env::args().skip(1).filter_map(|v| v.parse().ok()).collect();
    let (elem, its, tip, rmin) = (a.first().copied().unwrap_or(10.0), a.get(1).copied().unwrap_or(20.0), a.get(2).copied().unwrap_or(0.0), a.get(3).copied().unwrap_or(2.0));
    let (nx, ny) = ((5000.0 / elem) as usize, (800.0 / elem) as usize);
    let p = Problem { nx, ny, volfrac: 0.15, penal: 3.0, rmin, nu: 0.33, load: LoadCase::TipPad((100.0 / elem) as usize), max_iterations: its as usize, change_tol: 0.0, cg_tol: a.get(5).copied().unwrap_or(1e-6), cg_max: 5000, solver: Solver::Multigrid, beta_final: 16.0, beta_every: (its as usize / 8).max(5), tip_rows: tip as usize, coarse_budget: a.get(4).copied().unwrap_or(12e6) as usize, pad_rows: 0 };
    let t = Instant::now();
    let clock = || t.elapsed().as_secs_f64();
    let out = topo2d::optimise_timed(&p, &clock, |_, _| true);
    let total = t.elapsed().as_secs_f64();
    let ph = out.phases;
    println!("{nx}x{ny} cells, {its} its: total {total:.1} s = {:.2} s/it", total / its);
    println!("  preconditioner setup {:.1} s, CG solve {:.1} s, sensitivities {:.1} s, filter+OC update {:.1} s", ph.setup_s, ph.solve_s, ph.sensitivity_s, ph.update_s);
    let cg: Vec<usize> = out.history.iter().map(|i| i.cg_iterations).collect();
    println!("  CG steps per solve: first {:?} last {:?} max {} unconverged {}", &cg[..cg.len().min(5)], &cg[cg.len().saturating_sub(5)..], cg.iter().max().unwrap(), out.history.iter().filter(|i| !i.cg_converged).count());
}
