#![allow(clippy::unwrap_used, clippy::indexing_slicing)]
//! Look at the raw optimiser output of the beam problem: beam_look VOLFRAC [ELEM] [RMIN] [ITS]
use generative_design_engine::topo2d::{self, LoadCase, Problem, Solver};
fn main() {
    let a: Vec<f64> = std::env::args().skip(1).filter_map(|v| v.parse().ok()).collect();
    let (vf, elem, rmin, its) = (a.first().copied().unwrap_or(0.3), a.get(1).copied().unwrap_or(10.0), a.get(2).copied().unwrap_or(2.5), a.get(3).copied().unwrap_or(80.0));
    let (nx, ny) = ((5000.0 / elem) as usize, (400.0 / elem) as usize);
    let p = Problem { nx, ny, volfrac: vf, penal: 3.0, rmin, nu: 0.33, load: LoadCase::TipPad((100.0 / elem) as usize), max_iterations: its as usize, change_tol: 0.01, cg_tol: 1e-6, cg_max: 5000, solver: Solver::Multigrid, beta_final: a.get(4).copied().unwrap_or(16.0), beta_every: 40, tip_rows: 0, coarse_budget: 12_000_000, pad_rows: 0 };
    let out = topo2d::optimise(&p, |_, _| true);
    for (k, i) in out.history.iter().enumerate().filter(|(k, _)| k % 10 == 0 || *k + 1 == out.history.len()) {
        println!("it {k:3} compliance {:.4e} vol {:.3} change {:.3} cg {} ok {}", i.compliance, i.volume, i.change, i.cg_iterations, i.cg_converged);
    }
    let grey = out.density.iter().filter(|d| (0.1..0.9).contains(*d)).count() as f64 / out.density.len() as f64;
    println!("converged {} grey fraction {:.2}", out.converged, grey);
    for row in topo2d::ascii(&out.density, nx, ny, 200) {
        println!("{row}");
    }
}
