#![allow(clippy::indexing_slicing)] // serde_json::Value indexing returns Null, it does not panic
//! Native run of the same benchmark, for comparison with the WebAssembly timing in SolidCraft.
//! cargo run --release --example bench -- NX NY ITERATIONS [cantilever|mbb] [multigrid|jacobi]
use generative_design_engine::engine::{BenchParams, bench_report};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let num = |i: usize, d: usize| args.get(i).and_then(|s| s.parse().ok()).unwrap_or(d);
    let b = BenchParams { nx: num(0, 60), ny: num(1, 30), iterations: num(2, 50), load: args.get(3).cloned().unwrap_or_else(|| "cantilever".into()), solver: args.get(4).cloned().unwrap_or_else(|| "multigrid".into()), cg_tol: std::env::var("GDE_TOL").ok().and_then(|v| v.parse().ok()).unwrap_or(1e-6), ..Default::default() };
    let p = match b.problem::<()>() {
        Ok(p) => p,
        Err(e) => {
            eprintln!("bad parameters: {e:?}");
            std::process::exit(2);
        }
    };
    let t = Instant::now();
    let clock = || t.elapsed().as_secs_f64();
    let out = generative_design_engine::topo2d::optimise_timed(&p, &clock, |_, _| true);
    let ms = t.elapsed().as_secs_f64() * 1000.0;
    let r = bench_report(&p, &out, true);
    println!("{}x{} {} its {} ({} levels): {:.0} ms native, cg total {}, compliance {} -> {}", p.nx, p.ny, p.max_iterations, r["solver"].as_str().unwrap_or("?"), r["multigrid_levels"], ms, r["cg_iterations_total"], r["compliance_first"], r["compliance_last"]);
    if let Some(path) = std::env::var_os("GDE_DUMP") {
        let text: Vec<String> = out.density.iter().map(|d| format!("{d:.6}")).collect();
        let _ = std::fs::write(path, text.join("\n"));
    }
    let ph = out.phases;
    println!("  phases: setup {:.0} ms, solve {:.0} ms, sensitivities {:.0} ms, filter+update {:.0} ms", ph.setup_s * 1e3, ph.solve_s * 1e3, ph.sensitivity_s * 1e3, ph.update_s * 1e3);
    if std::env::var_os("NO_PICTURE").is_some() {
        return;
    }
    if let Some(rows) = r["picture"].as_array() {
        for row in rows {
            println!("  |{}|", row.as_str().unwrap_or(""));
        }
    }
}
