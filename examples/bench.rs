#![allow(clippy::indexing_slicing)] // serde_json::Value indexing returns Null, it does not panic
//! Native run of the same benchmark, for comparison with the WebAssembly timing in SolidCraft.
//! cargo run --release --example bench -- NX NY ITERATIONS [cantilever|mbb]
use generative_design_engine::engine::{BenchParams, bench_report};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let num = |i: usize, d: usize| args.get(i).and_then(|s| s.parse().ok()).unwrap_or(d);
    let b = BenchParams { nx: num(0, 60), ny: num(1, 30), iterations: num(2, 50), load: args.get(3).cloned().unwrap_or_else(|| "cantilever".into()), ..Default::default() };
    let p = match b.problem::<()>() {
        Ok(p) => p,
        Err(e) => {
            eprintln!("bad parameters: {e:?}");
            std::process::exit(2);
        }
    };
    let t = Instant::now();
    let out = generative_design_engine::topo2d::optimise(&p, |_, _| true);
    let ms = t.elapsed().as_secs_f64() * 1000.0;
    let r = bench_report(&p, &out, true);
    println!("{}x{} {} its: {:.0} ms native, cg total {}, compliance {} -> {}", p.nx, p.ny, p.max_iterations, ms, r["cg_iterations_total"], r["compliance_first"], r["compliance_last"]);
    if let Some(rows) = r["picture"].as_array() {
        for row in rows {
            println!("  |{}|", row.as_str().unwrap_or(""));
        }
    }
}
