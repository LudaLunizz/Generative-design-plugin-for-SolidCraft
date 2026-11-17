#![allow(clippy::unwrap_used, clippy::indexing_slicing)]
//! Times the pieces of one multigrid-CG step natively: stiffness product, V-cycle, setup.
use generative_design_engine::mg2d::{Banded, Cells, Hierarchy, Level, Settings};
use generative_design_engine::topo2d::element_stiffness;
use std::time::Instant;

fn main() {
    let (nx, ny) = (120, 60);
    let mut fixed = vec![false; 2 * (nx + 1) * (ny + 1)];
    for j in 0..=ny {
        fixed[2 * j] = true;
        fixed[2 * j + 1] = true;
    }
    let modulus: Vec<f64> = (0..nx * ny).map(|i| if i % 7 == 0 { 1e-9 } else { 0.5 }).collect();
    let level = || Level { nx, ny, cells: Cells::Scaled { ke: Box::new(element_stiffness(0.3)), modulus: modulus.clone() }, fixed: fixed.clone() };
    let fine = level();
    let n = fine.ndof();
    let u: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin()).collect();
    let mut y = vec![0.0; n];
    let reps = 2000;
    let t = Instant::now();
    for _ in 0..reps {
        fine.apply(&u, &mut y);
    }
    let apply_us = t.elapsed().as_secs_f64() * 1e6 / reps as f64;
    let ke = element_stiffness(0.3);
    let mut y2 = vec![0.0; n];
    let t = Instant::now();
    for _ in 0..reps {
        plain(nx, ny, &ke, &modulus, &u, &mut y2);
    }
    let plain_us = t.elapsed().as_secs_f64() * 1e6 / reps as f64;
    let mut y3 = y2.clone();
    for (v, f) in y3.iter_mut().zip(&fixed) { if *f { *v = 0.0; } }
    let diff = y3.iter().zip(&y).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    println!("plain indexing: {plain_us:.0} us ({:.1} ns per cell), max diff {diff:e}", plain_us * 1e3 / (nx * ny) as f64);
    let t = Instant::now();
    let h = Hierarchy::new(level(), Settings::default());
    let setup_us = t.elapsed().as_secs_f64() * 1e6;
    let t = Instant::now();
    for _ in 0..200 {
        let _ = h.precondition(&u);
    }
    let vcycle_us = t.elapsed().as_secs_f64() * 1e6 / 200.0;
    let t = Instant::now();
    let c1 = fine.coarsen().unwrap();
    let c2 = c1.coarsen().unwrap();
    let coarsen_us = t.elapsed().as_secs_f64() * 1e6;
    let t = Instant::now();
    let _f = Banded::factor(&c2);
    let factor_us = t.elapsed().as_secs_f64() * 1e6;
    let t = Instant::now();
    let _d = (fine.diagonal(), c1.diagonal(), c2.diagonal());
    let diag_us = t.elapsed().as_secs_f64() * 1e6;
    println!("setup pieces: coarsen x2 {coarsen_us:.0} us, banded factor ({} dofs) {factor_us:.0} us, diagonals {diag_us:.0} us", c2.ndof());
    println!("{nx}x{ny}: apply {apply_us:.0} us ({:.1} ns per cell), V-cycle {vcycle_us:.0} us (= {:.1} applies), setup {setup_us:.0} us, levels {}", apply_us * 1e3 / (nx * ny) as f64, vcycle_us / apply_us, h.depth());
}

#[allow(dead_code)]
pub fn plain(nx: usize, ny: usize, ke: &[[f64; 8]; 8], m: &[f64], u: &[f64], y: &mut [f64]) {
    y.iter_mut().for_each(|v| *v = 0.0);
    for ex in 0..nx {
        for ey in 0..ny {
            let d = generative_design_engine::mg2d::cell_dofs(ny, ex, ey);
            let mut ye = [0.0; 8];
            for j in 0..8 {
                let v = u[d[j]];
                for i in 0..8 {
                    ye[i] += ke[j][i] * v;
                }
            }
            let s = m[ex * ny + ey];
            for i in 0..8 {
                y[d[i]] += s * ye[i];
            }
        }
    }
}
