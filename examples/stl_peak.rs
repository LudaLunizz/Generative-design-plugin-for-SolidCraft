#![allow(clippy::unwrap_used, clippy::indexing_slicing)]
//! Stress check of a binary STL plate (fixed -x, load +x), at several resolutions.
//! Usage: stl_peak FILE.stl FORCE_X [RES...]
use generative_design_engine::stress2d::{Input, PlaneChoice, box_mesh, check};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let data = std::fs::read(&args[1]).unwrap();
    let force: f64 = args[2].parse().unwrap();
    let n = u32::from_le_bytes(data[80..84].try_into().unwrap()) as usize;
    let (mut pos, mut tri) = (Vec::new(), Vec::new());
    for k in 0..n {
        let base = 84 + k * 50 + 12;
        let mut ids = [0u32; 3];
        for (v, id) in ids.iter_mut().enumerate() {
            let o = base + v * 12;
            let c = |q: usize| f64::from(f32::from_le_bytes(data[o + 4 * q..o + 4 * q + 4].try_into().unwrap()));
            *id = pos.len() as u32;
            pos.push([c(0), c(1), c(2)]);
        }
        tri.push(ids);
    }
    let lo = (0..3).map(|k| pos.iter().map(|p: &[f64; 3]| p[k]).fold(f64::MAX, f64::min)).collect::<Vec<_>>();
    for p in &mut pos {
        for k in 0..3 {
            p[k] -= lo[k];
        }
    }
    println!("{n} triangles");
    let (_, _, faces) = box_mesh([120.0, 40.0, 5.0]);
    let resolutions: Vec<usize> = if args.len() > 3 { args[3..].iter().map(|r| r.parse().unwrap()).collect() } else { vec![100, 150, 200, 240, 300, 400] };
    for res in resolutions {
        let input = Input { positions: &pos, triangles: &tri, faces: &faces, plane: PlaneChoice::Auto, resolution: res, fixed: "-x", load: "+x", force: [force, 0.0, 0.0], e_mpa: 68_900.0, nu: 0.33 };
        let out = check(&input).unwrap();
        println!("res {res:3} cell {:.3} peak {:6.2} p99 {:6.2} warnings {}", out.cell, out.max_vm, out.p99_vm, out.warnings.len());
    }
}
