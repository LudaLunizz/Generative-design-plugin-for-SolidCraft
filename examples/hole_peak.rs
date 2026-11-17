#![allow(clippy::unwrap_used)]
//! Peak stress at a round hole against the grid: plate 120 x 40 x 5, hole 10 mm, 3000 N pull.
//! Textbook (Peterson, d/w = 0.25): Kt on the net section 2.42, net stress 20 MPa -> about 48.4.
use generative_design_engine::stress2d::{Input, PlaneChoice, box_mesh, check, plate_with_hole_mesh};

fn main() {
    let (_, _, faces) = box_mesh([120.0, 40.0, 5.0]);
    for shift in [0.0, 0.07, 0.13, 0.21] {
        println!("hole centre shifted by {shift} mm");
        for res in [100, 150, 200, 240, 300, 400] {
            let (pos, tri) = plate_with_hole_mesh([120.0, 40.0, 5.0], [60.0 + shift, 20.0 + shift * 0.6], 5.0, 96);
            let input = Input { positions: &pos, triangles: &tri, faces: &faces, plane: PlaneChoice::Auto, resolution: res, fixed: "-x", load: "+x", force: [3000.0, 0.0, 0.0], e_mpa: 68_900.0, nu: 0.33 };
            let out = check(&input).unwrap();
            println!("  res {res:3} cell {:.3}  peak {:6.2}  p99 {:6.2}  iters {}", out.cell, out.max_vm, out.p99_vm, out.cg_iterations);
        }
    }
}
