#![allow(clippy::unwrap_used)]
//! Peak stress of the 100 x 10 x 5 cantilever (beam theory 120 MPa at the root) against the grid.
use generative_design_engine::stress2d::{Input, PlaneChoice, box_mesh, check};
fn main() {
    let mesh = box_mesh([100.0, 10.0, 5.0]);
    for res in [100, 150, 200, 300, 400] {
        let input = Input { positions: &mesh.0, triangles: &mesh.1, faces: &mesh.2, plane: PlaneChoice::Auto, resolution: res, fixed: "-x", load: "+x", force: [0.0, -100.0, 0.0], e_mpa: 68_900.0, nu: 0.33 };
        let out = check(&input).unwrap();
        println!("res {res:3} peak {:6.2} at x={:.2} y={:.2}  p99 {:6.2} tip {:.4}", out.max_vm, out.max_vm_at[0], out.max_vm_at[1], out.p99_vm, out.max_disp);
    }
}
