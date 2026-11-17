//! Stress check of a flat part loaded in its own plane: plane stress on a grid.
//!
//! The body's triangle mesh is sampled on a square grid in its plane (the plane across its
//! thinnest direction unless chosen), at mid-thickness, by casting rays through the thickness.
//! Cells inside are material (E x thickness), the rest is near-void (1e-9 of it). Supports and
//! the load sit on the grid nodes along the chosen faces: planar side faces of the plate, or
//! cylinder faces whose axis runs through the thickness (holes, round edges). The system is solved
//! with the multigrid solver; stresses are taken at cell centres (von Mises, plane stress).
//!
//! Units: mm, N, MPa. Everything here is bounded and returns errors instead of panicking.

use crate::topo2d::{Solver, System, solve};

pub type P3 = [f64; 3];

fn comp(p: &P3, k: usize) -> f64 {
    p.get(k).copied().unwrap_or(0.0)
}

#[derive(Clone, Debug)]
pub enum FaceGeom {
    Plane { normal: P3 },
    Cylinder { axis: P3, origin: P3, radius: f64 },
    Other,
}

#[derive(Clone, Debug)]
pub struct Face {
    pub index: usize,
    pub name: String,
    pub centroid: P3,
    pub area: f64,
    pub geom: FaceGeom,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PlaneChoice {
    Auto,
    Xy,
    Xz,
    Yz,
}

pub struct Input<'a> {
    pub positions: &'a [P3],
    pub triangles: &'a [[u32; 3]],
    pub faces: &'a [Face],
    pub plane: PlaneChoice,
    /// Cells along the longer in-plane side.
    pub resolution: usize,
    pub fixed: &'a str,
    pub load: &'a str,
    /// Total force in N, component-local axes.
    pub force: P3,
    pub e_mpa: f64,
    pub nu: f64,
}

pub struct Output {
    pub axes: [usize; 3], // in-plane a, in-plane b, thickness w
    pub thickness: f64,
    pub cell: f64,
    pub nx: usize,
    pub ny: usize,
    pub solid_cells: usize,
    pub fixed_nodes: usize,
    pub load_nodes: usize,
    pub max_vm: f64,
    pub max_vm_at: P3,
    pub p99_vm: f64,
    pub max_disp: f64,
    pub max_disp_at: P3,
    pub cg_iterations: usize,
    pub converged: bool,
    pub warnings: Vec<String>,
    pub picture: Vec<String>,
}

const AXIS: [&str; 3] = ["x", "y", "z"];

/// Runs the check. `Err` carries a message for the user.
pub fn check(input: &Input) -> Result<Output, String> {
    if !(input.e_mpa.is_finite() && input.e_mpa > 0.0) || !(input.nu.is_finite() && (0.0..0.5).contains(&input.nu)) {
        return Err("the material needs a Young's modulus and a Poisson ratio from 0 to 0.5".into());
    }
    if !(10..=400).contains(&input.resolution) {
        return Err("resolution: expected 10 to 400 cells".into());
    }
    if input.positions.is_empty() || input.triangles.is_empty() {
        return Err("the body has no mesh".into());
    }
    if input.force.iter().any(|f| !f.is_finite()) {
        return Err("force: expected finite numbers".into());
    }
    let mut warnings = Vec::new();

    // Bounding box and axes.
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in input.positions {
        for k in 0..3 {
            let v = comp(p, k);
            if !v.is_finite() {
                return Err("the mesh has a non-finite point".into());
            }
            if let (Some(l), Some(h)) = (lo.get_mut(k), hi.get_mut(k)) {
                *l = l.min(v);
                *h = h.max(v);
            }
        }
    }
    let extent = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
    let w = match input.plane {
        PlaneChoice::Xy => 2,
        PlaneChoice::Xz => 1,
        PlaneChoice::Yz => 0,
        PlaneChoice::Auto => (0..3).min_by(|i, j| comp(&extent, *i).total_cmp(&comp(&extent, *j))).unwrap_or(2),
    };
    let (a, b) = match w {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    };
    let thickness = comp(&extent, w);
    let (la, lb) = (comp(&extent, a), comp(&extent, b));
    if !(thickness > 0.0 && la > 0.0 && lb > 0.0) {
        return Err("the body is flat or empty in some direction".into());
    }
    if thickness > 0.5 * la.min(lb) {
        warnings.push(format!(
            "the part is {thickness:.1} mm thick across {} but only {:.1} mm in its shorter side: it is not a thin plate, so a plane-stress result is rough",
            AXIS.get(w).unwrap_or(&"?"),
            la.min(lb)
        ));
    }
    let h = la.max(lb) / input.resolution as f64;
    let round8 = |n: usize| n.max(8).div_ceil(8) * 8;
    let nx = round8((la / h - 1e-9).ceil() as usize);
    let ny = round8((lb / h - 1e-9).ceil() as usize);
    if nx.saturating_mul(ny) > crate::topo2d::MAX_ELEMENTS {
        return Err("the grid would be too large; lower the resolution".into());
    }
    let (a0, b0) = (comp(&lo, a), comp(&lo, b));
    let wmid = comp(&lo, w) + 0.5 * thickness;

    // Rasterise: per cell, count mesh crossings above three heights (25%, 50%, 75% thickness).
    let levels = [comp(&lo, w) + 0.25 * thickness, wmid, comp(&lo, w) + 0.75 * thickness];
    let jitter = (h * 1.234_567e-4, h * 2.718_281e-4);
    let mut crossings = vec![[0u32; 3]; nx * ny];
    for tri in input.triangles {
        let [Some(p), Some(q), Some(r)] = tri.map(|i| input.positions.get(i as usize)) else {
            return Err("the mesh has a triangle with a bad vertex index".into());
        };
        let (pa, pb, pw) = (comp(p, a), comp(p, b), comp(p, w));
        let (qa, qb, qw) = (comp(q, a), comp(q, b), comp(q, w));
        let (ra, rb, rw) = (comp(r, a), comp(r, b), comp(r, w));
        let det = (qa - pa) * (rb - pb) - (ra - pa) * (qb - pb);
        if det.abs() < 1e-12 * la * lb {
            continue; // parallel to the thickness: rays along it never cross it
        }
        let ia0 = (((pa.min(qa).min(ra) - a0) / h - 0.5).floor().max(0.0)) as usize;
        let ia1 = ((((pa.max(qa).max(ra) - a0) / h - 0.5).ceil().max(0.0)) as usize).min(nx.saturating_sub(1));
        let ib0 = (((pb.min(qb).min(rb) - b0) / h - 0.5).floor().max(0.0)) as usize;
        let ib1 = ((((pb.max(qb).max(rb) - b0) / h - 0.5).ceil().max(0.0)) as usize).min(ny.saturating_sub(1));
        for i in ia0..=ia1 {
            for j in ib0..=ib1 {
                let ca = a0 + (i as f64 + 0.5) * h + jitter.0;
                let cb = b0 + (j as f64 + 0.5) * h + jitter.1;
                // Barycentric coordinates of the cell centre in the projected triangle.
                let s = ((ca - pa) * (rb - pb) - (ra - pa) * (cb - pb)) / det;
                let t = ((qa - pa) * (cb - pb) - (ca - pa) * (qb - pb)) / det;
                if s < 0.0 || t < 0.0 || s + t > 1.0 {
                    continue;
                }
                let hit = pw + s * (qw - pw) + t * (rw - pw);
                if let Some(c) = crossings.get_mut(i * ny + j) {
                    for (count, level) in c.iter_mut().zip(levels) {
                        if hit > level {
                            *count += 1;
                        }
                    }
                }
            }
        }
    }
    let solid: Vec<bool> = crossings.iter().map(|c| c[1] % 2 == 1).collect();
    let solid_cells = solid.iter().filter(|s| **s).count();
    if solid_cells == 0 {
        return Err("no material found in the body's plane; check the plane choice and that the body is closed".into());
    }
    let differ = crossings.iter().filter(|c| (c[0] % 2 == 1) != (c[1] % 2 == 1) || (c[2] % 2 == 1) != (c[1] % 2 == 1)).count();
    if differ as f64 > 0.02 * solid_cells as f64 {
        warnings.push(format!(
            "the outline changes through the thickness ({:.0}% of the cells differ between 25% and 75% depth): the part is not a flat plate of one thickness, so the result is approximate",
            100.0 * differ as f64 / solid_cells as f64
        ));
    }

    // Nodes that touch material.
    let nodes_y = ny + 1;
    let node_in_material = |i: usize, j: usize| {
        let cell = |ci: Option<usize>, cj: Option<usize>| match (ci, cj) {
            (Some(ci), Some(cj)) if ci < nx && cj < ny => solid.get(ci * ny + cj).copied().unwrap_or(false),
            _ => false,
        };
        cell(Some(i), Some(j)) || cell(i.checked_sub(1), Some(j)) || cell(Some(i), j.checked_sub(1)) || cell(i.checked_sub(1), j.checked_sub(1))
    };
    let node_pos = |i: usize, j: usize| (a0 + i as f64 * h, b0 + j as f64 * h);

    let select = |spec: &str, what: &str| -> Result<Vec<usize>, String> {
        let faces = resolve_faces(spec, input.faces, w, &lo, &hi, h, what)?;
        let mut nodes = Vec::new();
        for face in &faces {
            let before = nodes.len();
            match &face.geom {
                FaceGeom::Plane { normal } => {
                    let (na, nb) = (comp(normal, a), comp(normal, b));
                    let len = (na * na + nb * nb).sqrt();
                    if len < 0.95 {
                        return Err(format!(
                            "{what}: face {} is a top or bottom face of the plate; pick a side face (its edge in the outline), a hole, or a side like -x",
                            face.name
                        ));
                    }
                    let (na, nb) = (na / len, nb / len);
                    let (ca, cb) = (comp(&face.centroid, a), comp(&face.centroid, b));
                    let d = na * ca + nb * cb;
                    let s0 = -nb * ca + na * cb;
                    let half = 0.5 * face.area / thickness;
                    for i in 0..=nx {
                        for j in 0..=ny {
                            let (pa, pb) = node_pos(i, j);
                            if (na * pa + nb * pb - d).abs() <= 0.6 * h && (-nb * pa + na * pb - s0).abs() <= half + 0.6 * h && node_in_material(i, j) {
                                nodes.push(i * nodes_y + j);
                            }
                        }
                    }
                }
                FaceGeom::Cylinder { axis, origin, radius } => {
                    let norm = (axis.iter().map(|v| v * v).sum::<f64>()).sqrt();
                    if !(norm > 0.0 && (comp(axis, w) / norm).abs() > 0.95) {
                        return Err(format!("{what}: face {} is round but its axis does not run through the plate's thickness", face.name));
                    }
                    let (ca, cb) = (comp(origin, a), comp(origin, b));
                    for i in 0..=nx {
                        for j in 0..=ny {
                            let (pa, pb) = node_pos(i, j);
                            let dist = ((pa - ca).powi(2) + (pb - cb).powi(2)).sqrt();
                            if (dist - radius).abs() <= 0.75 * h && node_in_material(i, j) {
                                nodes.push(i * nodes_y + j);
                            }
                        }
                    }
                }
                FaceGeom::Other => return Err(format!("{what}: face {} is neither flat nor a cylinder; pick another face", face.name)),
            }
            if nodes.len() == before {
                return Err(format!("{what}: face {} touches no grid point at this resolution; raise the resolution", face.name));
            }
        }
        nodes.sort_unstable();
        nodes.dedup();
        Ok(nodes)
    };
    let fixed_nodes = select(input.fixed, "fixed")?;
    let load_nodes_all = select(input.load, "load")?;
    let fixed_set: std::collections::BTreeSet<usize> = fixed_nodes.iter().copied().collect();
    let load_nodes: Vec<usize> = load_nodes_all.iter().copied().filter(|n| !fixed_set.contains(n)).collect();
    if load_nodes.is_empty() {
        return Err("the load face lies entirely on the supports; nothing would move".into());
    }
    if load_nodes.len() < load_nodes_all.len() {
        warnings.push("part of the load face touches the supports; the load there goes straight into them".into());
    }
    let total = (input.force.iter().map(|f| f * f).sum::<f64>()).sqrt();
    if total == 0.0 {
        return Err("force: all components are zero".into());
    }
    let (fa, fb, fw) = (comp(&input.force, a), comp(&input.force, b), comp(&input.force, w));
    if fw.abs() > 0.01 * total {
        warnings.push(format!(
            "{:.1} N of the force points through the plate ({}); a plane-stress check ignores it (bending out of the plane is not covered)",
            fw,
            AXIS.get(w).unwrap_or(&"?")
        ));
    }
    if (fa * fa + fb * fb).sqrt() < 1e-9 * total {
        return Err("the force has no part in the plate's plane".into());
    }

    let ndof = 2 * (nx + 1) * nodes_y;
    let mut fixed = vec![false; ndof];
    for n in &fixed_nodes {
        for c in 0..2 {
            if let Some(slot) = fixed.get_mut(2 * n + c) {
                *slot = true;
            }
        }
    }
    let mut load = vec![0.0; ndof];
    let share = 1.0 / load_nodes.len() as f64;
    for n in &load_nodes {
        if let Some(slot) = load.get_mut(2 * n) {
            *slot += fa * share;
        }
        if let Some(slot) = load.get_mut(2 * n + 1) {
            *slot += fb * share;
        }
    }
    let et = input.e_mpa * thickness;
    let modulus: Vec<f64> = solid.iter().map(|s| if *s { et } else { et * 1e-9 }).collect();
    let system = System { nx, ny, nu: input.nu, modulus: &modulus, fixed, load: &load, solver: Solver::Multigrid, cg_tol: 1e-9, cg_max: 5000 };
    let analysis = solve(&system, None, &|| 0.0);
    if !analysis.converged {
        warnings.push("the solver did not fully converge; the supports may leave the part free to move or turn".into());
    }
    let u = &analysis.u;

    // Stresses at cell centres.
    let c = input.e_mpa / (1.0 - input.nu * input.nu);
    let corners = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];
    let mut vm = vec![f64::NAN; nx * ny];
    for i in 0..nx {
        for j in 0..ny {
            if !solid.get(i * ny + j).copied().unwrap_or(false) {
                continue;
            }
            let dofs = crate::mg2d::cell_dofs(ny, i, j);
            let ue = dofs.map(|d| u.get(d).copied().unwrap_or(0.0));
            let (mut ex, mut ey, mut gxy) = (0.0, 0.0, 0.0);
            for (n, (xi, eta)) in corners.iter().enumerate() {
                let (dx, dy) = (xi / (2.0 * h), eta / (2.0 * h));
                let ux = ue.get(2 * n).copied().unwrap_or(0.0);
                let uy = ue.get(2 * n + 1).copied().unwrap_or(0.0);
                ex += dx * ux;
                ey += dy * uy;
                gxy += dy * ux + dx * uy;
            }
            let sx = c * (ex + input.nu * ey);
            let sy = c * (ey + input.nu * ex);
            let txy = c * (1.0 - input.nu) / 2.0 * gxy;
            if let Some(slot) = vm.get_mut(i * ny + j) {
                *slot = (sx * sx - sx * sy + sy * sy + 3.0 * txy * txy).sqrt();
            }
        }
    }
    let to3 = |pa: f64, pb: f64| {
        let mut p = [0.0; 3];
        for (k, v) in [(a, pa), (b, pb), (w, wmid)] {
            if let Some(slot) = p.get_mut(k) {
                *slot = v;
            }
        }
        p
    };
    let mut max_vm = 0.0;
    let mut max_at = 0;
    for (k, v) in vm.iter().enumerate() {
        if v.is_finite() && *v > max_vm {
            max_vm = *v;
            max_at = k;
        }
    }
    let mut sorted: Vec<f64> = vm.iter().copied().filter(|v| v.is_finite()).collect();
    sorted.sort_by(f64::total_cmp);
    let p99 = sorted.get(((sorted.len() as f64 * 0.99) as usize).min(sorted.len().saturating_sub(1))).copied().unwrap_or(0.0);
    let (mi, mj) = (max_at / ny.max(1), max_at % ny.max(1));
    let max_vm_at = to3(a0 + (mi as f64 + 0.5) * h, b0 + (mj as f64 + 0.5) * h);

    let mut max_disp = 0.0;
    let mut max_disp_node = 0;
    for i in 0..=nx {
        for j in 0..=ny {
            if !node_in_material(i, j) {
                continue;
            }
            let n = i * nodes_y + j;
            let d = (u.get(2 * n).copied().unwrap_or(0.0).powi(2) + u.get(2 * n + 1).copied().unwrap_or(0.0).powi(2)).sqrt();
            if d > max_disp {
                max_disp = d;
                max_disp_node = n;
            }
        }
    }
    let (di, dj) = (max_disp_node / nodes_y, max_disp_node % nodes_y);
    let (dpa, dpb) = node_pos(di, dj);
    if !max_disp.is_finite() || max_disp > 0.2 * la.max(lb) {
        warnings.push("the part moves very far: check that the supports hold it (a single point or a part not connected to them lets it swing)".into());
    }
    if max_vm > 1.5 * p99 && sorted.len() >= 100 {
        warnings.push(
            "the peak stress sits in a few cells only (a sharp inside corner, or the edge of a support or load): such peaks grow with finer grids. Round the corner, or judge by the 99th percentile"
                .into(),
        );
    }

    Ok(Output {
        axes: [a, b, w],
        thickness,
        cell: h,
        nx,
        ny,
        solid_cells,
        fixed_nodes: fixed_nodes.len(),
        load_nodes: load_nodes.len(),
        max_vm,
        max_vm_at,
        p99_vm: p99,
        max_disp,
        max_disp_at: to3(dpa, dpb),
        cg_iterations: analysis.cg_iterations,
        converged: analysis.converged,
        warnings,
        picture: heat_picture(&vm, nx, ny, max_vm, 60),
    })
}

/// The faces a selector names: comma-separated face names, face numbers, or sides (+x, -y, ...).
fn resolve_faces<'f>(spec: &str, faces: &'f [Face], w: usize, lo: &P3, hi: &P3, h: f64, what: &str) -> Result<Vec<&'f Face>, String> {
    let mut out: Vec<&Face> = Vec::new();
    let parts: Vec<&str> = spec.split(',').map(str::trim).filter(|s| !s.is_empty()).collect();
    if parts.is_empty() {
        return Err(format!("{what}: name a face (as model.faces lists it, e.g. F1:-x), its number, or a side such as -x"));
    }
    if parts.len() > 32 {
        return Err(format!("{what}: at most 32 faces"));
    }
    for part in parts {
        let side = match part.as_bytes() {
            [sign @ (b'+' | b'-'), axis @ (b'x' | b'y' | b'z' | b'X' | b'Y' | b'Z')] => Some((*sign == b'+', (axis.to_ascii_lowercase() - b'x') as usize)),
            _ => None,
        };
        if let Some((plus, k)) = side {
            if k == w {
                return Err(format!("{what}: {part} is the plate's top or bottom; pick a side of the outline"));
            }
            let target = if plus { comp(hi, k) } else { comp(lo, k) };
            let sign = if plus { 1.0 } else { -1.0 };
            let found: Vec<&Face> = faces
                .iter()
                .filter(|f| match &f.geom {
                    FaceGeom::Plane { normal } => comp(normal, k) * sign > 0.99 && (comp(&f.centroid, k) - target).abs() <= 0.5 * h,
                    _ => false,
                })
                .collect();
            if found.is_empty() {
                return Err(format!("{what}: no flat face on the {part} side of the part"));
            }
            out.extend(found);
            continue;
        }
        let found = if let Ok(index) = part.parse::<usize>() {
            faces.iter().find(|f| f.index == index)
        } else {
            faces.iter().find(|f| f.name.eq_ignore_ascii_case(part))
        };
        match found {
            Some(f) => out.push(f),
            None => return Err(format!("{what}: no face {:?} on this body (model.faces lists them)", part.chars().take(64).collect::<String>())),
        }
    }
    Ok(out)
}

/// Text heat map: one character per block of cells (the block's highest stress), top row first.
fn heat_picture(vm: &[f64], nx: usize, ny: usize, max: f64, max_cols: usize) -> Vec<String> {
    let shades = ['.', ':', '-', '=', '+', '*', '#', '%', '@'];
    let step = nx.div_ceil(max_cols.max(1)).max(1);
    let ystep = 2 * step;
    let mut rows = Vec::new();
    let mut j = ny;
    while j > 0 {
        let j0 = j.saturating_sub(ystep);
        let mut line = String::new();
        let mut i = 0;
        while i < nx {
            let mut best: Option<f64> = None;
            for x in i..(i + step).min(nx) {
                for y in j0..j {
                    if let Some(v) = vm.get(x * ny + y).filter(|v| v.is_finite()) {
                        best = Some(best.map_or(*v, |b: f64| b.max(*v)));
                    }
                }
            }
            line.push(match best {
                None => ' ',
                Some(v) => {
                    let k = if max > 0.0 { ((v / max) * (shades.len() - 1) as f64).round() as usize } else { 0 };
                    shades.get(k.min(shades.len() - 1)).copied().unwrap_or('.')
                }
            });
            i += step;
        }
        rows.push(line.trim_end().to_string());
        j = j0;
    }
    while rows.last().is_some_and(String::is_empty) {
        rows.pop();
    }
    let first = rows.iter().position(|r| !r.is_empty()).unwrap_or(rows.len());
    rows.split_off(first)
}

/// A closed box mesh and its six faces, for tests and examples.
pub fn box_mesh(size: P3) -> (Vec<P3>, Vec<[u32; 3]>, Vec<Face>) {
    let [x, y, z] = size;
    let p = vec![[0.0, 0.0, 0.0], [x, 0.0, 0.0], [x, y, 0.0], [0.0, y, 0.0], [0.0, 0.0, z], [x, 0.0, z], [x, y, z], [0.0, y, z]];
    let t = vec![
        [0, 2, 1], [0, 3, 2], // bottom
        [4, 5, 6], [4, 6, 7], // top
        [0, 1, 5], [0, 5, 4], // -y
        [3, 7, 6], [3, 6, 2], // +y
        [0, 4, 7], [0, 7, 3], // -x
        [1, 2, 6], [1, 6, 5], // +x
    ];
    let face = |index: usize, name: &str, centroid: P3, area: f64, normal: P3| Face { index, name: name.into(), centroid, area, geom: FaceGeom::Plane { normal } };
    let faces = vec![
        face(0, "F1:+z", [x / 2.0, y / 2.0, z], x * y, [0.0, 0.0, 1.0]),
        face(1, "F1:-z", [x / 2.0, y / 2.0, 0.0], x * y, [0.0, 0.0, -1.0]),
        face(2, "F1:-y", [x / 2.0, 0.0, z / 2.0], x * z, [0.0, -1.0, 0.0]),
        face(3, "F1:+y", [x / 2.0, y, z / 2.0], x * z, [0.0, 1.0, 0.0]),
        face(4, "F1:-x", [0.0, y / 2.0, z / 2.0], y * z, [-1.0, 0.0, 0.0]),
        face(5, "F1:+x", [x, y / 2.0, z / 2.0], y * z, [1.0, 0.0, 0.0]),
    ];
    (p, t, faces)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
    use super::*;

    fn bar_input<'a>(mesh: &'a (Vec<P3>, Vec<[u32; 3]>, Vec<Face>), resolution: usize, force: P3) -> Input<'a> {
        Input {
            positions: &mesh.0,
            triangles: &mesh.1,
            faces: &mesh.2,
            plane: PlaneChoice::Auto,
            resolution,
            fixed: "-x",
            load: "+x",
            force,
            e_mpa: 68_900.0,
            nu: 0.33,
        }
    }

    /// Cantilever 100 x 10 mm, 5 mm thick, aluminium, 100 N down at the free end.
    /// Beam theory: bending stress at the root 6 F L / (t h^2) = 120 MPa; tip deflection
    /// F L^3 / (3 E I) + F L / (k G A) = 1.171 mm.
    #[test]
    fn cantilever_matches_beam_theory() {
        let mesh = box_mesh([100.0, 10.0, 5.0]);
        let out = check(&bar_input(&mesh, 200, [0.0, -100.0, 0.0])).unwrap();
        assert_eq!(out.axes, [0, 1, 2]);
        assert!((out.thickness - 5.0).abs() < 1e-12);
        assert!(out.converged);
        let (e, nu, l, h, t, f) = (68_900.0, 0.33, 100.0, 10.0, 5.0, 100.0);
        let i = t * h * h * h / 12.0;
        let g = e / (2.0 * (1.0 + nu));
        let tip = f * l * l * l / (3.0 * e * i) + f * l / (5.0 / 6.0 * g * h * t);
        let ratio_d = out.max_disp / tip;
        assert!((0.95..=1.02).contains(&ratio_d), "tip {} vs {tip}", out.max_disp);
        // The 99th percentile sits in the root's outer fibres, away from the corner peak.
        let sigma = 6.0 * f * l / (t * h * h);
        let ratio_s = out.p99_vm / sigma;
        assert!((0.85..=1.1).contains(&ratio_s), "p99 {} vs {sigma}", out.p99_vm);
        assert!(out.max_vm_at[0] < 5.0, "peak should be at the fixed end: {:?}", out.max_vm_at);
        assert_eq!(out.solid_cells, 200 * 20);
    }

    #[test]
    fn plane_is_the_thinnest_direction() {
        let mesh = box_mesh([100.0, 4.0, 20.0]);
        let mut input = bar_input(&mesh, 100, [0.0, 0.0, -50.0]);
        input.fixed = "-x";
        input.load = "+x";
        let out = check(&input).unwrap();
        assert_eq!(out.axes, [0, 2, 1]);
    }

    #[test]
    fn out_of_plane_force_is_warned_about() {
        let mesh = box_mesh([100.0, 10.0, 5.0]);
        let out = check(&bar_input(&mesh, 100, [0.0, -100.0, 30.0])).unwrap();
        assert!(out.warnings.iter().any(|w| w.contains("through the plate")));
    }

    #[test]
    fn bad_selections_are_explained() {
        let mesh = box_mesh([100.0, 10.0, 5.0]);
        for (fixed, load, needle) in [
            ("+z", "+x", "top or bottom"),
            ("F1:+z", "+x", "top or bottom"),
            ("nope", "+x", "no face"),
            ("", "+x", "name a face"),
            ("-x", "-x", "entirely on the supports"),
        ] {
            let mut input = bar_input(&mesh, 100, [0.0, -100.0, 0.0]);
            input.fixed = fixed;
            input.load = load;
            let err = check(&input).err().unwrap();
            assert!(err.contains(needle), "{fixed}/{load}: {err}");
        }
        let mut input = bar_input(&mesh, 100, [0.0, 0.0, 0.0]);
        input.load = "+x";
        assert!(check(&input).is_err());
    }

    #[test]
    fn faces_by_number_and_name_work() {
        let mesh = box_mesh([100.0, 10.0, 5.0]);
        let mut input = bar_input(&mesh, 100, [0.0, -100.0, 0.0]);
        input.fixed = "4";
        input.load = "f1:+X";
        assert!(check(&input).is_ok());
    }

    #[test]
    fn hole_support_selects_a_ring_of_nodes() {
        let mut mesh = box_mesh([100.0, 30.0, 5.0]);
        // A cylinder face at (15, 15), radius 3. The mesh has no hole, so the ring lies in the
        // material: enough to check the selection.
        mesh.2.push(Face {
            index: 9,
            name: "F2:bore:side:0".into(),
            centroid: [15.0, 15.0, 2.5],
            area: 94.0,
            geom: FaceGeom::Cylinder { axis: [0.0, 0.0, 1.0], origin: [15.0, 15.0, -200.0], radius: 3.0 },
        });
        let mut input = bar_input(&mesh, 100, [0.0, -100.0, 0.0]);
        input.fixed = "F2:bore:side:0";
        let out = check(&input).unwrap();
        assert!(out.fixed_nodes >= 8, "{}", out.fixed_nodes);
        assert!(out.converged);
    }

    #[test]
    fn hostile_inputs_never_panic() {
        let mesh = box_mesh([100.0, 10.0, 5.0]);
        let bad_tri = vec![[0u32, 1, 99]];
        let nan_pos = vec![[f64::NAN, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let mut inputs = Vec::new();
        let base = bar_input(&mesh, 100, [0.0, -100.0, 0.0]);
        inputs.push(Input { triangles: &bad_tri, ..bar_input(&mesh, 100, [0.0, -1.0, 0.0]) });
        inputs.push(Input { positions: &nan_pos, triangles: &bad_tri, ..bar_input(&mesh, 100, [0.0, -1.0, 0.0]) });
        inputs.push(Input { resolution: 0, ..bar_input(&mesh, 100, [0.0, -1.0, 0.0]) });
        inputs.push(Input { e_mpa: f64::NAN, ..bar_input(&mesh, 100, [0.0, -1.0, 0.0]) });
        inputs.push(Input { force: [f64::INFINITY, 0.0, 0.0], ..bar_input(&mesh, 100, [0.0, -1.0, 0.0]) });
        inputs.push(Input { positions: &[], ..bar_input(&mesh, 100, [0.0, -1.0, 0.0]) });
        for input in &inputs {
            assert!(check(input).is_err());
        }
        assert!(check(&base).is_ok());
    }
}
