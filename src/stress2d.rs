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
    const SS: usize = 4;
    let (snx, sny) = (nx * SS, ny * SS);
    let mut sub = vec![0u8; snx * sny];
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
        // Sub-samples (SS x SS per cell): parity of the crossings above the mid-thickness level.
        {
            let (sa0, sa1) = (ia0 * SS, ((ia1 + 1) * SS).min(snx));
            let (sb0, sb1) = (ib0 * SS, ((ib1 + 1) * SS).min(sny));
            let sh = h / SS as f64;
            for si in sa0..sa1 {
                let ca = a0 + (si as f64 + 0.5) * sh + jitter.0;
                for sj in sb0..sb1 {
                    let cb = b0 + (sj as f64 + 0.5) * sh + jitter.1;
                    let s = ((ca - pa) * (rb - pb) - (ra - pa) * (cb - pb)) / det;
                    let t = ((qa - pa) * (cb - pb) - (ca - pa) * (qb - pb)) / det;
                    if s < 0.0 || t < 0.0 || s + t > 1.0 {
                        continue;
                    }
                    if pw + s * (qw - pw) + t * (rw - pw) > wmid {
                        if let Some(c) = sub.get_mut(si * sny + sj) {
                            *c ^= 1;
                        }
                    }
                }
            }
        }
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
    let mut frac = vec![0.0f64; nx * ny];
    for i in 0..nx {
        for j in 0..ny {
            let mut n = 0u32;
            for di in 0..SS {
                for dj in 0..SS {
                    n += u32::from(sub.get((i * SS + di) * sny + j * SS + dj).copied().unwrap_or(0));
                }
            }
            if let Some(f) = frac.get_mut(i * ny + j) {
                *f = f64::from(n) / (SS * SS) as f64;
            }
        }
    }
    // A cell counts as material when at least half of it is; nodes touch material when any part is.
    let solid: Vec<bool> = frac.iter().map(|f| *f >= 0.5).collect();
    let present: Vec<bool> = frac.iter().map(|f| *f > 0.0).collect();
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
            (Some(ci), Some(cj)) if ci < nx && cj < ny => present.get(ci * ny + cj).copied().unwrap_or(false),
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
    // Each cell is as stiff as the share of it that is material, so a round edge is smooth, not stepped.
    let modulus: Vec<f64> = frac.iter().map(|f| et * f.max(1e-9)).collect();
    let system = System { nx, ny, nu: input.nu, modulus: &modulus, fixed, load: &load, solver: Solver::Multigrid, cg_tol: 1e-9, cg_max: 5000, coarse_budget: 0 };
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
    let (surface_vm, surface_at) = surface_peak(&frac, &vm, nx, ny, 1);
    let (surface_wide, _) = surface_peak(&frac, &vm, nx, ny, 2);
    let mut max_vm = 0.0;
    let mut max_at = 0;
    for (k, v) in vm.iter().enumerate() {
        if v.is_finite() && *v > max_vm {
            max_vm = *v;
            max_at = k;
        }
    }
    let cell_max = max_vm;
    if surface_vm > 0.0 {
        max_vm = surface_vm;
        max_at = surface_at;
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
    if surface_vm > 0.0 && (surface_vm - surface_wide) > 0.04 * surface_vm {
        warnings.push(format!(
            "the stress changes quickly over the few cells next to the peak (a small hole or round edge for this grid): the peak of {surface_vm:.1} MPa is probably low, by more than {:.0}%; raise the resolution",
            100.0 * (surface_vm - surface_wide) / surface_vm
        ));
    }
    if max_vm > 1.5 * p99 && sorted.len() >= 100 {
        warnings.push(
            "the peak stress sits in a few cells only. At a round hole or fillet it is the stress at the surface and settles with finer grids; at a sharp inside corner or the edge of a support or load it keeps growing with them: round the corner, or judge by the 99th percentile"
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
        picture: heat_picture(&vm, nx, ny, cell_max, 60),
    })
}

/// Area (in cells) of the square [-r, r]^2 on the side `p . n <= d` of a line.
fn clipped_area(n: (f64, f64), d: f64, r: f64) -> f64 {
    // Polygon clipping of the square against the half plane (Sutherland-Hodgman, one edge).
    let pts = [(-r, -r), (r, -r), (r, r), (-r, r)];
    let side = |p: (f64, f64)| p.0 * n.0 + p.1 * n.1 - d;
    let mut poly: Vec<(f64, f64)> = Vec::with_capacity(6);
    for (p, q) in pts.iter().zip(pts.iter().cycle().skip(1)) {
        let (sp, sq) = (side(*p), side(*q));
        if sp <= 0.0 {
            poly.push(*p);
        }
        if (sp <= 0.0) != (sq <= 0.0) {
            let t = sp / (sp - sq);
            poly.push((p.0 + t * (q.0 - p.0), p.1 + t * (q.1 - p.1)));
        }
    }
    let mut area = 0.0;
    for (p, q) in poly.iter().zip(poly.iter().cycle().skip(1)) {
        area += p.0 * q.1 - q.0 * p.1;
    }
    0.5 * area.abs()
}

/// The stress at the surface of the part, not at the centre of the nearest cell: the stress
/// changes fast next to a hole or a fillet (about 20 MPa per mm for a 5 mm hole), so a cell
/// centre 0 to 0.3 mm inside the surface reads several MPa low, and by a different amount for
/// every grid. For every full cell near an edge: a plane is fitted to the stress of the full cells
/// of its (2r+1)^2 patch, the surface is placed from the share of material in the patch (a straight
/// edge of that area, along the gradient of the share), and the plane is read there.
/// Returns (peak, cell index of the peak); cells away from edges count with their own value.
fn surface_peak(frac: &[f64], vm: &[f64], nx: usize, ny: usize, r: isize) -> (f64, usize) {
    let at = |i: isize, j: isize| -> (f64, f64) {
        if i < 0 || j < 0 || i >= nx as isize || j >= ny as isize {
            return (0.0, f64::NAN);
        }
        let k = i as usize * ny + j as usize;
        (frac.get(k).copied().unwrap_or(0.0), vm.get(k).copied().unwrap_or(f64::NAN))
    };
    let (mut best, mut best_at) = (0.0, 0);
    let half = r as f64 + 0.5;
    let full_area = (2.0 * half) * (2.0 * half);
    for i in 0..nx as isize {
        for j in 0..ny as isize {
            let (f0, v0) = at(i, j);
            if f0 < 0.999 || !v0.is_finite() {
                continue;
            }
            let mut value = v0;
            let (mut area, mut gx, mut gy) = (0.0, 0.0, 0.0);
            for di in -r..=r {
                for dj in -r..=r {
                    let f = at(i + di, j + dj).0;
                    area += f;
                    gx += di as f64 * f;
                    gy += dj as f64 * f;
                }
            }
            if area < full_area - 1e-9 {
                let glen = (gx * gx + gy * gy).sqrt();
                if glen > 1e-6 {
                    let n = (-gx / glen, -gy / glen); // out of the material
                    let (mut lo_d, mut hi_d) = (-half * 1.5, half * 1.5);
                    for _ in 0..40 {
                        let mid = 0.5 * (lo_d + hi_d);
                        if clipped_area(n, mid, half) < area {
                            lo_d = mid;
                        } else {
                            hi_d = mid;
                        }
                    }
                    let d = 0.5 * (lo_d + hi_d); // distance from the cell centre to the surface, in cells
                    // Least-squares plane through the full cells of the patch: v = a + b x + c y.
                    let (mut s1, mut sx, mut sy, mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                    let (mut sv, mut sxv, mut syv) = (0.0, 0.0, 0.0);
                    for di in -r..=r {
                        for dj in -r..=r {
                            let (f, v) = at(i + di, j + dj);
                            if f < 0.999 || !v.is_finite() {
                                continue;
                            }
                            let (x, y) = (di as f64, dj as f64);
                            s1 += 1.0;
                            sx += x;
                            sy += y;
                            sxx += x * x;
                            sxy += x * y;
                            syy += y * y;
                            sv += v;
                            sxv += x * v;
                            syv += y * v;
                        }
                    }
                    // Solve the 3x3 normal equations by Cramer's rule.
                    let det = s1 * (sxx * syy - sxy * sxy) - sx * (sx * syy - sxy * sy) + sy * (sx * sxy - sxx * sy);
                    if s1 >= 4.0 && det.abs() > 1e-9 && (0.0..=2.0 * half).contains(&d) {
                        let a = (sv * (sxx * syy - sxy * sxy) - sx * (sxv * syy - sxy * syv) + sy * (sxv * sxy - sxx * syv)) / det;
                        let b = (s1 * (sxv * syy - sxy * syv) - sv * (sx * syy - sxy * sy) + sy * (sx * syv - sxv * sy)) / det;
                        let c = (s1 * (sxx * syv - sxv * sxy) - sx * (sx * syv - sxv * sy) + sv * (sx * sxy - sxx * sy)) / det;
                        let extrapolated = a + (b * n.0 + c * n.1) * d;
                        if extrapolated.is_finite() {
                            value = extrapolated.max(v0);
                        }
                    }
                }
            }
            if value > best {
                best = value;
                best_at = i as usize * ny + j as usize;
            }
        }
    }
    (best, best_at)
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

/// A closed plate (size x, y, thickness z, corner at the origin) with one round hole through its
/// thickness, for tests and examples. Only the top and bottom faces matter for the grid check; the
/// side walls are left out of the mesh (they are parallel to the thickness and never sampled).
pub fn plate_with_hole_mesh(size: P3, centre: [f64; 2], radius: f64, segments: usize) -> (Vec<P3>, Vec<[u32; 3]>) {
    let [x, y, z] = size;
    let n = segments.max(8);
    // Outer boundary points along the same rays as the hole points, plus the four corners.
    let ray_hit = |ang: f64| -> [f64; 2] {
        let (dx, dy) = (ang.cos(), ang.sin());
        let mut t = f64::INFINITY;
        if dx > 1e-12 {
            t = t.min((x - centre[0]) / dx);
        } else if dx < -1e-12 {
            t = t.min(-centre[0] / dx);
        }
        if dy > 1e-12 {
            t = t.min((y - centre[1]) / dy);
        } else if dy < -1e-12 {
            t = t.min(-centre[1] / dy);
        }
        [centre[0] + t * dx, centre[1] + t * dy]
    };
    let mut angles: Vec<f64> = (0..n).map(|k| std::f64::consts::TAU * (k as f64 + 0.37) / n as f64).collect();
    for corner in [[x, y], [0.0, y], [0.0, 0.0], [x, 0.0]] {
        let mut a = (corner[1] - centre[1]).atan2(corner[0] - centre[0]);
        if a < 0.0 {
            a += std::f64::consts::TAU;
        }
        angles.push(a);
    }
    angles.sort_by(f64::total_cmp);
    // One ring of inner (hole) and outer points per angle; every ring step becomes two triangles.
    let m = angles.len();
    let mut pos = Vec::new();
    for zz in [0.0, z] {
        for a in &angles {
            pos.push([centre[0] + radius * a.cos(), centre[1] + radius * a.sin(), zz]);
        }
        for a in &angles {
            let o = ray_hit(*a);
            pos.push([o[0], o[1], zz]);
        }
    }
    let mut tri = Vec::new();
    for layer in 0..2usize {
        let base = layer * 2 * m;
        for k in 0..m {
            let k2 = (k + 1) % m;
            let (i0, i1, o0, o1) = ((base + k) as u32, (base + k2) as u32, (base + m + k) as u32, (base + m + k2) as u32);
            if layer == 0 {
                tri.push([i0, o1, o0]);
                tri.push([i0, i1, o1]);
            } else {
                tri.push([i0, o0, o1]);
                tri.push([i0, o1, i1]);
            }
        }
    }
    (pos, tri)
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

    /// Plate 120 x 40 x 5 with a 10 mm hole, 3000 N pull: net stress 20 MPa, Peterson's Kt on the net
    /// section for d/w = 0.25 is 2.42, so the stress at the hole's side is about 48.4 MPa. The peak at a
    /// round edge must not jump by several MPa when the hole moves a fraction of a cell on the grid.
    #[test]
    fn hole_peak_is_stable_and_near_theory() {
        let (_, _, faces) = box_mesh([120.0, 40.0, 5.0]);
        let mut peaks = Vec::new();
        for shift in [0.0, 0.07, 0.13, 0.21] {
            let (pos, tri) = plate_with_hole_mesh([120.0, 40.0, 5.0], [60.0 + shift, 20.0 + 0.6 * shift], 5.0, 96);
            let input = Input { positions: &pos, triangles: &tri, faces: &faces, plane: PlaneChoice::Auto, resolution: 300, fixed: "-x", load: "+x", force: [3000.0, 0.0, 0.0], e_mpa: 68_900.0, nu: 0.33 };
            let out = check(&input).unwrap();
            assert!(out.converged);
            assert!((out.p99_vm - 22.0).abs() < 1.5, "p99 {}", out.p99_vm);
            assert!((out.max_vm_at[1] - 20.0).abs() > 3.0 && (out.max_vm_at[0] - 60.0).abs() < 3.0, "peak should be at the hole's side: {:?}", out.max_vm_at);
            peaks.push(out.max_vm);
        }
        let (lo, hi) = peaks.iter().fold((f64::MAX, f64::MIN), |(l, h), v| (l.min(*v), h.max(*v)));
        assert!(lo > 0.90 * 48.4 && hi < 1.04 * 48.4, "peaks {peaks:?}");
        assert!(hi - lo < 0.06 * hi, "peak jumps with the hole position: {peaks:?}");
    }

    #[test]
    fn coarse_grid_at_a_hole_warns_that_the_peak_is_low() {
        let (_, _, faces) = box_mesh([120.0, 40.0, 5.0]);
        let (pos, tri) = plate_with_hole_mesh([120.0, 40.0, 5.0], [60.0, 20.0], 5.0, 96);
        let run = |resolution| {
            let input = Input { positions: &pos, triangles: &tri, faces: &faces, plane: PlaneChoice::Auto, resolution, fixed: "-x", load: "+x", force: [3000.0, 0.0, 0.0], e_mpa: 68_900.0, nu: 0.33 };
            check(&input).unwrap()
        };
        let coarse = run(100);
        assert!(coarse.warnings.iter().any(|w| w.contains("probably low")), "{:?}", coarse.warnings);
        let fine = run(400);
        assert!(!fine.warnings.iter().any(|w| w.contains("probably low")), "{:?}", fine.warnings);
        assert!(fine.max_vm > coarse.max_vm);
    }

    /// A straight edge on a cell boundary must not be moved or sharpened by the surface reading.
    #[test]
    fn surface_peak_leaves_flat_fields_alone() {
        let (nx, ny) = (12, 10);
        let frac = vec![1.0; nx * ny];
        let vm = vec![7.0; nx * ny];
        let (peak, _) = surface_peak(&frac, &vm, nx, ny, 1);
        assert!((peak - 7.0).abs() < 1e-9);
        // A half-covered edge column with a stress rising 2 MPa per cell towards the edge: the surface
        // sits 0.5 cell beyond the last full column, so the reading is 1 MPa above that column.
        let mut frac = vec![1.0; nx * ny];
        let mut vm = vec![0.0; nx * ny];
        for i in 0..nx {
            for j in 0..ny {
                vm[i * ny + j] = 10.0 + 2.0 * i as f64;
                if i == nx - 1 {
                    frac[i * ny + j] = 0.0;
                }
            }
        }
        let (peak, _) = surface_peak(&frac, &vm, nx, ny, 1);
        // Last full column i = nx-2 has vm 10 + 2 (nx-2) = 30; its surface is 0.5 cell on: 31.
        assert!((peak - 31.0).abs() < 0.3, "peak {peak}");
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
