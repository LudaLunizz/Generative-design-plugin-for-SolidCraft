//! 3D voxel stress check of an STL part, outside SolidCraft (no plugin limits).
//!
//! stress3d PART.stl --material ID --fix SIDE --pin AXIS,C1,C2,R --case NAME:FX,FY,FZ[:MX,MY,MZ] ...
//!          [--voxel MM] [--sf N] [--out DIR] [--jacobi]
//!
//!   --fix SIDE       supports: the part's xmin/xmax/ymin/ymax/zmin/zmax face (material within half
//!                    a voxel of it is held)
//!   --pin A,C1,C2,R  load region: the bore along axis A (x, y or z) through the other two
//!                    coordinates C1, C2, radius R: material nodes from R - h/2 to R + h from its axis
//!   --case           a load case: force (N) spread over the pin, optional couple (N mm) about it
//!   --material       a library id (materials/library.json); --sf the safety factor needed (2)
//!   --jacobi         solve with plain Jacobi-preconditioned CG instead of multigrid (for comparison)
//!   --out DIR        writes report.json and per case three heat-map PNGs (largest stress along x, y, z)

#![allow(clippy::indexing_slicing)] // command-line tool: every index below follows a length check (numbers() returns exactly n values)

use generative_design_engine::library::Library;
use generative_design_engine::stress3d::{Case, Region, Setup, analyse, png, projection, read_stl, voxelize};
use serde_json::json;
use std::time::Instant;

fn main() {
    if let Err(e) = run() {
        eprintln!("stress3d: {e}");
        std::process::exit(2);
    }
}

fn numbers(s: &str, n: usize, what: &str) -> Result<Vec<f64>, String> {
    let v: Vec<f64> = s.split(',').map(|x| x.trim().parse::<f64>().map_err(|_| format!("{what}: {x:?} is not a number"))).collect::<Result<_, _>>()?;
    if v.len() != n || v.iter().any(|x| !x.is_finite()) {
        return Err(format!("{what}: expected {n} finite numbers"));
    }
    Ok(v)
}

fn axis(s: &str) -> Result<usize, String> {
    match s {
        "x" => Ok(0),
        "y" => Ok(1),
        "z" => Ok(2),
        _ => Err(format!("axis {s:?}: expected x, y or z")),
    }
}

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut stl = None;
    let (mut voxel, mut material, mut fix, mut pin, mut out, mut sf) = (2.5, String::from("petg-printed"), None, None, None, 2.0);
    let mut cases = Vec::new();
    let mut multigrid = true;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let mut value = |name: &str| it.next().cloned().ok_or(format!("{name} needs a value"));
        match a.as_str() {
            "--jacobi" => multigrid = false,
            "--voxel" => voxel = value("--voxel")?.parse().map_err(|_| "--voxel: a number in mm")?,
            "--material" => material = value("--material")?,
            "--fix" => fix = Some(value("--fix")?),
            "--pin" => pin = Some(value("--pin")?),
            "--out" => out = Some(value("--out")?),
            "--sf" => sf = value("--sf")?.parse().map_err(|_| "--sf: a number")?,
            "--case" => {
                let v = value("--case")?;
                let parts: Vec<&str> = v.split(':').collect();
                let (name, force, moment) = match parts.as_slice() {
                    [n, f] => (n.to_string(), numbers(f, 3, "force")?, vec![0.0; 3]),
                    [n, f, m] => (n.to_string(), numbers(f, 3, "force")?, numbers(m, 3, "moment")?),
                    _ => return Err("--case NAME:FX,FY,FZ[:MX,MY,MZ]".into()),
                };
                cases.push(Case { name, force: [force[0], force[1], force[2]], moment: [moment[0], moment[1], moment[2]] });
            }
            s if s.starts_with("--") => return Err(format!("unknown option {s}")),
            s => stl = Some(s.to_string()),
        }
    }
    let stl = stl.ok_or("give an STL file")?;
    if cases.is_empty() {
        return Err("give at least one --case".into());
    }
    if !(0.2..=50.0).contains(&voxel) {
        return Err("--voxel: 0.2 to 50 mm".into());
    }
    let library = Library::builtin()?;
    let m = library.find(&material).ok_or(format!("material {material:?} is not in the library"))?;
    let s = m.strengths();
    let (e_gpa, nu) = (s.young_modulus_gpa.ok_or("material has no stiffness")?, s.poisson_ratio.ok_or("material has no Poisson ratio")?);
    let (strength, basis) = match (s.yield_mpa, s.tensile_mpa) {
        (Some(y), _) => (y, "yield"),
        (None, Some(t)) => (t, "tensile strength"),
        _ => return Err("material has no strength".into()),
    };
    let fix = fix.ok_or("give --fix")?;
    let (fix_axis, fix_max) = match fix.as_str() {
        "xmin" => (0, false),
        "xmax" => (0, true),
        "ymin" => (1, false),
        "ymax" => (1, true),
        "zmin" => (2, false),
        "zmax" => (2, true),
        _ => return Err("--fix: xmin, xmax, ymin, ymax, zmin or zmax".into()),
    };
    let pin = pin.ok_or("give --pin")?;
    let pp: Vec<&str> = pin.split(',').collect();
    let [a, c1, c2, r] = pp.as_slice() else { return Err("--pin AXIS,C1,C2,R".into()) };
    let pa = axis(a)?;
    let nums = numbers(&format!("{c1},{c2},{r}"), 3, "--pin")?;
    let mut point = [0.0; 3];
    let others: Vec<usize> = (0..3).filter(|k| *k != pa).collect();
    for (k, v) in others.iter().zip(&nums) {
        point[*k] = *v;
    }
    let mut dir = [0.0; 3];
    dir[pa] = 1.0;
    let radius = nums[2];

    let t0 = Instant::now();
    let bytes = std::fs::read(&stl).map_err(|e| format!("{stl}: {e}"))?;
    let (positions, triangles) = read_stl(&bytes)?;
    let (grid, solid) = voxelize(&positions, &triangles, voxel)?;
    let material_cells = solid.iter().filter(|v| **v).count();
    eprintln!("{stl}: {} triangles, {}x{}x{} voxels of {voxel} mm, {material_cells} with material ({:.0} cm3)", triangles.len(), grid.n[0], grid.n[1], grid.n[2], material_cells as f64 * voxel.powi(3) / 1000.0);
    let fixed = Region::Side { axis: fix_axis, max: fix_max, depth: 0.5 * voxel };
    let load = Region::Ring { point, dir, r_in: radius - 0.5 * voxel, r_out: radius + voxel };
    let setup = Setup { grid: &grid, solid: &solid, e_mpa: e_gpa * 1000.0, nu, fixed: &fixed, load: &load, cg_tol: 1e-8, cg_max: 50_000, exclude_layers: 2, cell_scale: None, multigrid };
    let result = analyse(&setup, &cases, |msg| eprintln!("  {msg} ..."))?;
    let secs = t0.elapsed().as_secs_f64();

    println!("{} with {} (E {e_gpa} GPa, nu {nu}, {basis} {strength} MPa), voxels {voxel} mm, {:.1} s", stl, m.name, secs);
    println!("held: {fix} ({} nodes), load: pin {pin} ({} nodes){}", result.fixed_nodes, result.load_nodes, if result.cells_dropped > 0 { format!(", {} loose voxels left out", result.cells_dropped) } else { String::new() });
    println!("{:<8} {:>9} {:>10} {:>9} {:>7} {:>10} {:>6}  away peak at (mm) / overall peak at (mm)", "case", "peak MPa", "away MPa", "p99 MPa", "SF", "max d mm", "CG");
    let mut cases_json = Vec::new();
    for c in &result.cases {
        // The safety factor uses the peak away from the load and support zones (see the doc).
        let found = if c.max_vm_away > 0.0 { strength / c.max_vm_away } else { f64::INFINITY };
        println!(
            "{:<8} {:>9.3} {:>10.3} {:>9.3} {:>7.1} {:>10.4} {:>6}  ({:.1}, {:.1}, {:.1}) / ({:.1}, {:.1}, {:.1}){}",
            c.name, c.max_vm, c.max_vm_away, c.p99_vm, found, c.max_disp, c.cg_iterations,
            c.max_vm_away_at[0], c.max_vm_away_at[1], c.max_vm_away_at[2], c.max_vm_at[0], c.max_vm_at[1], c.max_vm_at[2],
            if found >= sf { "" } else { "  BELOW the needed safety factor" }
        );
        cases_json.push(json!({
            "name": c.name, "peak_von_mises_mpa": c.max_vm, "peak_at_mm": c.max_vm_at, "peak_away_from_load_and_supports_mpa": c.max_vm_away,
            "peak_away_at_mm": c.max_vm_away_at, "p99_von_mises_mpa": c.p99_vm,
            "safety_factor": found, "passes": found >= sf, "max_displacement_mm": c.max_disp, "max_displacement_at_mm": c.max_disp_at,
            "cg_iterations": c.cg_iterations, "converged": c.converged,
        }));
        if let Some(dir) = &out {
            std::fs::create_dir_all(dir).map_err(|e| format!("{dir}: {e}"))?;
            for (ax, name) in [(0, "along-x"), (1, "along-y"), (2, "along-z")] {
                let (w, h, rgb) = projection(&grid, &c.vm, ax, c.max_vm_away, 4); // colours up to the peak away from the load; above it shows red
                let path = format!("{dir}/{}-{name}.png", c.name);
                std::fs::write(&path, png(w, h, &rgb)).map_err(|e| format!("{path}: {e}"))?;
            }
        }
    }
    for w in &result.warnings {
        println!("warning: {w}");
    }
    if let Some(dir) = &out {
        let report = json!({
            "stl": stl, "triangles": triangles.len(), "voxel_mm": voxel, "grid": grid.n, "material_voxels": result.cells_material,
            "dropped_voxels": result.cells_dropped, "material": m.name, "e_gpa": e_gpa, "nu": nu, "strength_mpa": strength, "strength_basis": basis,
            "safety_factor_needed": sf, "fixed": fix, "fixed_nodes": result.fixed_nodes, "pin": pin, "load_nodes": result.load_nodes,
            "cases": cases_json, "warnings": result.warnings, "seconds": secs,
        });
        let path = format!("{dir}/report.json");
        std::fs::write(&path, serde_json::to_string_pretty(&report).unwrap_or_default()).map_err(|e| format!("{path}: {e}"))?;
        println!("written: {dir}/report.json and {} heat maps", 3 * result.cases.len());
    }
    Ok(())
}
