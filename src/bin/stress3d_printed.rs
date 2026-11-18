//! 3D voxel stress check of a PRINTED part: solid walls and top/bottom skins around a lighter infill
//! core, with the layers stacked along z (the part stands on its zmin face when printed).
//!
//! stress3d_printed PART.stl --material ID --fix SIDE --pin AXIS,C1,C2,R --case NAME:FX,FY,FZ[:MX,MY,MZ] ...
//!          [--voxel MM] [--walls 8] [--wall-first 0.42] [--wall-step 0.45] [--top-layers 8]
//!          [--bottom-layers 3] [--layer 0.2] [--infill 0.25] [--infill-exp 1] [--z-tensile 30]
//!          [--z-long 10] [--density 1.24] [--out DIR]
//!
//! The walls and skins are the full material (stiffness 1); the infill core is a homogenised cell whose
//! stiffness is infill^infill-exp (exp 1 = stretch-dominated, the upper bound; exp 2 = bending-dominated,
//! the lower bound). Stresses in the walls are the material's own; the stress in an infill strut is
//! estimated as the cell average divided by the infill fraction.

#![allow(clippy::indexing_slicing)] // command-line tool: every index below follows a length check

use generative_design_engine::library::Library;
use generative_design_engine::stress3d::{Case, Region, Setup, analyse, png, projection, read_stl, voxelize};
use serde_json::json;
use std::time::Instant;

fn main() {
    if let Err(e) = run() {
        eprintln!("stress3d_printed: {e}");
        std::process::exit(2);
    }
}

fn numbers(s: &str, n: usize, what: &str) -> Result<Vec<f64>, String> {
    let v: Result<Vec<f64>, _> = s.split(',').map(|x| x.trim().parse::<f64>()).collect();
    let v = v.map_err(|_| format!("{what}: numbers separated by commas"))?;
    if v.len() != n {
        return Err(format!("{what}: {n} numbers"));
    }
    Ok(v)
}

fn pn(v: String, name: &str) -> Result<f64, String> {
    v.parse().map_err(|_| format!("{name}: a number"))
}

fn axis(s: &str) -> Result<usize, String> {
    match s {
        "x" => Ok(0),
        "y" => Ok(1),
        "z" => Ok(2),
        _ => Err("axis: x, y or z".into()),
    }
}

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut stl = None;
    let (mut voxel, mut material, mut fix, mut pin, mut out) = (1.5, String::from("pla-printed"), None, None, None);
    let (mut walls, mut wall_first, mut wall_step, mut top_layers, mut bottom_layers, mut layer) = (8.0, 0.42, 0.45, 8.0, 3.0, 0.2);
    let (mut infill, mut infill_exp, mut z_tensile, mut z_long, mut density) = (0.25, 1.0, 30.0, 10.0, 1.24);
    let mut mass_only = false;
    let mut cases = Vec::new();
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let mut value = |name: &str| it.next().cloned().ok_or(format!("{name} needs a value"));
        match a.as_str() {
            "--voxel" => voxel = pn(value("--voxel")?, "--voxel")?,
            "--material" => material = value("--material")?,
            "--fix" => fix = Some(value("--fix")?),
            "--pin" => pin = Some(value("--pin")?),
            "--out" => out = Some(value("--out")?),
            "--walls" => walls = pn(value("--walls")?, "--walls")?,
            "--wall-first" => wall_first = pn(value("--wall-first")?, "--wall-first")?,
            "--wall-step" => wall_step = pn(value("--wall-step")?, "--wall-step")?,
            "--top-layers" => top_layers = pn(value("--top-layers")?, "--top-layers")?,
            "--bottom-layers" => bottom_layers = pn(value("--bottom-layers")?, "--bottom-layers")?,
            "--layer" => layer = pn(value("--layer")?, "--layer")?,
            "--infill" => infill = pn(value("--infill")?, "--infill")?,
            "--infill-exp" => infill_exp = pn(value("--infill-exp")?, "--infill-exp")?,
            "--z-tensile" => z_tensile = pn(value("--z-tensile")?, "--z-tensile")?,
            "--z-long" => z_long = pn(value("--z-long")?, "--z-long")?,
            "--mass-only" => mass_only = true,
            "--density" => density = pn(value("--density")?, "--density")?,
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
    if cases.is_empty() && !mass_only {
        return Err("give at least one --case".into());
    }
    if !(0.2..=50.0).contains(&voxel) || !(0.01..=1.0).contains(&infill) {
        return Err("--voxel 0.2 to 50 mm; --infill 0.01 to 1".into());
    }
    let library = Library::builtin()?;
    let m = library.find(&material).ok_or(format!("material {material:?} is not in the library"))?;
    let s = m.strengths();
    let (e_gpa, nu) = (s.young_modulus_gpa.ok_or("material has no stiffness")?, s.poisson_ratio.ok_or("material has no Poisson ratio")?);
    let strength = s.tensile_mpa.or(s.yield_mpa).ok_or("material has no strength")?;
    let fix = fix.unwrap_or_else(|| "zmin".into());
    let (fix_axis, fix_max) = match fix.as_str() {
        "xmin" => (0, false),
        "xmax" => (0, true),
        "ymin" => (1, false),
        "ymax" => (1, true),
        "zmin" => (2, false),
        "zmax" => (2, true),
        _ => return Err("--fix: xmin, xmax, ymin, ymax, zmin or zmax".into()),
    };
    let pin = pin.unwrap_or_else(|| "x,0,0,1".into());
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
    let [nx, ny, nz] = grid.n;
    let h = voxel;

    // Walls: material within `wall_mm` of any surface; skins: within top_mm below an upward surface or bottom_mm above a
    // downward one. A cell the wall only partly covers (walls thinner than a voxel) gets the covered fraction of the full
    // stiffness plus the rest as infill (equal-strain mixing along the wall), so thin walls need no thin voxels.
    let wall_mm = wall_first + wall_step * (walls - 1.0);
    let top_mm = top_layers * layer;
    let bottom_mm = bottom_layers * layer;
    let is_solid = |i: i64, j: i64, k: i64| -> bool {
        i >= 0 && j >= 0 && k >= 0 && (i as usize) < nx && (j as usize) < ny && (k as usize) < nz && solid[grid.cell_index(i as usize, j as usize, k as usize)]
    };
    let reach = ((wall_mm + h) / h).ceil() as i64 + 1;
    let up = (top_mm / h).ceil().max(1.0) as i64 + 1;
    let down = (bottom_mm / h).ceil().max(1.0) as i64 + 1;
    let e_core = infill.powf(infill_exp);
    let mut kind = vec![0u8; grid.cells()]; // 0 void, 1 shell (wall or skin covers part of the cell), 2 infill core
    let mut cover = vec![0.0f64; grid.cells()];
    let (mut n_shell, mut n_core) = (0usize, 0usize);
    let mut mass_cm3 = 0.0; // material volume (cm3 of solid plastic)
    let mut total_cm3 = 0.0;
    let cell_cm3 = h.powi(3) / 1000.0;
    for i in 0..nx as i64 {
        for j in 0..ny as i64 {
            for k in 0..nz as i64 {
                if !is_solid(i, j, k) {
                    continue;
                }
                let c = grid.cell_index(i as usize, j as usize, k as usize);
                let mut dmin = f64::INFINITY;
                for di in -reach..=reach {
                    for dj in -reach..=reach {
                        for dk in -reach..=reach {
                            let d = ((di * di + dj * dj + dk * dk) as f64).sqrt() * h;
                            if d < dmin && !is_solid(i + di, j + dj, k + dk) {
                                dmin = d;
                            }
                        }
                    }
                }
                let mut cv = if dmin.is_finite() { ((wall_mm - (dmin - h).max(0.0)) / h).clamp(0.0, 1.0) } else { 0.0 };
                if let Some(s_up) = (1..=up).find(|d| !is_solid(i, j, k + d)) {
                    cv = cv.max(((top_mm - (s_up - 1) as f64 * h) / h).clamp(0.0, 1.0));
                }
                if let Some(s_dn) = (1..=down).find(|d| !is_solid(i, j, k - d)) {
                    cv = cv.max(((bottom_mm - (s_dn - 1) as f64 * h) / h).clamp(0.0, 1.0));
                }
                cover[c] = cv;
                if cv > 0.0 {
                    kind[c] = 1;
                    n_shell += 1;
                } else {
                    kind[c] = 2;
                    n_core += 1;
                }
                total_cm3 += cell_cm3;
                mass_cm3 += cell_cm3 * (cv + (1.0 - cv) * infill);
            }
        }
    }
    let mass = density * mass_cm3;
    eprintln!(
        "{stl}: {nx}x{ny}x{nz} voxels of {h} mm; walls {wall_mm:.2} mm, skins {top_mm:.1}/{bottom_mm:.1} mm: {n_shell} shell + {n_core} core voxels; {total_cm3:.0} cm3 total, predicted mass {mass:.1} g (all solid {:.0} g)",
        density * total_cm3
    );
    if mass_only {
        println!("walls {walls} ({wall_mm:.2} mm) infill {:.1}%: predicted mass {mass:.1} g (all solid {:.1} g)", infill * 100.0, density * total_cm3);
        return Ok(());
    }
    let cell_scale: Vec<f64> = cover.iter().map(|cv| cv + (1.0 - cv) * e_core).collect();
    let solid_used: Vec<bool> = kind.iter().map(|k| *k != 0).collect();
    let (v_solid, v_core) = (n_shell as f64 * cell_cm3, n_core as f64 * cell_cm3);

    let fixed = Region::Side { axis: fix_axis, max: fix_max, depth: 0.5 * voxel };
    let load = Region::Ring { point, dir, r_in: radius - 0.5 * voxel, r_out: radius + voxel };
    let setup = Setup { grid: &grid, solid: &solid_used, e_mpa: e_gpa * 1000.0, nu, fixed: &fixed, load: &load, cg_tol: 1e-8, cg_max: 80_000, exclude_layers: 2, cell_scale: Some(&cell_scale) };
    let result = analyse(&setup, &cases, |msg| eprintln!("  {msg} ..."))?;
    let secs = t0.elapsed().as_secs_f64();

    // Zones left out of the peaks: around the bore (how the load goes in) and at the held face.
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for (c, kd) in kind.iter().enumerate() {
        if *kd != 0 {
            let p = grid.cell_centre(c);
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
    }
    let skip = |p: [f64; 3]| -> bool {
        let rad = (0..3).filter(|k| *k != pa).map(|k| (p[k] - point[k]).powi(2)).sum::<f64>().sqrt();
        let held = if fix_max { hi[fix_axis] - p[fix_axis] } else { p[fix_axis] - lo[fix_axis] };
        rad < radius + 3.0 * h || held < 3.0 * h
    };

    println!(
        "{stl}: {} (E {e_gpa} GPa, in-layer {strength} MPa, across layers {z_tensile} MPa short-term / {z_long} MPa long-term), voxels {h} mm, {secs:.0} s",
        m.name
    );
    println!("walls {walls} ({wall_mm:.2} mm), top {top_layers} / bottom {bottom_layers} layers, infill {:.0}% (core stiffness x{e_core:.4}), predicted mass {mass:.0} g, all solid {:.0} g", infill * 100.0, density * (v_solid + v_core));
    let mut cases_json = Vec::new();
    for c in &result.cases {
        let (mut wall_vm, mut wall_at) = (0.0f64, [0.0; 3]);
        let (mut wall_sz, mut wall_sz_at) = (0.0f64, [0.0; 3]);
        let (mut wall_ix, mut wall_ix_at) = (0.0f64, [0.0; 3]);
        let (mut core_vm, mut core_at) = (0.0f64, [0.0; 3]);
        let mut core_ix = 0.0f64;
        for (cell, st) in c.stress.iter().enumerate() {
            if kind[cell] == 0 || st[0].is_nan() {
                continue;
            }
            let p = grid.cell_centre(cell);
            if skip(p) {
                continue;
            }
            let f = if kind[cell] == 1 { 1.0 / cell_scale[cell] } else { 1.0 }; // wall material's own stress = cell average / stiffness factor
            let vm = c.vm[cell] * f;
            let ix = (st[2].max(0.0).powi(2) + st[4].powi(2) + st[5].powi(2)).sqrt() * f; // across-layer demand: tension along z plus the two shears on the layers
            if kind[cell] == 1 {
                if vm > wall_vm {
                    wall_vm = vm;
                    wall_at = p;
                }
                if st[2] * f > wall_sz {
                    wall_sz = st[2] * f;
                    wall_sz_at = p;
                }
                if ix > wall_ix {
                    wall_ix = ix;
                    wall_ix_at = p;
                }
            } else {
                if vm > core_vm {
                    core_vm = vm;
                    core_at = p;
                }
                core_ix = core_ix.max(ix);
            }
        }
        let strut = |v: f64| v / infill;
        let sf = |allow: f64, v: f64| if v > 0.0 { allow / v } else { f64::INFINITY };
        // The weakest link: walls in the layer plane (vm), walls across layers (ix), infill struts (as walls across layers, the cautious reading).
        let sf_wall_vm = sf(strength, wall_vm);
        let sf_wall_ix_short = sf(z_tensile, wall_ix);
        let sf_wall_ix_long = sf(z_long, wall_ix);
        let sf_core_short = sf(z_tensile, strut(core_ix));
        let sf_core_long = sf(z_long, strut(core_ix));
        println!("{:<8} walls: vm {:.3} MPa (SF {:.1}) at ({:.0}, {:.0}, {:.0}); across layers {:.3} MPa (SF {:.1} short, {:.1} long) at ({:.0}, {:.0}, {:.0}); tension along z {:.3} MPa; infill strut est. {:.2} MPa (SF {:.1} short, {:.1} long); hinge moves {:.3} mm; CG {}",
            c.name, wall_vm, sf_wall_vm, wall_at[0], wall_at[1], wall_at[2], wall_ix, sf_wall_ix_short, sf_wall_ix_long, wall_ix_at[0], wall_ix_at[1], wall_ix_at[2], wall_sz, strut(core_ix), sf_core_short, sf_core_long, c.max_disp, c.cg_iterations);
        cases_json.push(json!({
            "name": c.name, "wall_von_mises_mpa": wall_vm, "wall_von_mises_at_mm": wall_at, "sf_wall_in_layer": sf_wall_vm,
            "wall_across_layer_demand_mpa": wall_ix, "wall_across_layer_at_mm": wall_ix_at, "sf_wall_across_short": sf_wall_ix_short, "sf_wall_across_long": sf_wall_ix_long,
            "wall_tension_z_mpa": wall_sz, "wall_tension_z_at_mm": wall_sz_at,
            "core_average_von_mises_mpa": core_vm, "core_average_at_mm": core_at, "core_strut_estimate_mpa": strut(core_ix),
            "sf_core_across_short": sf_core_short, "sf_core_across_long": sf_core_long,
            "max_displacement_mm": c.max_disp, "max_displacement_at_mm": c.max_disp_at, "cg_iterations": c.cg_iterations, "converged": c.converged,
        }));
        if let Some(dir) = &out {
            std::fs::create_dir_all(dir).map_err(|e| format!("{dir}: {e}"))?;
            let peak = wall_vm.max(core_vm).max(1e-9);
            for (ax, name) in [(0, "along-x"), (1, "along-y"), (2, "along-z")] {
                let (w, hh, rgb) = projection(&grid, &c.vm, ax, peak, 4);
                let path = format!("{dir}/{}-{name}.png", c.name);
                std::fs::write(&path, png(w, hh, &rgb)).map_err(|e| format!("{path}: {e}"))?;
            }
        }
    }
    for w in &result.warnings {
        println!("warning: {w}");
    }
    if let Some(dir) = &out {
        let report = json!({
            "stl": stl, "voxel_mm": h, "grid": grid.n, "material": m.name, "e_gpa": e_gpa, "nu": nu, "in_layer_strength_mpa": strength,
            "across_layer_short_mpa": z_tensile, "across_layer_long_mpa": z_long, "walls": walls, "wall_mm": wall_mm, "top_layers": top_layers,
            "bottom_layers": bottom_layers, "layer_mm": layer, "infill": infill, "infill_exp": infill_exp, "core_stiffness_factor": e_core,
            "voxels": {"shell": n_shell, "core": n_core}, "predicted_mass_g": mass, "all_solid_mass_g": density * (v_solid + v_core),
            "fixed": fix, "pin": pin, "cases": cases_json, "warnings": result.warnings, "seconds": secs,
        });
        let path = format!("{dir}/report.json");
        std::fs::write(&path, serde_json::to_string_pretty(&report).unwrap_or_default()).map_err(|e| format!("{path}: {e}"))?;
        println!("written: {dir}/report.json and {} heat maps", 3 * result.cases.len());
    }
    Ok(())
}
