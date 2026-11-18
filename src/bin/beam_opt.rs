#![allow(clippy::unwrap_used, clippy::indexing_slicing)]
//! Lightest plate-truss for a long cantilever beam, outside SolidCraft (no plugin limits).
//!
//! beam_opt [key=value ...] --out DIR
//!   length=5000 depth=400 elem=10 load=98.1 pad=100 defl=25 sf=2 lsf=3 rmin=3 its=80
//!   gpu=1 (solve on the GPU; needs a build with --features gpu, else the CPU is used)
//!   jobs=10 (threads at most; default 10)
//!   e=68900 nu=0.33 rho=2.70 yield=276 thick=8,10,12,15,20,25,30 weight=1 tol=
use generative_design_engine::beam::{self, BeamSpec};
use std::io::Write;

/// A cleaned design (solid cells, tip attached) or why the optimiser failed, per volume fraction.
type Topology = (f64, Result<(Vec<bool>, bool), String>);

fn main() {
    let mut kv = std::collections::HashMap::new();
    let mut out = String::from("beam-out");
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        if a == "--out" {
            out = args.next().unwrap_or(out);
        } else if let Some((k, v)) = a.split_once('=') {
            kv.insert(k.to_string(), v.to_string());
        }
    }
    let num = |k: &str, d: f64| kv.get(k).and_then(|v| v.parse::<f64>().ok()).unwrap_or(d);
    let s = BeamSpec {
        length: num("length", 5000.0),
        depth: num("depth", 400.0),
        tip_depth: num("tip", num("depth", 400.0)),
        min_member_mm: num("minmember", 0.0),
        elem: num("elem", 10.0),
        e_mpa: num("e", 68_900.0),
        nu: num("nu", 0.33),
        density_g_cm3: num("rho", 2.70),
        yield_mpa: num("yield", 276.0),
        load_n: num("load", 98.1),
        pad_mm: num("pad", 100.0),
        safety: num("sf", 2.0),
        buckling_safety: num("lsf", 3.0),
        max_defl_mm: num("defl", 25.0),
        rmin_cells: num("rmin", 3.0),
        iterations: num("its", 80.0) as usize,
        self_weight: num("weight", 1.0) != 0.0,
        gpu: num("gpu", 0.0) != 0.0,
        member_len_mm: num("member", 1.5 * 0.5 * (num("depth", 400.0) + num("tip", num("depth", 400.0)))),
    };
    if s.gpu && !cfg!(feature = "gpu") {
        eprintln!("gpu=1 needs a build with --features gpu; solving on the CPU");
    }
    if let Err(e) = s.check() {
        eprintln!("bad parameters: {e}");
        std::process::exit(2);
    }
    let thick: Vec<f64> = kv.get("thick").map(|v| v.split(',').filter_map(|x| x.parse().ok()).collect()).unwrap_or_else(|| vec![6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0]);
    std::fs::create_dir_all(&out).unwrap();
    println!("beam {} mm long, {} mm deep at the root, {} at the tip, min member {} mm, cells {} x {} of {} mm, load {} N, deflection <= {} mm, SF {} (stability {})", s.length, s.depth, s.tip_depth, s.min_member_mm, s.nx(), s.ny(), s.elem, s.load_n, s.max_defl_mm, s.safety, s.buckling_safety);
    let vfs: Vec<f64> = kv.get("vf").map(|v| v.split(',').filter_map(|x| x.parse().ok()).collect()).unwrap_or_else(|| vec![0.06, 0.08, 0.10, 0.12, 0.15, 0.18, 0.22, 0.27, 0.33, 0.40]);
    // One optimisation per volume fraction, all at once.
    let jobs = (num("jobs", 10.0) as usize).clamp(1, 10);
    let mut tops: Vec<Topology> = Vec::new();
    for chunk in vfs.chunks(jobs) {
        let done: Vec<Topology> = std::thread::scope(|sc| {
            let handles: Vec<_> = chunk.iter().map(|&vf| { let s = &s; sc.spawn(move || (vf, beam::topology(s, vf))) }).collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });
        tops.extend(done);
    }
    let mut best: Option<beam::Design> = None;
    println!("\n  vf in  t mm   vol%    mass kg  defl mm (no weight)  p99 MPa  member mm  strut margin  Pcr/P  verdict");
    for (vf, r) in &tops {
        let (solid, top_ok) = match r {
            Ok(v) => v,
            Err(e) => {
                println!("  {vf:.2}  error: {e}");
                continue;
            }
        };
        let evals: Vec<(f64, Result<beam::Design, String>)> = std::thread::scope(|sc| {
            let hs: Vec<_> = thick.iter().map(|&t| { let (s, solid) = (&s, solid); sc.spawn(move || (t, beam::evaluate(s, solid, *top_ok, t))) }).collect();
            hs.into_iter().map(|h| h.join().unwrap()).collect()
        });
        for (t, e) in evals {
            let Ok(d) = e else { continue };
            let ratio = d.ltb_pcr_n / s.load_n;
            let why = if !d.connected { "tip not connected" } else if d.tip_defl_mm > s.max_defl_mm { "too soft" } else if d.vm_p99 * s.safety > s.yield_mpa { "overstressed" } else if d.member_margin < s.safety { "struts buckle" } else if ratio < s.buckling_safety { "tips sideways" } else { "OK" };
            println!("  {vf:.2}  {t:5.1}  {:5.1}  {:8.2}  {:6.1} ({:5.1})     {:7.1}  {:8.1}  {:10.1}  {:6.1}  {why}", d.volfrac * 100.0, d.mass_kg, d.tip_defl_mm, d.tip_defl_no_weight_mm, d.vm_p99, d.member_width_mm, d.member_margin, ratio);
            if why == "OK" && best.as_ref().is_none_or(|b| d.mass_kg < b.mass_kg) {
                best = Some(d);
            }
        }
    }
    let Some(d) = best else {
        println!("\nno thickness gives a design that is also laterally stable");
        std::process::exit(1);
    };
    let loops = beam::outline(&d, s.elem, num("tol", 0.7 * s.elem));
    let cells_area = d.solid.iter().filter(|b| **b).count() as f64 * s.elem * s.elem;
    let poly_area: f64 = loops.iter().map(|l| beam::area(l)).sum();
    println!("\nbest: thickness {} mm, {:.2} kg, {} loops, {} vertices, outline area {:+.2} % against the cells", d.thickness, d.mass_kg, loops.len(), loops.iter().map(Vec::len).sum::<usize>(), (poly_area / cells_area - 1.0) * 100.0);
    let json = serde_json::json!({
        "length": s.length, "depth": s.depth, "tip_depth": s.tip_depth, "min_member_mm": s.min_member_mm, "thickness": d.thickness, "elem": s.elem,
        "mass_kg": d.mass_kg, "volume_fraction": d.volfrac, "tip_deflection_mm": d.tip_defl_mm,
        "tip_deflection_no_weight_mm": d.tip_defl_no_weight_mm, "vm_p99_mpa": d.vm_p99, "vm_max_mpa": d.vm_max,
        "self_weight_n": d.weight_n, "member_width_mm": d.member_width_mm, "member_margin": d.member_margin, "ltb_pcr_n": d.ltb_pcr_n, "outline_area_mm2": poly_area,
        "loops": loops.iter().map(|l| l.iter().map(|(x, y)| [x, y]).collect::<Vec<_>>()).collect::<Vec<_>>(),
    });
    std::fs::File::create(format!("{out}/design.json")).unwrap().write_all(serde_json::to_string(&json).unwrap().as_bytes()).unwrap();
    std::fs::write(format!("{out}/design.svg"), beam::svg(&loops, s.length, s.depth)).unwrap();
    println!("wrote {out}/design.json and design.svg");
}
