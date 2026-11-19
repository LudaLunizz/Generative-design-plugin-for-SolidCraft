//! The plate optimiser's extra refine pass, outside SolidCraft: the same optimisation as the
//! plugin's `optimise-plate`, but at up to 800 cells across, solved on the GPU, with the inside
//! corners rounded. It writes the build steps (sketch, extrude, material) as a script that
//! `solidcraft-cli run SCRIPT` or the MCP `batch {path}` runs in the open design.
//!
//! plate_refine [key=value ...] --out SCRIPT.json
//!   width=200 height=100 thickness=10 material=aluminium-6061-t6 load=cantilever share=0.4
//!   cells=480 iterations=150 min_member=3 (mm; the filter radius, so it holds at any grid)
//!   round=4 (fillet radius where members meet, mm) force=500
//!   name="Refined plate"
//!   gpu=1 (default; needs a build with --features gpu, else the CPU multigrid is used)
//!   gpu_max=80 (largest GPU busy share, percent)
//!
//! The reply (stdout) is the same JSON report the plugin gives, plus the script's path.
#![allow(clippy::unwrap_used)] // command-line tool: a failure to write the script is fatal anyway

use generative_design_engine::engine::{Failure, Host, Level, PlateLimits, PlateParams, optimise_plate};
use generative_design_engine::library::Library;
use generative_design_engine::topo2d::Solver;
use serde_json::{Value, json};

/// Records the edit commands instead of running them; answers the reads as an empty design.
struct Recorder {
    steps: Vec<Value>,
    last_progress: f32,
}

impl Host for Recorder {
    type Error = String;
    fn run(&mut self, command: &str, args: &Value) -> Result<Value, Failure<String>> {
        match command {
            "document.inspect" => Ok(json!({ "bodies": [] })),
            "material.list" => Ok(json!({ "materials": [] })),
            "material.create" | "material.assign" | "sketch.create" | "sketch.finish" | "sketch.line" | "sketch.spline.fit_point" => {
                self.steps.push(json!({ "command": command, "params": args }));
                Ok(json!({}))
            }
            "solid.extrude" => {
                self.steps.push(json!({ "command": command, "params": args }));
                Ok(json!({ "created_bodies": [args.get("body_name").cloned().unwrap_or(Value::Null)] }))
            }
            other => Err(Failure::Host(format!("plate_refine does not run {other}"))),
        }
    }
    fn log(&mut self, _: Level, message: &str) {
        eprintln!("{message}");
    }
    fn progress(&mut self, fraction: f32, message: &str) {
        if fraction - self.last_progress >= 0.1 || fraction >= 1.0 {
            eprintln!("{:3.0} %  {message}", fraction * 100.0);
            self.last_progress = fraction;
        }
    }
    fn cancelled(&mut self) -> bool {
        false
    }
}

fn main() {
    let mut kv = std::collections::HashMap::new();
    let mut out = String::from("plate-refine.json");
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        if a == "--out" {
            out = args.next().unwrap_or(out);
        } else if let Some((k, v)) = a.split_once('=') {
            kv.insert(k.to_string(), v.to_string());
        } else {
            eprintln!("unknown argument {a:?} (use key=value or --out FILE)");
            std::process::exit(2);
        }
    }
    let num = |k: &str, d: f64| -> f64 {
        match kv.get(k) {
            None => d,
            Some(v) => v.parse().unwrap_or_else(|_| {
                eprintln!("{k}={v}: not a number");
                std::process::exit(2)
            }),
        }
    };
    let text = |k: &str, d: &str| kv.get(k).cloned().unwrap_or_else(|| d.to_string());
    let gpu = num("gpu", 1.0) != 0.0;
    if gpu && !cfg!(feature = "gpu") {
        eprintln!("gpu=1 needs a build with --features gpu; solving on the CPU");
    }
    #[cfg(feature = "gpu")]
    generative_design_engine::gpu::set_max_busy(num("gpu_max", 80.0).round().clamp(1.0, 100.0) as u32);

    let p = PlateParams {
        width_mm: num("width", 200.0),
        height_mm: num("height", 100.0),
        thickness_mm: num("thickness", 10.0),
        material: text("material", "aluminium-6061-t6"),
        load_case: text("load", "cantilever"),
        force_n: num("force", 500.0),
        volume_share: num("share", 0.4),
        cells: num("cells", 480.0) as usize,
        iterations: num("iterations", 150.0) as usize,
        // The plugin counts the filter in cells; here it is mm, so finer grids keep members printable.
        filter_radius: (num("min_member", 3.0) / (num("width", 200.0) / num("cells", 480.0))).max(1.0),
        round_mm: num("round", 4.0),
        name: text("name", "Refined plate"),
    };
    let limits = PlateLimits {
        max_cells: 400_000,
        max_work: 400_000 * 300,
        max_across: 800,
        max_iterations: 300,
        max_filter: 24.0,
        solver: if gpu { Solver::Gpu } else { Solver::Multigrid },
    };
    let library = Library::builtin().unwrap_or_else(|e| {
        eprintln!("material library: {e}");
        std::process::exit(1)
    });
    let t0 = std::time::Instant::now();
    let mut host = Recorder { steps: Vec::new(), last_progress: -1.0 };
    match optimise_plate(&mut host, &library, p, limits) {
        Ok(mut report) => {
            std::fs::write(&out, serde_json::to_string_pretty(&json!({ "commands": host.steps })).unwrap()).unwrap();
            if let Some(o) = report.as_object_mut() {
                o.insert("script".into(), json!(out));
                o.insert("seconds".into(), json!((t0.elapsed().as_secs_f64() * 10.0).round() / 10.0));
                o.insert("solver".into(), json!(if gpu && cfg!(feature = "gpu") { "GPU multigrid" } else { "CPU multigrid" }));
            }
            println!("{}", serde_json::to_string_pretty(&report).unwrap());
        }
        Err(Failure::Input { field, reason }) => {
            eprintln!("{field}: {reason}");
            std::process::exit(2);
        }
        Err(Failure::Host(e)) | Err(Failure::Internal(e)) => {
            eprintln!("failed: {e}");
            std::process::exit(1);
        }
        Err(Failure::Cancelled) => std::process::exit(1),
    }
}
