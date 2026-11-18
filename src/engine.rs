//! The plugin's entries, written against [`Host`] so they run the same in SolidCraft (through the
//! WIT bindings in `component.rs`) and in native tests (through a fake host).
//!
//! Inputs are hostile: every parameter object refuses unknown fields, numbers are range-checked,
//! names are bounded, and nothing here panics.

use crate::library::{Library, Material};
use serde::Deserialize;
use serde_json::{Value, json};
use std::collections::BTreeSet;

/// The largest number of bodies one call takes (SolidCraft's own commands cap lists at 10 000).
pub const MAX_BODIES: usize = 10_000;



#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Level {
    Debug,
    Info,
    Warn,
    Error,
}

/// Why an entry failed. `Host` carries the app's own error unchanged.
#[derive(Debug, PartialEq)]
pub enum Failure<E> {
    Input { field: String, reason: String },
    Host(E),
    Cancelled,
    Internal(String),
}

impl<E> Failure<E> {
    fn input(field: &str, reason: impl Into<String>) -> Self {
        Self::Input { field: field.into(), reason: reason.into() }
    }
}

pub trait Host {
    type Error;
    /// Runs one app command. The result is the command's JSON reply.
    fn run(&mut self, command: &str, args: &Value) -> Result<Value, Failure<Self::Error>>;
    fn log(&mut self, level: Level, message: &str);
    fn progress(&mut self, fraction: f32, message: &str);
    fn cancelled(&mut self) -> bool;
}

/// Runs entry `entry` with JSON `params`, returning the JSON result.
pub fn execute<H: Host>(host: &mut H, entry: &str, params: &str) -> Result<String, Failure<H::Error>> {
    let library = Library::builtin().map_err(Failure::Internal)?;
    let result = match entry {
        "materials" => materials(&library, parse(params)?),
        "apply-material" => apply_material(host, &library, parse(params)?)?,
        "material-report" => material_report(host, &library, parse(params)?)?,
        "solver-benchmark" => solver_benchmark(host, parse(params)?)?,
        "stress-check" => stress_check(host, &library, parse(params)?)?,
        other => return Err(Failure::input("entry", format!("unknown entry {:?}", truncate(other, 64)))),
    };
    serde_json::to_string(&result).map_err(|e| Failure::Internal(format!("could not encode the result: {e}")))
}

fn parse<T: for<'de> Deserialize<'de>, E>(params: &str) -> Result<T, Failure<E>> {
    let text = if params.trim().is_empty() { "{}" } else { params };
    let value: Value = serde_json::from_str(text).map_err(|e| Failure::input("params", e.to_string()))?;
    if !value.is_object() {
        return Err(Failure::input("params", "expected a JSON object"));
    }
    serde_json::from_value(value).map_err(|e| Failure::input("params", e.to_string()))
}

fn truncate(s: &str, max: usize) -> String {
    s.chars().take(max).collect()
}

fn check_safety_factor<E>(sf: f64) -> Result<f64, Failure<E>> {
    if sf.is_finite() && (1.0..=10.0).contains(&sf) {
        Ok(sf)
    } else {
        Err(Failure::input("safety_factor", "expected a number from 1 to 10"))
    }
}

fn check_bodies<E>(bodies: &[String]) -> Result<(), Failure<E>> {
    if bodies.len() > MAX_BODIES {
        return Err(Failure::input("bodies", format!("at most {MAX_BODIES} bodies")));
    }
    for b in bodies {
        if b.trim().is_empty() || b.chars().count() > 256 || b.chars().any(char::is_control) {
            return Err(Failure::input("bodies", "each body name is 1-256 printable characters"));
        }
    }
    let unique: BTreeSet<_> = bodies.iter().collect();
    if unique.len() != bodies.len() {
        return Err(Failure::input("bodies", "a body is named twice"));
    }
    Ok(())
}

/// `body` (one name, from the GUI's text field) and `bodies` (a list, from scripts) together.
fn merge_body<E>(mut bodies: Vec<String>, body: &str) -> Result<Vec<String>, Failure<E>> {
    let body = body.trim();
    if !body.is_empty() {
        if !bodies.is_empty() {
            return Err(Failure::input("body", "give body or bodies, not both"));
        }
        bodies.push(body.to_string());
    }
    Ok(bodies)
}

fn round(v: f64, digits: i32) -> f64 {
    let f = 10f64.powi(digits);
    (v * f).round() / f
}

// ---------------------------------------------------------------------------------------------
// materials: list the library (read, no app commands)

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields, default)]
struct MaterialsParams {
    family: FamilyFilter,
    safety_factor: Option<f64>,
}

#[derive(Debug, Default, Deserialize, Clone, Copy, PartialEq)]
#[serde(rename_all = "snake_case")]
enum FamilyFilter {
    #[default]
    Any,
    Metal,
    Polymer,
}

fn materials(library: &Library, p: MaterialsParams) -> Value {
    let sf = p.safety_factor.filter(|v| v.is_finite() && (1.0..=10.0).contains(v)).unwrap_or(2.0);
    let rows: Vec<Value> = library
        .materials
        .iter()
        .filter(|m| match p.family {
            FamilyFilter::Any => true,
            FamilyFilter::Metal => m.family.as_str() == "metal",
            FamilyFilter::Polymer => m.family.as_str() == "polymer",
        })
        .map(|m| material_row(library, m, sf))
        .collect();
    json!({ "safety_factor": sf, "count": rows.len(), "materials": rows })
}

fn material_row(library: &Library, m: &Material, sf: f64) -> Value {
    let s = m.strengths();
    let allowable = s.allowable_mpa(sf);
    let sources: BTreeSet<&str> = m.properties.iter().map(|p| p.source_id.as_str()).collect();
    let illustrative = sources
        .iter()
        .any(|id| library.sources.get(*id).is_some_and(|src| src.kind == crate::library::SourceKind::Illustrative));
    json!({
        "id": m.id,
        "name": m.name,
        "family": m.family,
        "designation": m.designation,
        "density_g_cm3": s.density_kg_m3.map(|d| round(d / 1000.0, 4)),
        "young_modulus_gpa": s.young_modulus_gpa.map(|v| round(v, 2)),
        "poisson_ratio": s.poisson_ratio,
        "yield_mpa": s.yield_mpa.map(|v| round(v, 1)),
        "tensile_mpa": s.tensile_mpa.map(|v| round(v, 1)),
        "shear_mpa": s.shear_mpa.map(|v| round(v, 1)),
        "allowable_mpa": allowable.map(|(a, _)| round(a, 1)),
        "allowable_basis": allowable.map(|(_, b)| b),
        "illustrative_values": illustrative,
        "sources": sources,
        "notes": m.notes,
    })
}

// ---------------------------------------------------------------------------------------------
// apply-material: create the library material in the design (once) and assign it (edit)

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ApplyParams {
    material: String,
    #[serde(default)]
    bodies: Vec<String>,
    /// One body by name (the GUI dialog draws text, not lists). Empty: see `bodies`.
    #[serde(default)]
    body: String,
}

fn apply_material<H: Host>(host: &mut H, library: &Library, mut p: ApplyParams) -> Result<Value, Failure<H::Error>> {
    p.bodies = merge_body(std::mem::take(&mut p.bodies), &p.body)?;
    let material = library.find(p.material.trim()).ok_or_else(|| {
        let ids: Vec<&str> = library.materials.iter().map(|m| m.id.as_str()).collect();
        Failure::input("material", format!("not in the library; use one of: {}", ids.join(", ")))
    })?;
    check_bodies(&p.bodies)?;

    let existing_bodies = body_names(host)?;
    if existing_bodies.is_empty() {
        return Err(Failure::input("bodies", "the design has no bodies"));
    }
    let bodies = if p.bodies.is_empty() {
        existing_bodies
    } else {
        let known: BTreeSet<&String> = existing_bodies.iter().collect();
        let missing: Vec<&str> = p.bodies.iter().filter(|b| !known.contains(b)).map(String::as_str).take(10).collect();
        if !missing.is_empty() {
            return Err(Failure::input("bodies", format!("no body named {}", missing.join(", "))));
        }
        p.bodies
    };

    if host.cancelled() {
        return Err(Failure::Cancelled);
    }
    let created = ensure_material(host, library, material)?;
    if host.cancelled() {
        return Err(Failure::Cancelled);
    }
    host.progress(0.8, "Assigning the material");
    let assigned = host.run("material.assign", &json!({ "bodies": bodies, "material": material.name }))?;
    let assigned_bodies = assigned.get("bodies").cloned().unwrap_or(Value::Array(Vec::new()));
    host.log(Level::Info, &format!("{} assigned to {} bodies", material.name, assigned_bodies.as_array().map_or(0, Vec::len)));
    Ok(json!({
        "material": material.name,
        "material_id": material.id,
        "created": created,
        "bodies": assigned_bodies,
    }))
}

fn body_names<H: Host>(host: &mut H) -> Result<Vec<String>, Failure<H::Error>> {
    let doc = host.run("document.inspect", &json!({ "bodies": "names", "timeline": false }))?;
    let names = doc
        .get("bodies")
        .and_then(Value::as_array)
        .ok_or_else(|| Failure::Internal("document.inspect gave no body list".into()))?
        .iter()
        .filter_map(|b| b.as_str().map(str::to_owned).or_else(|| b.get("name").and_then(Value::as_str).map(str::to_owned)))
        .collect();
    Ok(names)
}

/// Makes sure the design has `material` with exactly the library's card. Returns whether it was
/// created now. A design material of the same name with other values is never overwritten.
fn ensure_material<H: Host>(host: &mut H, library: &Library, material: &Material) -> Result<bool, Failure<H::Error>> {
    let list = host.run("material.list", &json!({}))?;
    let existing = list
        .get("materials")
        .and_then(Value::as_array)
        .and_then(|all| all.iter().find(|m| m.get("name").and_then(Value::as_str).is_some_and(|n| n.trim().eq_ignore_ascii_case(&material.name))));
    if let Some(existing) = existing {
        if same_material(existing, library, material) {
            host.log(Level::Info, &format!("{} is already in the design", material.name));
            return Ok(false);
        }
        return Err(Failure::input(
            "material",
            format!(
                "the design already has a material named {:?} with other values; rename or delete it first, then apply again",
                material.name
            ),
        ));
    }
    let params = material.create_params(library).ok_or_else(|| Failure::Internal(format!("{} has no density", material.id)))?;
    host.progress(0.4, "Creating the material");
    host.run("material.create", &params)?;
    host.log(Level::Info, &format!("Created material {}", material.name));
    Ok(true)
}

/// Whether a `material.list` row carries the same density and card as the library entry.
fn same_material(row: &Value, library: &Library, material: &Material) -> bool {
    let Some(density) = material.density_kg_m3() else { return false };
    let row_density = row.get("density_g_cm3").and_then(Value::as_f64);
    if row_density.is_none_or(|d| (d * 1000.0 - density).abs() > 1e-6 * density) {
        return false;
    }
    let Some(card) = row.get("physical_card") else { return false };
    let ours = material.card(library);
    let props = |c: &Value| -> Option<Vec<crate::library::Property>> { serde_json::from_value(c.get("properties")?.clone()).ok() };
    card.get("designation") == ours.get("designation")
        && card.get("family") == ours.get("family")
        && props(card).is_some()
        && props(card) == props(&ours)
}

// ---------------------------------------------------------------------------------------------
// material-report: mass, material and allowable stress per body (read)

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ReportParams {
    #[serde(default)]
    bodies: Vec<String>,
    #[serde(default)]
    body: String,
    #[serde(default = "default_sf")]
    safety_factor: f64,
}

fn default_sf() -> f64 {
    2.0
}

fn material_report<H: Host>(host: &mut H, library: &Library, mut p: ReportParams) -> Result<Value, Failure<H::Error>> {
    p.bodies = merge_body(std::mem::take(&mut p.bodies), &p.body)?;
    let sf = check_safety_factor(p.safety_factor)?;
    check_bodies(&p.bodies)?;
    let args = if p.bodies.is_empty() { json!({}) } else { json!({ "bodies": p.bodies }) };
    let measured = host.run("inspect.measure", &args)?;
    let rows = measured
        .get("bodies")
        .and_then(Value::as_array)
        .ok_or_else(|| Failure::Internal("inspect.measure gave no body list".into()))?;

    let mut total_mass = 0.0;
    let mut report = Vec::with_capacity(rows.len());
    for b in rows {
        let name = b.get("name").and_then(Value::as_str).unwrap_or("?");
        if let Some(error) = b.get("error").and_then(Value::as_str) {
            report.push(json!({ "body": name, "error": error }));
            continue;
        }
        let material = b.get("material").and_then(Value::as_str);
        let mass = b.get("mass_g").and_then(Value::as_f64);
        total_mass += mass.unwrap_or(0.0);
        // Strengths come from the library, matched by the material's name. SolidCraft 0.42 does
        // not class material.list or material.inspect as read-only, so a read entry can't read the
        // design's card. The density SolidCraft measured with guards against a same-named
        // material whose values were changed in the design.
        let entry = material.and_then(|m| library.materials.iter().find(|l| l.name.eq_ignore_ascii_case(m.trim())));
        let density = b.get("density_g_cm3").and_then(Value::as_f64);
        let (s, note) = match entry {
            None => (Default::default(), Some("not a library material, so no strength values; apply one with apply-material")),
            Some(l) => match (l.density_kg_m3(), density) {
                (Some(lib), Some(d)) if (d * 1000.0 - lib).abs() <= 1e-6 * lib => (l.strengths(), None),
                _ => (Default::default(), Some("this material's density differs from the library entry of that name; strengths left out")),
            },
        };
        let allowable = s.allowable_mpa(sf);
        report.push(json!({
            "body": name,
            "material": material,
            "library_id": entry.map(|l| l.id.as_str()),
            "closed": b.get("closed"),
            "volume_mm3": b.get("volume_mm3"),
            "mass_g": mass,
            "young_modulus_gpa": s.young_modulus_gpa.map(|v| round(v, 2)),
            "yield_mpa": s.yield_mpa.map(|v| round(v, 1)),
            "tensile_mpa": s.tensile_mpa.map(|v| round(v, 1)),
            "allowable_mpa": allowable.map(|(a, _)| round(a, 1)),
            "allowable_basis": allowable.map(|(_, b)| b),
            "note": note,
        }));
    }
    Ok(json!({
        "safety_factor": sf,
        "body_count": report.len(),
        "total_mass_g": round(total_mass, 3),
        "bodies": report,
        "stress_analysis": "not done: this report gives allowable stress only; no solver has run",
    }))
}

// ---------------------------------------------------------------------------------------------
// solver-benchmark: run 2D SIMP on a standard test case and report the work done (read, no app
// commands). Time is measured outside (solidcraft-cli --timings): without the clock permission
// the plugin's clock is a fixed epoch.

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct BenchParams {
    pub nx: usize,
    pub ny: usize,
    pub iterations: usize,
    pub load: String,
    pub volfrac: f64,
    pub rmin: f64,
    pub cg_tol: f64,
    pub cg_max: usize,
    pub picture: bool,
    pub solver: String,
    /// Profiling: run only one part `repeat` times ("apply", "setup", "vcycle") instead of SIMP.
    pub part: String,
    pub repeat: usize,
}

impl Default for BenchParams {
    fn default() -> Self {
        Self { nx: 60, ny: 30, iterations: 50, load: "cantilever".into(), volfrac: 0.5, rmin: 1.5, cg_tol: 1e-6, cg_max: 20_000, picture: true, solver: "multigrid".into(), part: String::new(), repeat: 1 }
    }
}

impl BenchParams {
    pub fn problem<E>(&self) -> Result<crate::topo2d::Problem, Failure<E>> {
        use crate::topo2d::{LoadCase, Problem, Solver};
        let solver = match self.solver.as_str() {
            "multigrid" => Solver::Multigrid,
            "jacobi" => Solver::Jacobi,
            _ => return Err(Failure::input("solver", "expected multigrid or jacobi")),
        };
        let load = match self.load.as_str() {
            "cantilever" => LoadCase::Cantilever,
            "mbb" => LoadCase::Mbb,
            _ => return Err(Failure::input("load", "expected cantilever or mbb")),
        };
        let p = Problem {
            nx: self.nx,
            ny: self.ny,
            volfrac: self.volfrac,
            penal: 3.0,
            rmin: self.rmin,
            nu: 0.3,
            load,
            max_iterations: self.iterations,
            change_tol: 0.0,
            cg_tol: self.cg_tol,
            cg_max: self.cg_max,
            solver,
            beta_final: 0.0,
            beta_every: 30,
            tip_rows: 0,
            coarse_budget: 0,
            pad_rows: 0,
        };
        p.check().map_err(|reason| Failure::input("params", reason))?;
        Ok(p)
    }
}

/// The benchmark's JSON report for a finished run.
pub fn bench_report(p: &crate::topo2d::Problem, out: &crate::topo2d::Outcome, picture: bool) -> Value {
    let cg: Vec<usize> = out.history.iter().map(|i| i.cg_iterations).collect();
    let total_cg: usize = cg.iter().sum();
    let elements = p.nx * p.ny;
    let compliance: Vec<f64> = out.history.iter().map(|i| round(i.compliance, 4)).collect();
    json!({
        "elements": elements,
        "dofs": p.ndof(),
        "design_iterations": out.history.len(),
        "cg_iterations_total": total_cg,
        "cg_iterations_per_design_iteration": cg,
        "all_solves_converged": out.history.iter().all(|i| i.cg_converged),
        "solver": match p.solver { crate::topo2d::Solver::Jacobi => "jacobi", crate::topo2d::Solver::Multigrid => "multigrid", crate::topo2d::Solver::Gpu => "gpu" },
        "multigrid_levels": out.history.first().map(|i| i.levels),
        // One stiffness product touches every element once: 8x8 multiply-adds per element.
        "element_products": (total_cg + out.history.len()) * elements,
        "compliance_first": compliance.first(),
        "compliance_last": compliance.last(),
        "volume_last": out.history.last().map(|i| round(i.volume, 4)),
        "change_last": out.history.last().map(|i| round(i.change, 4)),
        "cancelled": out.cancelled,
        "picture": if picture { json!(crate::topo2d::ascii(&out.density, p.nx, p.ny, 60)) } else { Value::Null },
    })
}

fn solver_benchmark<H: Host>(host: &mut H, b: BenchParams) -> Result<Value, Failure<H::Error>> {
    let p = b.problem()?;
    if !b.part.is_empty() {
        if !(1..=10_000).contains(&b.repeat) {
            return Err(Failure::input("repeat", "expected 1 to 10000"));
        }
        return crate::topo2d::profile_part(&p, &b.part, b.repeat).map_err(|reason| Failure::input("part", reason));
    }
    let total = p.max_iterations;
    let out = crate::topo2d::optimise(&p, |it, i| {
        host.progress((it + 1) as f32 / total as f32, &format!("Iteration {} of {total}, compliance {:.3}", it + 1, i.compliance));
        !host.cancelled()
    });
    if out.cancelled {
        return Err(Failure::Cancelled);
    }
    Ok(bench_report(&p, &out, b.picture))
}

// ---------------------------------------------------------------------------------------------
// stress-check: plane-stress analysis of a flat body under one load (read)

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct StressParams {
    pub body: String,
    pub fixed: String,
    pub load: String,
    pub force_x: f64,
    pub force_y: f64,
    pub force_z: f64,
    pub safety_factor: f64,
    pub material: String,
    pub plane: String,
    pub resolution: usize,
    pub picture: bool,
}

impl Default for StressParams {
    fn default() -> Self {
        Self {
            body: String::new(),
            fixed: String::new(),
            load: String::new(),
            force_x: 0.0,
            force_y: -100.0,
            force_z: 0.0,
            safety_factor: 2.0,
            material: "body".into(),
            plane: "auto".into(),
            resolution: 240,
            picture: true,
        }
    }
}

fn stress_check<H: Host>(host: &mut H, library: &Library, p: StressParams) -> Result<Value, Failure<H::Error>> {
    use crate::stress2d::{Face, FaceGeom, Input, PlaneChoice, check};
    let sf = check_safety_factor(p.safety_factor)?;
    if !(10..=400).contains(&p.resolution) {
        return Err(Failure::input("resolution", "expected 10 to 400 cells"));
    }
    let plane = match p.plane.as_str() {
        "auto" => PlaneChoice::Auto,
        "xy" => PlaneChoice::Xy,
        "xz" => PlaneChoice::Xz,
        "yz" => PlaneChoice::Yz,
        _ => return Err(Failure::input("plane", "expected auto, xy, xz or yz")),
    };
    for (field, value) in [("fixed", &p.fixed), ("load", &p.load)] {
        if value.trim().is_empty() || value.chars().count() > 2048 || value.chars().any(char::is_control) {
            return Err(Failure::input(field, "name a face (e.g. F1:-x), its number, or a side such as -x; several with commas"));
        }
    }
    let force = [p.force_x, p.force_y, p.force_z];
    if force.iter().any(|f| !f.is_finite() || f.abs() > 1e9) {
        return Err(Failure::input("force", "expected numbers up to 1e9 N"));
    }
    if p.body.chars().count() > 256 || p.body.chars().any(char::is_control) {
        return Err(Failure::input("body", "at most 256 printable characters"));
    }

    // The body and its material.
    host.progress(0.05, "Reading the body");
    let args = if p.body.trim().is_empty() { json!({}) } else { json!({ "bodies": [p.body.trim()] }) };
    let measured = host.run("inspect.measure", &args)?;
    let rows: Vec<&Value> = measured.get("bodies").and_then(Value::as_array).map(|r| r.iter().collect()).unwrap_or_default();
    let row = match rows.as_slice() {
        [one] => *one,
        [] => return Err(Failure::input("body", "the design has no bodies")),
        many => {
            let names: Vec<&str> = many.iter().filter_map(|b| b.get("name").and_then(Value::as_str)).take(10).collect();
            return Err(Failure::input("body", format!("the design has several bodies; name one: {}", names.join(", "))));
        }
    };
    let body = row.get("name").and_then(Value::as_str).unwrap_or_default().to_string();
    if let Some(error) = row.get("error").and_then(Value::as_str) {
        return Err(Failure::input("body", format!("{body} cannot be measured: {error}")));
    }
    if row.get("closed").and_then(Value::as_bool) == Some(false) {
        return Err(Failure::input("body", format!("{body} is not closed, so it has no inside to analyse")));
    }
    let body_material = row.get("material").and_then(Value::as_str).unwrap_or("");
    let (material, source) = if p.material == "body" {
        let entry = library.materials.iter().find(|m| m.name.eq_ignore_ascii_case(body_material.trim())).ok_or_else(|| {
            Failure::input(
                "material",
                format!("{body}'s material ({body_material}) has no strength values; apply a library material first, or pick one here"),
            )
        })?;
        let density = row.get("density_g_cm3").and_then(Value::as_f64);
        if let (Some(lib), Some(d)) = (entry.density_kg_m3(), density)
            && (d * 1000.0 - lib).abs() > 1e-6 * lib
        {
            return Err(Failure::input("material", format!("{body}'s material is named {} but has other values than the library's; pick one here", entry.name)));
        }
        (entry, "the body's material")
    } else {
        let entry = library.find(&p.material).ok_or_else(|| Failure::input("material", "not in the library"))?;
        (entry, "chosen for this check")
    };
    let s = material.strengths();
    let (Some(e_gpa), Some(nu)) = (s.young_modulus_gpa, s.poisson_ratio) else {
        return Err(Failure::input("material", format!("{} has no stiffness values", material.name)));
    };
    let (strength, basis) = match (s.yield_mpa, s.tensile_mpa) {
        (Some(y), _) => (y, "yield_strength"),
        (None, Some(t)) => (t, "tensile_strength"),
        _ => return Err(Failure::input("material", format!("{} has no strength values", material.name))),
    };

    // Geometry and faces.
    let longest = row
        .get("bbox")
        .and_then(|b| Some((b.get("min")?.as_array()?.clone(), b.get("max")?.as_array()?.clone())))
        .map(|(lo, hi)| lo.iter().zip(&hi).filter_map(|(l, h)| Some(h.as_f64()? - l.as_f64()?)).fold(0.0, f64::max))
        .unwrap_or(100.0);
    if host.cancelled() {
        return Err(Failure::Cancelled);
    }
    host.progress(0.15, "Reading the mesh and faces");
    let mesh = host.run("model.mesh", &json!({ "body": body, "tolerance": (longest / 400.0).max(0.05) }))?;
    let positions: Vec<crate::stress2d::P3> = mesh
        .get("positions")
        .and_then(Value::as_array)
        .map(|ps| ps.iter().filter_map(|p| Some([p.get(0)?.as_f64()?, p.get(1)?.as_f64()?, p.get(2)?.as_f64()?])).collect())
        .unwrap_or_default();
    let triangles: Vec<[u32; 3]> = mesh
        .get("triangles")
        .and_then(Value::as_array)
        .map(|ts| {
            ts.iter()
                .filter_map(|t| {
                    let i = |k: usize| t.get(k)?.as_u64().and_then(|v| u32::try_from(v).ok());
                    Some([i(0)?, i(1)?, i(2)?])
                })
                .collect()
        })
        .unwrap_or_default();
    let faces_json = host.run("model.faces", &json!({ "body": body }))?;
    let p3 = |v: Option<&Value>| -> Option<crate::stress2d::P3> {
        let a = v?.as_array()?;
        Some([a.first()?.as_f64()?, a.get(1)?.as_f64()?, a.get(2)?.as_f64()?])
    };
    let faces: Vec<Face> = faces_json
        .get("faces")
        .and_then(Value::as_array)
        .map(|fs| {
            fs.iter()
                .filter_map(|f| {
                    let geom = if let Some(n) = p3(f.get("plane_normal")) {
                        FaceGeom::Plane { normal: n }
                    } else if let Some(c) = f.get("cylinder") {
                        FaceGeom::Cylinder { axis: p3(c.get("axis"))?, origin: p3(c.get("origin"))?, radius: c.get("radius")?.as_f64()? }
                    } else {
                        FaceGeom::Other
                    };
                    Some(Face {
                        index: usize::try_from(f.get("index")?.as_u64()?).ok()?,
                        name: f.get("name")?.as_str()?.to_string(),
                        centroid: p3(f.get("centroid"))?,
                        area: f.get("area")?.as_f64()?,
                        geom,
                    })
                })
                .collect()
        })
        .unwrap_or_default();
    if host.cancelled() {
        return Err(Failure::Cancelled);
    }
    host.progress(0.3, "Solving");
    let input = Input {
        positions: &positions,
        triangles: &triangles,
        faces: &faces,
        plane,
        resolution: p.resolution,
        fixed: &p.fixed,
        load: &p.load,
        force,
        e_mpa: e_gpa * 1000.0,
        nu,
    };
    let out = check(&input).map_err(|reason| Failure::input("check", reason))?;
    host.progress(1.0, "Done");

    let found_sf = if out.max_vm > 0.0 { strength / out.max_vm } else { f64::INFINITY };
    let passes = found_sf >= sf;
    let axis = |k: usize| ["x", "y", "z"].get(k).copied().unwrap_or("?");
    let plane_name = format!("{}{}", axis(out.axes[0]), axis(out.axes[1]));
    let pt = |v: [f64; 3]| v.map(|c| round(c, 2));
    let verdict = format!(
        "{}: peak {:.1} MPa against {:.0} MPa {} gives a safety factor of {:.2} (needed {:.2})",
        if passes { "Passes" } else { "Does not pass" },
        out.max_vm,
        strength,
        if basis == "yield_strength" { "yield" } else { "tensile strength" },
        found_sf,
        sf
    );
    host.log(Level::Info, &verdict);
    Ok(json!({
        "verdict": verdict,
        "passes": passes,
        "body": body,
        "material": material.name,
        "material_source": source,
        "strength_mpa": strength,
        "strength_basis": basis,
        "safety_factor_needed": sf,
        "safety_factor_found": if found_sf.is_finite() { json!(round(found_sf, 3)) } else { Value::Null },
        "max_von_mises_mpa": round(out.max_vm, 3),
        "max_von_mises_at_mm": pt(out.max_vm_at),
        "p99_von_mises_mpa": round(out.p99_vm, 3),
        "max_displacement_mm": round(out.max_disp, 5),
        "max_displacement_at_mm": pt(out.max_disp_at),
        "plane": plane_name,
        "thickness_mm": round(out.thickness, 4),
        "grid": { "cells_x": out.nx, "cells_y": out.ny, "cell_mm": round(out.cell, 4), "material_cells": out.solid_cells },
        "supports": { "faces": p.fixed.trim(), "grid_nodes": out.fixed_nodes },
        "load": { "faces": p.load.trim(), "grid_nodes": out.load_nodes, "force_n": force },
        "solver": { "cg_iterations": out.cg_iterations, "converged": out.converged },
        "warnings": out.warnings,
        "picture": if p.picture { json!(out.picture) } else { Value::Null },
        "method": "2D plane stress on a grid through the part's plane: linear elastic, small deflections, rigid supports, the load spread evenly over the load face. The peak is read at the surface of the part (a plane fitted to the stress of the cells next to an edge, evaluated on the edge), so it does not depend on where a round edge falls on the grid. A design check, not a certified analysis.",
    }))
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
    use super::*;

    /// A fake SolidCraft holding bodies and design materials, recording every command.
    #[derive(Default)]
    struct Fake {
        bodies: Vec<(String, f64, String)>, // name, volume mm3, material
        materials: Vec<Value>,              // design materials as material.list rows
        calls: Vec<String>,
        cancel: bool,
    }

    impl Fake {
        fn with_bodies(names: &[&str]) -> Self {
            Self { bodies: names.iter().map(|n| (n.to_string(), 1000.0, "Default".to_string())).collect(), ..Default::default() }
        }
        fn density(&self, material: &str) -> f64 {
            match material {
                "Default" => 1.29,
                other => self
                    .materials
                    .iter()
                    .find(|m| m["name"] == other)
                    .and_then(|m| m["density_g_cm3"].as_f64())
                    .unwrap_or(1.0),
            }
        }
    }

    impl Host for Fake {
        type Error = String;
        fn run(&mut self, command: &str, args: &Value) -> Result<Value, Failure<String>> {
            self.calls.push(command.to_string());
            match command {
                "document.inspect" => Ok(json!({ "bodies": self.bodies.iter().map(|b| b.0.clone()).collect::<Vec<_>>() })),
                "material.list" => {
                    let mut all = vec![json!({"name": "Default", "density_g_cm3": 1.29}), json!({"name": "Steel", "density_g_cm3": 7.85})];
                    all.extend(self.materials.iter().cloned());
                    Ok(json!({ "materials": all }))
                }
                "material.create" => {
                    let mut row = args.clone();
                    row["user"] = json!(true);
                    self.materials.push(row);
                    Ok(json!({}))
                }
                "material.assign" => {
                    let names: Vec<String> = serde_json::from_value(args["bodies"].clone()).unwrap();
                    let material = args["material"].as_str().unwrap().to_string();
                    for b in &mut self.bodies {
                        if names.contains(&b.0) {
                            b.2 = material.clone();
                        }
                    }
                    Ok(json!({ "bodies": names, "material": material }))
                }
                "inspect.measure" => {
                    let wanted: Option<Vec<String>> = args.get("bodies").map(|b| serde_json::from_value(b.clone()).unwrap());
                    let rows: Vec<Value> = self
                        .bodies
                        .iter()
                        .filter(|b| wanted.as_ref().is_none_or(|w| w.contains(&b.0)))
                        .map(|b| {
                            let d = self.density(&b.2);
                            json!({"name": b.0, "volume_mm3": b.1, "mass_g": b.1 * d / 1000.0, "density_g_cm3": d, "material": b.2, "closed": true,
                                   "bbox": {"min": [0.0, 0.0, 0.0], "max": [100.0, 10.0, 5.0]}})
                        })
                        .collect();
                    Ok(json!({ "bodies": rows }))
                }
                "material.inspect" => {
                    let name = args["material"].as_str().unwrap();
                    Ok(self.materials.iter().find(|m| m["name"] == name).cloned().unwrap_or(json!({"name": name, "density_g_cm3": 1.29})))
                }
                "model.mesh" => {
                    let (p, tr, _) = crate::stress2d::box_mesh([100.0, 10.0, 5.0]);
                    Ok(json!({ "body": args["body"], "positions": p, "triangles": tr }))
                }
                "model.faces" => {
                    let (_, _, faces) = crate::stress2d::box_mesh([100.0, 10.0, 5.0]);
                    let list: Vec<Value> = faces
                        .iter()
                        .map(|f| {
                            let normal = match &f.geom {
                                crate::stress2d::FaceGeom::Plane { normal } => json!(normal),
                                _ => Value::Null,
                            };
                            json!({"index": f.index, "name": f.name, "centroid": f.centroid, "area": f.area, "plane_normal": normal})
                        })
                        .collect();
                    Ok(json!({ "faces": list }))
                }
                other => Err(Failure::Host(format!("unexpected command {other}"))),
            }
        }
        fn log(&mut self, _: Level, _: &str) {}
        fn progress(&mut self, _: f32, _: &str) {}
        fn cancelled(&mut self) -> bool {
            self.cancel
        }
    }

    fn run(host: &mut Fake, entry: &str, params: Value) -> Result<Value, Failure<String>> {
        execute(host, entry, &params.to_string()).map(|s| serde_json::from_str(&s).unwrap())
    }

    #[test]
    fn materials_lists_the_library_without_app_commands() {
        let mut host = Fake::default();
        let out = run(&mut host, "materials", json!({})).unwrap();
        assert!(out["count"].as_u64().unwrap() >= 10);
        assert!(host.calls.is_empty());
        let metals = run(&mut host, "materials", json!({"family": "metal"})).unwrap();
        assert!(metals["materials"].as_array().unwrap().iter().all(|m| m["family"] == "metal"));
    }

    #[test]
    fn empty_params_mean_defaults() {
        let mut host = Fake::default();
        assert!(execute(&mut host, "materials", "").is_ok());
        assert!(execute(&mut host, "materials", "  ").is_ok());
    }

    #[test]
    fn apply_creates_once_then_reuses() {
        let mut host = Fake::with_bodies(&["Bracket", "Plate"]);
        let out = run(&mut host, "apply-material", json!({"material": "aluminium-6061-t6"})).unwrap();
        assert_eq!(out["created"], true);
        assert_eq!(out["bodies"], json!(["Bracket", "Plate"]));
        assert_eq!(host.bodies[0].2, "Aluminium 6061-T6");

        let again = run(&mut host, "apply-material", json!({"material": "aluminium-6061-t6", "bodies": ["Plate"]})).unwrap();
        assert_eq!(again["created"], false);
        assert_eq!(host.calls.iter().filter(|c| *c == "material.create").count(), 1);
    }

    #[test]
    fn apply_refuses_a_same_named_material_with_other_values() {
        let mut host = Fake::with_bodies(&["Bracket"]);
        host.materials.push(json!({"name": "aluminium 6061-t6", "density_g_cm3": 2.5, "user": true}));
        let err = run(&mut host, "apply-material", json!({"material": "aluminium-6061-t6"})).unwrap_err();
        assert!(matches!(err, Failure::Input { ref field, .. } if field == "material"));
        assert!(!host.calls.contains(&"material.create".to_string()));
        assert!(!host.calls.contains(&"material.assign".to_string()));
    }

    #[test]
    fn apply_checks_input_before_editing() {
        let mut host = Fake::with_bodies(&["Bracket"]);
        for params in [
            json!({"material": "unobtainium"}),
            json!({"material": "steel-s235jr", "bodies": ["Nope"]}),
            json!({"material": "steel-s235jr", "bodies": ["Bracket", "Bracket"]}),
            json!({"material": "steel-s235jr", "bodies": [""]}),
            json!({"material": "steel-s235jr", "bodies": ["a\u{7}b"]}),
            json!({"material": "steel-s235jr", "typo": 1}),
            json!({"bodies": ["Bracket"]}),
            json!({"material": 5}),
        ] {
            let err = run(&mut host, "apply-material", params.clone()).unwrap_err();
            assert!(matches!(err, Failure::Input { .. }), "{params}: {err:?}");
        }
        assert!(!host.calls.iter().any(|c| c.starts_with("material.create") || c == "material.assign"));
    }

    #[test]
    fn apply_on_empty_design_is_refused() {
        let mut host = Fake::default();
        let err = run(&mut host, "apply-material", json!({"material": "steel-s235jr"})).unwrap_err();
        assert!(matches!(err, Failure::Input { .. }));
    }

    #[test]
    fn cancel_stops_before_editing() {
        let mut host = Fake::with_bodies(&["Bracket"]);
        host.cancel = true;
        assert_eq!(run(&mut host, "apply-material", json!({"material": "steel-s235jr"})).unwrap_err(), Failure::Cancelled);
        assert!(!host.calls.contains(&"material.create".to_string()));
    }

    #[test]
    fn report_gives_mass_and_allowable_stress() {
        let mut host = Fake::with_bodies(&["Bracket", "Plate"]);
        run(&mut host, "apply-material", json!({"material": "steel-s235jr", "bodies": ["Bracket"]})).unwrap();
        let out = run(&mut host, "material-report", json!({"safety_factor": 1.5})).unwrap();
        let bracket = &out["bodies"][0];
        assert_eq!(bracket["material"], "Steel S235JR");
        assert_eq!(bracket["library_id"], "steel-s235jr");
        assert!((bracket["mass_g"].as_f64().unwrap() - 7.85).abs() < 1e-9);
        assert!((bracket["allowable_mpa"].as_f64().unwrap() - 156.7).abs() < 1e-9);
        assert_eq!(bracket["allowable_basis"], "yield_strength");
        let plate = &out["bodies"][1];
        assert!(plate["allowable_mpa"].is_null());
        assert!(plate["note"].is_string());
        assert!(out["stress_analysis"].as_str().unwrap().starts_with("not done"));
    }

    #[test]
    fn report_leaves_out_strengths_of_a_changed_same_named_material() {
        let mut host = Fake::with_bodies(&["Bracket"]);
        run(&mut host, "apply-material", json!({"material": "steel-s235jr"})).unwrap();
        host.materials[0]["density_g_cm3"] = json!(7.0);
        let out = run(&mut host, "material-report", json!({})).unwrap();
        assert!(out["bodies"][0]["allowable_mpa"].is_null());
        assert!(out["bodies"][0]["note"].as_str().unwrap().contains("differs"));
    }

    #[test]
    fn report_only_reads() {
        let mut host = Fake::with_bodies(&["Bracket"]);
        run(&mut host, "material-report", json!({})).unwrap();
        assert_eq!(host.calls, ["inspect.measure"]);
    }

    #[test]
    fn report_checks_safety_factor() {
        let mut host = Fake::with_bodies(&["Bracket"]);
        for sf in [json!(0.5), json!(11), json!("2"), json!(null)] {
            assert!(run(&mut host, "material-report", json!({"safety_factor": sf})).is_err(), "{sf}");
        }
    }

    #[test]
    fn unknown_entry_and_bad_json_are_input_errors() {
        let mut host = Fake::default();
        assert!(matches!(execute(&mut host, "optimize", "{}"), Err(Failure::Input { .. })));
        assert!(matches!(execute(&mut host, "materials", "{"), Err(Failure::Input { .. })));
        assert!(matches!(execute(&mut host, "materials", "[]"), Err(Failure::Input { .. })));
    }

    #[test]
    fn host_errors_pass_through_unchanged() {
        struct Broken;
        impl Host for Broken {
            type Error = &'static str;
            fn run(&mut self, _: &str, _: &Value) -> Result<Value, Failure<&'static str>> {
                Err(Failure::Host("denied"))
            }
            fn log(&mut self, _: Level, _: &str) {}
            fn progress(&mut self, _: f32, _: &str) {}
            fn cancelled(&mut self) -> bool {
                false
            }
        }
        assert_eq!(execute(&mut Broken, "material-report", "{}"), Err(Failure::Host("denied")));
    }

    #[test]
    fn benchmark_runs_and_checks_its_input() {
        let mut host = Fake::default();
        let out = run(&mut host, "solver-benchmark", json!({"nx": 20, "ny": 10, "iterations": 5})).unwrap();
        assert_eq!(out["design_iterations"], 5);
        assert_eq!(out["elements"], 200);
        assert!(host.calls.is_empty());
        for bad in [json!({"nx": 1}), json!({"load": "x"}), json!({"nx": 2000, "ny": 2000}), json!({"volfrac": 2})] {
            assert!(matches!(run(&mut host, "solver-benchmark", bad.clone()), Err(Failure::Input { .. })), "{bad}");
        }
    }

    #[test]
    fn stress_check_on_a_cantilever() {
        let mut host = Fake::with_bodies(&["Bar"]);
        run(&mut host, "apply-material", json!({"material": "aluminium-6061-t6", "body": "Bar"})).unwrap();
        let out = run(&mut host, "stress-check", json!({"fixed": "-x", "load": "+x", "force_y": -100, "resolution": 200})).unwrap();
        let p99 = out["p99_von_mises_mpa"].as_f64().unwrap();
        assert!((100.0..=135.0).contains(&p99), "p99 {p99}");
        assert_eq!(out["material"], "Aluminium 6061-T6");
        assert_eq!(out["plane"], "xy");
        assert_eq!(out["strength_basis"], "yield_strength");
        assert!(out["verdict"].as_str().unwrap().contains("safety factor"));
        assert!(host.calls.ends_with(&["inspect.measure".to_string(), "model.mesh".into(), "model.faces".into()]));
    }

    #[test]
    fn stress_check_needs_strength_values_and_one_body() {
        let mut host = Fake::with_bodies(&["Bar"]);
        let err = run(&mut host, "stress-check", json!({"fixed": "-x", "load": "+x"})).unwrap_err();
        assert!(matches!(err, Failure::Input { ref field, ref reason } if field == "material" && reason.contains("apply a library material")), "{err:?}");
        let ok = run(&mut host, "stress-check", json!({"fixed": "-x", "load": "+x", "material": "steel-s235jr"})).unwrap();
        assert_eq!(ok["material_source"], "chosen for this check");

        let mut two = Fake::with_bodies(&["A", "B"]);
        let err = run(&mut two, "stress-check", json!({"fixed": "-x", "load": "+x", "material": "steel-s235jr"})).unwrap_err();
        assert!(matches!(err, Failure::Input { ref reason, .. } if reason.contains("A, B")), "{err:?}");
    }

    #[test]
    fn stress_check_refuses_bad_input_before_reading() {
        let mut host = Fake::with_bodies(&["Bar"]);
        for bad in [
            json!({"fixed": "", "load": "+x"}),
            json!({"fixed": "-x", "load": "+x", "resolution": 5}),
            json!({"fixed": "-x", "load": "+x", "plane": "ab"}),
            json!({"fixed": "-x", "load": "+x", "force_x": 1e20}),
            json!({"fixed": "-x", "load": "+x", "safety_factor": 0}),
            json!({"fixed": "-x", "load": "+x", "unknown": 1}),
        ] {
            assert!(matches!(run(&mut host, "stress-check", bad.clone()), Err(Failure::Input { .. })), "{bad}");
        }
        assert!(host.calls.is_empty());
    }

    #[test]
    fn body_field_and_bodies_list_are_alternatives() {
        let mut host = Fake::with_bodies(&["A", "B"]);
        let out = run(&mut host, "apply-material", json!({"material": "steel-s235jr", "body": "B"})).unwrap();
        assert_eq!(out["bodies"], json!(["B"]));
        assert!(run(&mut host, "apply-material", json!({"material": "steel-s235jr", "body": "B", "bodies": ["A"]})).is_err());
        let report = run(&mut host, "material-report", json!({"body": "A"})).unwrap();
        assert_eq!(report["body_count"], 1);
    }

    #[test]
    fn hostile_params_never_panic() {
        let mut host = Fake::with_bodies(&["Bracket"]);
        let long = "x".repeat(100_000);
        let samples = [
            "null".to_string(),
            "1e999".to_string(),
            "{\"bodies\": null}".to_string(),
            "{\"safety_factor\": 1e308}".to_string(),
            format!("{{\"material\": \"{long}\"}}"),
            format!("{{\"bodies\": [\"{long}\"]}}"),
            "{\"bodies\": [1,2,3]}".to_string(),
            "\u{0}".to_string(),
            "[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[".to_string(),
        ];
        for entry in ["materials", "apply-material", "material-report", "solver-benchmark", "stress-check", "", &long] {
            for s in &samples {
                let _ = execute(&mut host, entry, s);
            }
        }
    }
}
