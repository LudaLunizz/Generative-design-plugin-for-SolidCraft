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
}

fn apply_material<H: Host>(host: &mut H, library: &Library, p: ApplyParams) -> Result<Value, Failure<H::Error>> {
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
    #[serde(default = "default_sf")]
    safety_factor: f64,
}

fn default_sf() -> f64 {
    2.0
}

fn material_report<H: Host>(host: &mut H, library: &Library, p: ReportParams) -> Result<Value, Failure<H::Error>> {
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
                            json!({"name": b.0, "volume_mm3": b.1, "mass_g": b.1 * d / 1000.0, "density_g_cm3": d, "material": b.2, "closed": true})
                        })
                        .collect();
                    Ok(json!({ "bodies": rows }))
                }
                "material.inspect" => {
                    let name = args["material"].as_str().unwrap();
                    Ok(self.materials.iter().find(|m| m["name"] == name).cloned().unwrap_or(json!({"name": name, "density_g_cm3": 1.29})))
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
        for entry in ["materials", "apply-material", "material-report", "", &long] {
            for s in &samples {
                let _ = execute(&mut host, entry, s);
            }
        }
    }
}
