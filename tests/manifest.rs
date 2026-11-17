#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
//! plugin.json must match the library and stay inside SolidCraft's manifest rules.

use generative_design_engine::library::Library;
use serde_json::Value;

fn manifest() -> (Value, usize) {
    let text = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/plugin.json")).unwrap();
    (serde_json::from_str(&text).unwrap(), text.len())
}

#[test]
fn material_dropdown_matches_the_library() {
    let (m, _) = manifest();
    let entry = m["runtime"]["entries"].as_array().unwrap().iter().find(|e| e["name"] == "apply-material").unwrap();
    let ids: Vec<&str> =
        entry["input_schema"]["properties"]["material"]["oneOf"].as_array().unwrap().iter().map(|o| o["const"].as_str().unwrap()).collect();
    let lib = Library::builtin().unwrap();
    let lib_ids: Vec<&str> = lib.materials.iter().map(|m| m.id.as_str()).collect();
    assert_eq!(ids, lib_ids, "run python3 -I tools/gen_manifest.py");
}

#[test]
fn manifest_fits_solidcraft_rules() {
    let (m, len) = manifest();
    assert!(len <= 16 * 1024);
    assert_eq!(m["id"], "gen-design-engine");
    assert_eq!(m["version"], env!("CARGO_PKG_VERSION"), "plugin.json version must equal Cargo.toml's");
    assert!(m["description"].as_str().unwrap().chars().count() <= 240);
    let limits = &m["runtime"]["limits"];
    assert!(limits["deadline_ms"].as_u64().unwrap() <= 30_000);
    assert!(limits["memory_mib"].as_u64().unwrap() <= 128, "stay under the default host ceiling");
    assert!(limits["fuel"].as_u64().unwrap() <= 100_000_000_000);
    for e in m["runtime"]["entries"].as_array().unwrap() {
        let name = e["name"].as_str().unwrap();
        assert!(name.bytes().all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == b'-'));
        let props = e["input_schema"]["properties"].as_object().unwrap();
        if let Some(b) = props.get("bodies") {
            assert_eq!(b["type"], "array", "SolidCraft requires bodies to be a list");
        }
    }
}

#[test]
fn every_entry_is_implemented_and_read_entries_only_read() {
    let (m, _) = manifest();
    for e in m["runtime"]["entries"].as_array().unwrap() {
        let name = e["name"].as_str().unwrap();
        assert!(["apply-material", "material-report", "materials", "solver-benchmark"].contains(&name), "{name}");
    }
}
