//! The built-in material library: `materials/library.json`, compiled into the plugin.
//!
//! Each material becomes a SolidCraft physical material card (`material.create`'s
//! `physical_card`, schema_version 1, SI units). The types mirror that public command schema so
//! a bad entry fails the tests here instead of failing in the app.

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet};

const LIBRARY_JSON: &str = include_str!("../materials/library.json");

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Library {
    pub schema: u32,
    pub sources: BTreeMap<String, Source>,
    pub materials: Vec<Material>,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Source {
    pub title: String,
    pub kind: SourceKind,
    #[serde(default)]
    pub license: Option<String>,
    #[serde(default)]
    pub reference: Option<String>,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum SourceKind {
    UserSupplied,
    Measurement,
    PublicReference,
    Illustrative,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Family {
    Wood,
    Polymer,
    Glass,
    Metal,
    Other,
}

impl Family {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Wood => "wood",
            Self::Polymer => "polymer",
            Self::Glass => "glass",
            Self::Metal => "metal",
            Self::Other => "other",
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Material {
    pub id: String,
    pub name: String,
    pub color: String,
    pub family: Family,
    pub designation: String,
    pub notes: String,
    pub properties: Vec<Property>,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Property {
    pub source_id: String,
    pub data: PropertyData,
    #[serde(default)]
    pub conditions: Conditions,
}

#[derive(Clone, Debug, Default, Deserialize, Serialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct Conditions {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub temperature_k: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub moisture_dry_basis_kg_kg: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strain_rate_per_s: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub state: Option<String>,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Axes {
    Material123,
    WoodLrt,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum PropertyData {
    Density { kg_m3: f64 },
    Elasticity { model: Elasticity },
    TensileStrength { pa: f64 },
    CompressiveStrength { pa: f64 },
    YieldStrength { pa: f64 },
    ShearStrength { pa: f64 },
    ThermalConductivity { w_m_k: [f64; 3], axes: Axes },
    ThermalExpansion { per_k: [f64; 3], axes: Axes },
    SpecificHeat { j_kg_k: f64 },
}

impl PropertyData {
    pub fn kind(&self) -> &'static str {
        match self {
            Self::Density { .. } => "density",
            Self::Elasticity { .. } => "elasticity",
            Self::TensileStrength { .. } => "tensile_strength",
            Self::CompressiveStrength { .. } => "compressive_strength",
            Self::YieldStrength { .. } => "yield_strength",
            Self::ShearStrength { .. } => "shear_strength",
            Self::ThermalConductivity { .. } => "thermal_conductivity",
            Self::ThermalExpansion { .. } => "thermal_expansion",
            Self::SpecificHeat { .. } => "specific_heat",
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
#[serde(tag = "symmetry", rename_all = "snake_case", deny_unknown_fields)]
pub enum Elasticity {
    Isotropic { young_modulus_pa: f64, poisson_ratio: f64 },
}

/// The numbers a design check needs, in engineering units.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Strengths {
    pub density_kg_m3: Option<f64>,
    pub young_modulus_gpa: Option<f64>,
    pub poisson_ratio: Option<f64>,
    pub yield_mpa: Option<f64>,
    pub tensile_mpa: Option<f64>,
    pub shear_mpa: Option<f64>,
}

impl Strengths {
    /// The stress a design may reach under `safety_factor`, and what it is based on: yield when
    /// the material has one, else tensile strength (brittle or printed materials).
    pub fn allowable_mpa(&self, safety_factor: f64) -> Option<(f64, &'static str)> {
        if !safety_factor.is_finite() || safety_factor <= 0.0 {
            return None;
        }
        self.yield_mpa
            .map(|y| (y / safety_factor, "yield_strength"))
            .or_else(|| self.tensile_mpa.map(|t| (t / safety_factor, "tensile_strength")))
    }
}

impl Material {
    pub fn density_kg_m3(&self) -> Option<f64> {
        self.properties.iter().find_map(|p| match p.data {
            PropertyData::Density { kg_m3 } => Some(kg_m3),
            _ => None,
        })
    }

    pub fn strengths(&self) -> Strengths {
        let mut s = Strengths::default();
        for property in &self.properties {
            match &property.data {
                PropertyData::Density { kg_m3 } => s.density_kg_m3 = Some(*kg_m3),
                PropertyData::Elasticity { model: Elasticity::Isotropic { young_modulus_pa, poisson_ratio } } => {
                    s.young_modulus_gpa = Some(young_modulus_pa / 1e9);
                    s.poisson_ratio = Some(*poisson_ratio);
                }
                PropertyData::YieldStrength { pa } => s.yield_mpa = Some(pa / 1e6),
                PropertyData::TensileStrength { pa } => s.tensile_mpa = Some(pa / 1e6),
                PropertyData::ShearStrength { pa } => s.shear_mpa = Some(pa / 1e6),
                _ => {}
            }
        }
        s
    }

    /// The `physical_card` for `material.create`: schema_version 1 with only the sources this
    /// material cites, in the order they are first cited.
    pub fn card(&self, library: &Library) -> Value {
        let mut seen = BTreeSet::new();
        let mut sources = Vec::new();
        for property in &self.properties {
            if seen.insert(property.source_id.as_str())
                && let Some(source) = library.sources.get(&property.source_id)
            {
                sources.push(json!({
                    "id": property.source_id,
                    "title": source.title,
                    "kind": source.kind,
                    "license": source.license,
                    "reference": source.reference,
                }));
            }
        }
        json!({
            "schema_version": 1,
            "designation": self.designation,
            "family": self.family,
            "sources": sources,
            "properties": self.properties,
        })
    }

    /// The `material.create` parameters for this material.
    pub fn create_params(&self, library: &Library) -> Option<Value> {
        let density = self.density_kg_m3()?;
        Some(json!({
            "name": self.name,
            "density_g_cm3": density / 1000.0,
            "color": self.color,
            "physical_card": self.card(library),
        }))
    }
}

impl Library {
    /// The library compiled into the plugin.
    pub fn builtin() -> Result<Self, String> {
        let library: Self = serde_json::from_str(LIBRARY_JSON).map_err(|e| format!("the built-in material library does not parse: {e}"))?;
        library.validate()?;
        Ok(library)
    }

    pub fn find(&self, id: &str) -> Option<&Material> {
        self.materials.iter().find(|m| m.id == id)
    }

    /// The same checks SolidCraft makes on a card, plus the library's own rules, so a bad entry is
    /// caught before it reaches the app.
    pub fn validate(&self) -> Result<(), String> {
        if self.schema != 1 {
            return Err(format!("library schema {} is not supported (expected 1)", self.schema));
        }
        for (id, source) in &self.sources {
            text(id, 128, "source id")?;
            text(&source.title, 256, &format!("source {id} title"))?;
            if let Some(reference) = &source.reference {
                text(reference, 2048, &format!("source {id} reference"))?;
            }
            if let Some(license) = &source.license {
                text(license, 512, &format!("source {id} license"))?;
            }
        }
        let mut ids = BTreeSet::new();
        let mut names = BTreeSet::new();
        for m in &self.materials {
            let at = |what: &str| format!("material {}: {what}", m.id);
            if m.id.is_empty() || m.id.len() > 64 || !m.id.bytes().all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == b'-') {
                return Err(at("id must be 1-64 lowercase letters, digits or hyphens"));
            }
            if !ids.insert(m.id.as_str()) {
                return Err(at("duplicate id"));
            }
            text(&m.name, 128, &at("name"))?;
            if !names.insert(m.name.to_lowercase()) {
                return Err(at("duplicate name"));
            }
            if BUILTIN_NAMES.iter().any(|b| b.eq_ignore_ascii_case(m.name.trim())) {
                return Err(at("name equals one of SolidCraft's built-in materials"));
            }
            text(&m.designation, 256, &at("designation"))?;
            text(&m.notes, 1024, &at("notes"))?;
            if !valid_color(&m.color) {
                return Err(at("color must be #rrggbb"));
            }
            if m.properties.is_empty() || m.properties.len() > 16 {
                return Err(at("needs 1-16 properties"));
            }
            let mut kinds = BTreeSet::new();
            for p in &m.properties {
                if !kinds.insert(p.data.kind()) {
                    return Err(at(&format!("duplicate property {}", p.data.kind())));
                }
                if !self.sources.contains_key(&p.source_id) {
                    return Err(at(&format!("{} cites unknown source {}", p.data.kind(), p.source_id)));
                }
                check_property(&p.data).map_err(|e| at(&e))?;
                check_conditions(&p.conditions).map_err(|e| at(&e))?;
            }
            let density = m.density_kg_m3().ok_or_else(|| at("needs a density"))?;
            // material.create takes 0 < density_g_cm3 <= 1000.
            if density / 1000.0 > 1000.0 {
                return Err(at("density above 1000 g/cm3"));
            }
            let s = m.strengths();
            if s.yield_mpa.is_none() && s.tensile_mpa.is_none() {
                return Err(at("needs a yield or tensile strength"));
            }
            if let (Some(y), Some(t)) = (s.yield_mpa, s.tensile_mpa)
                && y > t
            {
                return Err(at("yield strength above tensile strength"));
            }
        }
        Ok(())
    }
}

/// SolidCraft's built-in material names (0.42): a design material may not reuse them.
pub const BUILTIN_NAMES: &[&str] = &[
    "Default",
    "Steel",
    "Stainless Steel",
    "Aluminum",
    "Brass",
    "Copper",
    "Titanium",
    "Cast Iron",
    "ABS Plastic",
    "PLA",
    "Nylon",
    "Polycarbonate",
    "Oak",
    "Glass",
];

fn text(value: &str, limit: usize, field: &str) -> Result<(), String> {
    if value.trim().is_empty() || value.len() > limit || value.chars().any(char::is_control) {
        return Err(format!("{field}: expected nonempty text of at most {limit} bytes without control characters"));
    }
    Ok(())
}

fn valid_color(c: &str) -> bool {
    c.len() == 7 && c.starts_with('#') && c.bytes().skip(1).all(|b| b.is_ascii_hexdigit())
}

fn positive(v: f64, high: f64, field: &str) -> Result<(), String> {
    if !v.is_finite() || v <= 0.0 || v > high {
        return Err(format!("{field}: expected a finite value in (0, {high}]"));
    }
    Ok(())
}

fn check_property(data: &PropertyData) -> Result<(), String> {
    match data {
        PropertyData::Density { kg_m3 } => positive(*kg_m3, 1e6, "density"),
        PropertyData::Elasticity { model: Elasticity::Isotropic { young_modulus_pa, poisson_ratio } } => {
            positive(*young_modulus_pa, 1e15, "young_modulus_pa")?;
            if !poisson_ratio.is_finite() || *poisson_ratio <= -1.0 || *poisson_ratio >= 0.5 {
                return Err("poisson_ratio: isotropic stability requires -1 < nu < 0.5".into());
            }
            Ok(())
        }
        PropertyData::TensileStrength { pa }
        | PropertyData::CompressiveStrength { pa }
        | PropertyData::YieldStrength { pa }
        | PropertyData::ShearStrength { pa } => positive(*pa, 1e15, data.kind()),
        PropertyData::ThermalConductivity { w_m_k, .. } => w_m_k.iter().try_for_each(|k| positive(*k, 1e9, "thermal_conductivity")),
        PropertyData::ThermalExpansion { per_k, .. } => {
            if per_k.iter().all(|a| a.is_finite() && (-1.0..=1.0).contains(a)) {
                Ok(())
            } else {
                Err("thermal_expansion: expected values in [-1, 1] per K".into())
            }
        }
        PropertyData::SpecificHeat { j_kg_k } => positive(*j_kg_k, 1e9, "specific_heat"),
    }
}

fn check_conditions(c: &Conditions) -> Result<(), String> {
    if let Some(t) = c.temperature_k {
        positive(t, 1e5, "temperature_k")?;
    }
    if let Some(m) = c.moisture_dry_basis_kg_kg
        && (!m.is_finite() || m < 0.0)
    {
        return Err("moisture_dry_basis_kg_kg: expected a finite value >= 0".into());
    }
    if let Some(r) = c.strain_rate_per_s {
        positive(r, 1e12, "strain_rate_per_s")?;
    }
    if let Some(s) = &c.state {
        text(s, 128, "state")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
    use super::*;

    #[test]
    fn builtin_library_parses_and_validates() {
        let lib = Library::builtin().unwrap();
        assert!(lib.materials.len() >= 10);
    }

    #[test]
    fn every_material_has_density_and_a_strength() {
        let lib = Library::builtin().unwrap();
        for m in &lib.materials {
            let s = m.strengths();
            assert!(s.density_kg_m3.is_some(), "{}", m.id);
            assert!(m.strengths().allowable_mpa(2.0).is_some(), "{}", m.id);
        }
    }

    #[test]
    fn every_cited_source_exists_and_every_source_is_cited() {
        let lib = Library::builtin().unwrap();
        let cited: BTreeSet<_> = lib.materials.iter().flat_map(|m| m.properties.iter().map(|p| p.source_id.clone())).collect();
        let declared: BTreeSet<_> = lib.sources.keys().cloned().collect();
        assert_eq!(cited, declared);
    }

    #[test]
    fn card_lists_only_cited_sources() {
        let lib = Library::builtin().unwrap();
        let m = lib.find("steel-s235jr").unwrap();
        let card = m.card(&lib);
        let ids: Vec<_> = card["sources"].as_array().unwrap().iter().map(|s| s["id"].as_str().unwrap()).collect();
        assert_eq!(ids, ["eurocode-steel", "en10025-2"]);
        assert_eq!(card["schema_version"], 1);
    }

    #[test]
    fn create_params_use_g_per_cm3() {
        let lib = Library::builtin().unwrap();
        let p = lib.find("aluminium-6061-t6").unwrap().create_params(&lib).unwrap();
        assert!((p["density_g_cm3"].as_f64().unwrap() - 2.7).abs() < 1e-12);
    }

    #[test]
    fn allowable_prefers_yield_then_tensile() {
        let lib = Library::builtin().unwrap();
        let (a, basis) = lib.find("aluminium-6061-t6").unwrap().strengths().allowable_mpa(2.0).unwrap();
        assert_eq!((a, basis), (138.0, "yield_strength"));
        let (a, basis) = lib.find("pla-printed").unwrap().strengths().allowable_mpa(2.0).unwrap();
        assert_eq!((a, basis), (25.5, "tensile_strength"));
        assert!(Strengths::default().allowable_mpa(f64::NAN).is_none());
    }

    #[test]
    fn validation_catches_bad_entries() {
        let mut lib = Library::builtin().unwrap();
        lib.materials[0].properties[1].data =
            PropertyData::Elasticity { model: Elasticity::Isotropic { young_modulus_pa: 1e9, poisson_ratio: 0.5 } };
        assert!(lib.validate().is_err());

        let mut lib = Library::builtin().unwrap();
        lib.materials[0].name = "steel".into();
        assert!(lib.validate().unwrap_err().contains("built-in"));

        let mut lib = Library::builtin().unwrap();
        lib.materials[1].id = lib.materials[0].id.clone();
        assert!(lib.validate().unwrap_err().contains("duplicate id"));
    }
}
