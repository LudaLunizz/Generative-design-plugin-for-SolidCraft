#!/usr/bin/env python3
"""Write plugin.json from materials/library.json, so the Material dropdown always matches the library.
Run from anywhere: python3 -I tools/gen_manifest.py"""
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
VERSION = "0.2.1"

lib = json.loads((ROOT / "materials" / "library.json").read_text(encoding="utf-8"))
one_of = [{"const": m["id"], "title": m["name"]} for m in lib["materials"]]
bodies = {"type": "array", "title": "Bodies", "items": {"type": "string"}, "default": [],
          "description": "Body names; empty means every body in the design"}
sf = {"type": "number", "title": "Safety factor", "default": 2.0, "minimum": 1.0, "maximum": 10.0,
      "description": "Allowable stress = yield strength (or tensile strength when there is no yield) / safety factor"}
manifest = {
    "id": "gen-design-engine",
    "name": "Generative Design Engine",
    "version": VERSION,
    "description": "Engineering material library with sourced cards (density, stiffness, strength), applied to bodies, "
                   "plus a mass and allowable-stress report. Generative design comes in a later version.",
    "capabilities": ["commands@0.1"],
    "runtime": {
        "module": "plugin.wasm",
        "entries": [
            {"name": "apply-material", "effect": "edit",
             "slots": [{"kind": "ribbon", "name": "Manage"}],
             "input_schema": {"type": "object", "required": ["material"], "properties": {
                 "material": {"type": "string", "title": "Material", "default": "aluminium-6061-t6", "oneOf": one_of,
                              "description": "Creates this material with its engineering card in the design (once) and assigns it"},
                 "bodies": bodies}}},
            {"name": "material-report", "effect": "read",
             "slots": [{"kind": "ribbon", "name": "Inspect"}],
             "input_schema": {"type": "object", "properties": {"bodies": bodies, "safety_factor": sf}}},
            {"name": "materials", "effect": "read",
             "slots": [{"kind": "ribbon", "name": "Manage"}],
             "input_schema": {"type": "object", "properties": {
                 "family": {"type": "string", "title": "Family", "default": "any", "oneOf": [
                     {"const": "any", "title": "All materials"}, {"const": "metal", "title": "Metals"},
                     {"const": "polymer", "title": "Plastics"}]},
                 "safety_factor": sf}}},
            {"name": "solver-benchmark", "effect": "read",
             "input_schema": {"type": "object", "properties": {
                 "nx": {"type": "integer", "title": "Cells along X", "default": 60, "minimum": 2, "maximum": 2000},
                 "ny": {"type": "integer", "title": "Cells along Y", "default": 30, "minimum": 2, "maximum": 2000},
                 "iterations": {"type": "integer", "title": "Design iterations", "default": 50, "minimum": 1, "maximum": 1000},
                 "load": {"type": "string", "title": "Test case", "default": "cantilever", "oneOf": [
                     {"const": "cantilever", "title": "Cantilever (left edge fixed, load at right middle)"},
                     {"const": "mbb", "title": "MBB beam (half, symmetric)"}]},
                 "volfrac": {"type": "number", "title": "Material share", "default": 0.5, "minimum": 0.05, "maximum": 0.95},
                 "rmin": {"type": "number", "title": "Filter radius (cells)", "default": 1.5, "minimum": 1, "maximum": 10},
                 "cg_tol": {"type": "number", "title": "Solver tolerance", "default": 1e-6, "minimum": 1e-12, "maximum": 0.1},
                 "cg_max": {"type": "integer", "title": "Solver step limit", "default": 20000, "minimum": 1, "maximum": 100000},
                 "picture": {"type": "boolean", "title": "Text picture of the result", "default": True}}}},
        ],
        "commands": sorted(["document.inspect", "inspect.measure", "material.assign", "material.create",
                            "material.list"]),
        "limits": {"fuel": 100000000000, "memory_mib": 64, "deadline_ms": 30000},
    },
}
text = json.dumps(manifest, indent=2, ensure_ascii=False) + "\n"
if len(text.encode()) > 16 * 1024:
    sys.exit("plugin.json would exceed SolidCraft's 16 KiB manifest cap")
if len(manifest["description"]) > 240:
    sys.exit("description over 240 characters")
(ROOT / "plugin.json").write_text(text, encoding="utf-8")
print(f"plugin.json: {len(text.encode())} bytes, {len(one_of)} materials")
