#!/usr/bin/env python3
"""Write plugin.json from materials/library.json, so the Material dropdown always matches the library.
Run from anywhere: python3 -I tools/gen_manifest.py"""
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
VERSION = "0.4.1"

lib = json.loads((ROOT / "materials" / "library.json").read_text(encoding="utf-8"))
one_of = [{"const": m["id"], "title": m["name"]} for m in lib["materials"]]
body = {"type": "string", "title": "Body", "default": "",
        "description": "A body name; empty means every body in the design (scripts may pass a bodies list instead)"}
sf = {"type": "number", "title": "Safety factor", "default": 2.0, "minimum": 1.0, "maximum": 10.0,
      "description": "Allowable stress = yield strength (or tensile strength when there is no yield) / safety factor"}
manifest = {
    "id": "gen-design-engine",
    "name": "Generative Design Engine",
    "version": VERSION,
    "description": "Engineering material library (sourced cards), material report, a 2D stress check of flat parts "
                   "with safety factor, and a plate optimiser (sketch + extrude). Third-party plugin, not affiliated with SolidCraft.",
    "capabilities": ["commands@0.1"],
    "runtime": {
        "module": "plugin.wasm",
        "entries": [
            {"name": "apply-material", "effect": "edit",
             "slots": [{"kind": "ribbon", "name": "Manage"}],
             "input_schema": {"type": "object", "required": ["material"], "properties": {
                 "material": {"type": "string", "title": "Material", "default": "aluminium-6061-t6", "oneOf": one_of,
                              "description": "Creates this material with its engineering card in the design (once) and assigns it"},
                 "body": body}}},
            {"name": "material-report", "effect": "read",
             "slots": [{"kind": "ribbon", "name": "Inspect"}],
             "input_schema": {"type": "object", "properties": {"body": body, "safety_factor": sf}}},
            {"name": "materials", "effect": "read",
             "slots": [{"kind": "ribbon", "name": "Manage"}],
             "input_schema": {"type": "object", "properties": {
                 "family": {"type": "string", "title": "Family", "default": "any", "oneOf": [
                     {"const": "any", "title": "All materials"}, {"const": "metal", "title": "Metals"},
                     {"const": "polymer", "title": "Plastics"}]},
                 "safety_factor": sf}}},
            {"name": "stress-check", "effect": "read",
             "slots": [{"kind": "ribbon", "name": "Inspect"}],
             "input_schema": {"type": "object", "required": ["fixed", "load"], "properties": {
                 "fixed": {"type": "string", "title": "Fixed faces",
                           "description": "Where the part is held: a face name (F1:-x, F2:bore:side:0), its number, or a side (-x, +y ...); several with commas"},
                 "load": {"type": "string", "title": "Load faces", "description": "Where the force acts, as for Fixed faces; the force is spread evenly over them"},
                 "force_x": {"type": "number", "title": "Force X", "unit": "N", "default": 0},
                 "force_y": {"type": "number", "title": "Force Y", "unit": "N", "default": -100},
                 "force_z": {"type": "number", "title": "Force Z", "unit": "N", "default": 0},
                 "body": {"type": "string", "title": "Body", "default": "", "description": "The body to check; empty when the design has one body"},
                 "material": {"type": "string", "title": "Material", "default": "body", "oneOf":
                     [{"const": "body", "title": "The body's own material"}] + one_of},
                 "safety_factor": {"type": "number", "title": "Safety factor needed", "default": 2.0, "minimum": 1.0, "maximum": 10.0},
                 "plane": {"type": "string", "title": "Plate plane", "default": "auto", "oneOf": [
                     {"const": "auto", "title": "Automatic (across the thinnest direction)"},
                     {"const": "xy", "title": "XY (thickness along Z)"}, {"const": "xz", "title": "XZ (thickness along Y)"},
                     {"const": "yz", "title": "YZ (thickness along X)"}]},
                 "resolution": {"type": "integer", "title": "Cells along the longer side", "default": 240, "minimum": 10, "maximum": 400},
                 "picture": {"type": "boolean", "title": "Text heat map", "default": True}}}},
            {"name": "optimise-plate", "effect": "edit",
             "slots": [{"kind": "ribbon", "name": "Manage"}],
             "input_schema": {"type": "object", "properties": {
                 "width_mm": {"type": "number", "title": "Width", "unit": "mm", "default": 200, "minimum": 1, "maximum": 10000,
                              "description": "Along X; the fixed edge is at X = 0"},
                 "height_mm": {"type": "number", "title": "Height", "unit": "mm", "default": 100, "minimum": 1, "maximum": 10000,
                               "description": "Along Y; rounded to whole grid cells"},
                 "thickness_mm": {"type": "number", "title": "Thickness", "unit": "mm", "default": 10, "minimum": 0.1, "maximum": 1000},
                 "material": {"type": "string", "title": "Material", "default": "aluminium-6061-t6", "oneOf": one_of},
                 "load_case": {"type": "string", "title": "Load case", "default": "cantilever", "oneOf": [
                     {"const": "cantilever", "title": "Cantilever (left edge fixed, load at right middle)"},
                     {"const": "mbb", "title": "MBB half beam (symmetric, load at top left)"}]},
                 "force_n": {"type": "number", "title": "Force", "unit": "N", "default": 500, "minimum": 0.001, "maximum": 1000000000,
                             "description": "Used for the deflection estimate only; the shape does not depend on it"},
                 "volume_share": {"type": "number", "title": "Material share", "default": 0.4, "minimum": 0.1, "maximum": 0.9},
                 "cells": {"type": "integer", "title": "Cells across the width", "default": 80, "minimum": 20, "maximum": 120},
                 "iterations": {"type": "integer", "title": "Design iterations", "default": 60, "minimum": 10, "maximum": 150},
                 "filter_radius": {"type": "number", "title": "Smallest member (cells)", "default": 2, "minimum": 1, "maximum": 4},
                 "name": {"type": "string", "title": "Body name", "default": "Optimised plate"}}}},
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
                 "solver": {"type": "string", "title": "Solver", "default": "multigrid", "oneOf": [
                     {"const": "multigrid", "title": "Multigrid (fast)"}, {"const": "jacobi", "title": "Jacobi (simple, for comparison)"}]},
                 "picture": {"type": "boolean", "title": "Text picture of the result", "default": True}}}},
        ],
        "commands": sorted(["document.inspect", "inspect.measure", "material.assign", "material.create",
                            "material.list", "model.faces", "model.mesh", "sketch.create", "sketch.finish",
                            "sketch.line", "solid.extrude"]),
        "limits": {"fuel": 100000000000, "memory_mib": 64, "deadline_ms": 30000},
    },
}
# One line per entry (compact) keeps the file readable and well under SolidCraft's 16 KiB cap;
# the material dropdown is repeated in three entries.
entries = manifest["runtime"]["entries"]
manifest["runtime"]["entries"] = [f"@@{k}@@" for k in range(len(entries))]
text = json.dumps(manifest, indent=2, ensure_ascii=False)
for k, e in enumerate(entries):
    text = text.replace(json.dumps(f"@@{k}@@"), json.dumps(e, separators=(",", ":"), ensure_ascii=False))
text += "\n"
manifest["runtime"]["entries"] = entries
json.loads(text)  # must still be valid JSON
if len(text.encode()) > 16 * 1024:
    sys.exit("plugin.json would exceed SolidCraft's 16 KiB manifest cap")
if len(manifest["description"]) > 240:
    sys.exit("description over 240 characters")
(ROOT / "plugin.json").write_text(text, encoding="utf-8")
print(f"plugin.json: {len(text.encode())} bytes, {len(one_of)} materials")
