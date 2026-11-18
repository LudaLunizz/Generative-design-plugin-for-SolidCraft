#!/usr/bin/env python3
"""Turn beam_opt's design.json into a SolidCraft command script (batch / solidcraft-cli run).

usage: python3 -I tools/beam_to_script.py design.json script.json [body-name]

The side view is sketched on the plane with sketch x = world X (along the beam) and sketch
y = world Z (up), and extruded symmetrically about Y = 0 to the plate thickness.
"""
import json, sys

d = json.load(open(sys.argv[1]))
name = sys.argv[3] if len(sys.argv) > 3 else "Beam"
cmds = [{"command": "sketch.create", "params": {"plane": {"origin": [0, 0, 0], "x_dir": [1, 0, 0], "y_dir": [0, 0, 1]}, "name": name + " side view"}}]
for loop in d["loops"]:
    cmds.append({"command": "sketch.line", "params": {"points": [[round(x, 2), round(y, 2)] for x, y in loop], "closed": True}})
cmds.append({"command": "sketch.finish", "params": {}})
cmds.append({"command": "solid.extrude", "params": {"distance": d["thickness"], "direction": "symmetric", "name": name, "body_name": name}})
json.dump(cmds, open(sys.argv[2], "w"))
print(len(cmds), "steps,", sum(len(l) for l in d["loops"]), "vertices")
