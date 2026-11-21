# Wall bracket runs: the data behind the calculations

These are the runs of `eso3d`, the 3D "start with a solid and take out what is not needed" tool, on one
job: a 3D-printed (PETG) wall bracket that carries **50 N down and 15 N sideways** on a small pad at its tip,
150 mm out from the wall. The start is a solid frustum (80 x 100 mm at the wall, 20 x 20 mm at the tip); the
tool removes the material that carries the least stress, round by round, and checks the result with the 3D solver.
Every picture below is the solver's own output, drawn by the tool's own renderer.

> **A design aid, not a certified analysis.** The numbers come from a voxel finite-element model with
> library material values (PETG, 47 MPa allowable) and can be wrong. Nobody should rely on them for a
> real mounting without a qualified check. This project is third party and not affiliated with SolidCraft
> (see the main [README](../../README.md)).

## How to read the pictures

- Colour is **von Mises stress**: blue = 0, green, yellow, red = the **peak of that picture**. The bar under
  each picture runs from 0 to that peak, so **colours are not comparable between pictures**: use the numbers
  in the table (a red picture with a low peak is a calmer part than a green picture with a high one).
- The peak is taken *away from the load pad and the supports*, where a point load would give a meaningless spike.
- Each run has two load cases (**down** 50 N and **side** 15 N) in three views (isometric, from the side,
  from above). The checks use 1.25 mm voxels.
- "SF" is the safety factor against 47 MPa. From run 3 on the target was SF 4: a run stops, or a bar grows, when it would drop below.

## The runs

| Run | What changed | Mass | SF down / side | Tip moves down / side |
|---|---|---|---|---|
| 1 | first design: voxel staircase from the pyramid | 78.5 g (61.8 cm³) | 17.4 / 14.4 (2.5 mm grid) | 0.58 / 1.30 mm |
| 2 | run 1 smoothed: **failed** the finer check, thin diagonals vanished in the blur | | fails: 0.3 (side) | 843 mm down (load path broken) |
| 3 | minimum member 5 mm; smoothing that keeps the plate and pad exact | 92.5 g (72.9 cm³) | 5.5 / 8.0 | 1.3 / 3.8 mm |
| 4 | run 3 + two 4.5 mm countersunk screw holes | ~92 g (72.6 cm³) | 5.5 / 8.0 | 1.3 / 3.8 mm |
| 5 | run 3 + holes for 5 mm screws in 6 mm wall plugs (5.5 / 10.4 mm) | ~92 g (72.5 cm³) | 5.5 / 8.0 | 1.3 / 3.8 mm |
| [6](#run-6) | **round bars**: the 13 pillars of run 5 rebuilt as round bars with filleted joints | 93.7 g (73.8 cm³) | 7.6 / 15.0 | 0.92 / 0.67 mm |
| [7](#run-7) | two-way bar sizing: thinner where the stress is low, thicker where it is high | 74.3 g (58.5 cm³) | 4.3 / 6.8 | 2.4 / 1.6 mm |
| [8](#run-8) | run 7 + a slim frame instead of the full wall plate, **held only at the screws** | 57.9 g (45.6 cm³) | 4.1 / 6.8 | 6.0 / 3.3 mm |
| [9](#run-9) | run 8 + tapered bars | 57.3 g (45.1 cm³) | 4.1 / 5.7 | 7.8 / 3.5 mm |
| 10 | run 8 + the wall pushing back (`hold=contact`): **interim, still running**. Round 0 (before any bar was resized): 46.8 cm³, strips 12.0 mm | ~59 g (46.8 cm³) | 4.4 / 8.0 | 2.6 / 1.2 mm |

Run 10 so far (checked 21 Nov, 02:38, round 2 of up to 10): round 0 passes with SF 4.4 / 8.0 at about the
mass of run 8 and a tip that moves less than half as far (2.6 against 6.0 mm), because the wall now pushes back
where the plate presses. Round 1 went lighter (32.7 cm³) and dropped under the target on the down load (SF 3.5), so
the run is not finished and is **not a result yet**; the final design and its pictures follow when it ends.

Things to know when comparing:

- Runs 6 and 7 are held on the **whole wall face**; runs 8 and 9 only **at the two screws** with no wall
  pushing back, which is the safe side but makes the tip move more. Like for like, run 6 held at the screws
  only gives SF 4.7 / 9.3 and a tip movement of 4.3 / 2.4 mm at 94 g, against run 8's 58 g at SF 4.1 / 6.8.
- Run 8 stays the best result so far: tapering (run 9) gained about 1 %, because 11 of its 13 bars already
  sit at the 5 mm minimum, so there is nothing to taper.
- Runs 1 to 5 only have small voxel previews, so the gallery starts at run 6, where the surface is smooth
  and the stress is drawn on the real shape.

The full text output of every run (what was run, the per-bar sizes and stresses, the screw loads) is in
[`data/`](data/); the command that made a run is on its first lines.

## Run 6

Round bars with 3 mm fillets, bends rounded into arcs, hole rims rounded. 73.8 cm³, closed mesh of 489k triangles.

| Isometric | From the side | From above |
|---|---|---|
| ![run 6 down, isometric](run6-down-iso.jpg) | ![run 6 down, side view](run6-down-side.jpg) | ![run 6 down, from above](run6-down-top.jpg) |
| ![run 6 side load, isometric](run6-side-iso.jpg) | ![run 6 side load, side view](run6-side-side.jpg) | ![run 6 side load, from above](run6-side-top.jpg) |

Top row: 50 N down (peak 6.20 MPa). Bottom row: 15 N sideways (peak 3.12 MPa).

## Run 7

Every bar sized both ways until it just passes SF 4; no bar thinner than 5 mm. 58.5 cm³, a fifth lighter than run 6.

| Isometric | From the side | From above |
|---|---|---|
| ![run 7 down, isometric](run7-down-iso.jpg) | ![run 7 down, side view](run7-down-side.jpg) | ![run 7 down, from above](run7-down-top.jpg) |
| ![run 7 side load, isometric](run7-side-iso.jpg) | ![run 7 side load, side view](run7-side-side.jpg) | ![run 7 side load, from above](run7-side-top.jpg) |

Top row: 50 N down (peak 10.88 MPa). Bottom row: 15 N sideways (peak 6.93 MPa).

## Run 8

The full wall plate becomes a 5 mm frame: pads round the screws and the bar roots, strips between them (21.6 mm
wide, sized too). 45.6 cm³ / 57.9 g. The hot spots next to the screw pads are the plate bending between the
screws (partly an artefact of how sharp the model makes the edge of the clamped zone).

| Isometric | From the side | From above |
|---|---|---|
| ![run 8 down, isometric](run8-down-iso.jpg) | ![run 8 down, side view](run8-down-side.jpg) | ![run 8 down, from above](run8-down-top.jpg) |
| ![run 8 side load, isometric](run8-side-iso.jpg) | ![run 8 side load, side view](run8-side-side.jpg) | ![run 8 side load, from above](run8-side-top.jpg) |

Top row: 50 N down (peak 11.40 MPa). Bottom row: 15 N sideways (peak 6.95 MPa).

## Run 9

Run 8 with every point along each bar sized on its own. Only two bars taper (bar 2 from 6.8 to 5.0 mm, bar 7
from 5.0 to 6.4 mm); the gain is about 1 % and the tip moves more, so run 8 remains the better one.

| Isometric | From the side | From above |
|---|---|---|
| ![run 9 down, isometric](run9-down-iso.jpg) | ![run 9 down, side view](run9-down-side.jpg) | ![run 9 down, from above](run9-down-top.jpg) |
| ![run 9 side load, isometric](run9-side-iso.jpg) | ![run 9 side load, side view](run9-side-side.jpg) | ![run 9 side load, from above](run9-side-top.jpg) |

Top row: 50 N down (peak 11.56 MPa). Bottom row: 15 N sideways (peak 8.27 MPa).

## Not shown here

- The stress colours as Blender-ready `.ply` files (6 to 10 MB each) are not in the repository; the tool
  writes them with every `bars=1` and `check=` run.
- The screw pull-out numbers (runs 4 and 5) assumed a full rigid plate; they do not apply to the frame of
  runs 8 and 9. The wall plug is the weak point of the whole mounting (top screw about 106 N pull against
  650 N recommended for a fischer SX Plus 6 x 30 in concrete).
