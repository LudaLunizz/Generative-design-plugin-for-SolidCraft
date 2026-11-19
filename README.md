# Generative Design Engine (SolidCraft plugin)

A [SolidCraft](https://github.com/tangency-labs/solidcraft) command plugin (`craft:plugin@0.1.0`,
capability `commands@0.1`). Version 0.4: an engineering **material library** that becomes
SolidCraft's own physical material cards, a per-body **material report**, a **stress check**
for flat parts, and a **plate optimiser** that designs a lightweight plate and draws it as a sketch
and an extrude. Next to the plugin there are three command-line tools (3D stress check,
printed-part stress check, beam optimiser; see below).

Independent project: not made, endorsed or supported by Tangency Labs or the SolidCraft
maintainers. "SolidCraft" is their name. Version 0.4.1 is tested with **SolidCraft 0.44.1**
(`tests/cli-test.sh`); 0.3.x also ran on 0.43.0.

> **A design aid, not a certified analysis.** The numbers can be wrong. A qualified person must
> check them before anyone relies on them for a design, manufacturing or safety decision. The
> material values are listed below with their sources: they were checked against public data
> sheets in 0.4.0, but not against the (paywalled) standards themselves.

## Entries

| Entry | Effect | What it does |
|---|---|---|
| `gen-design-engine.apply-material` | edit | Creates the chosen library material in the design (once, with its full engineering card) and assigns it to a body (empty = every body; scripts may pass a `bodies` list). One undo step. |
| `gen-design-engine.material-report` | read | Per body: material, volume, mass, stiffness, yield/tensile strength and the allowable stress at a safety factor. |
| `gen-design-engine.stress-check` | read | Plane-stress analysis of a flat body: fixed faces, load faces and a force in N; gives peak von Mises stress and where, the safety factor against the material's strength (pass/fail), the largest displacement, and a text heat map. |
| `gen-design-engine.optimise-plate` | edit | Topology-optimises a flat plate (cantilever or MBB half beam) for stiffness at a chosen material share, then draws the result as a sketch of closed outlines and extrudes it into a new body with the chosen library material. One undo step. Reports mass, outlines, holes and a rough load-point deflection. |
| `gen-design-engine.materials` | read | Lists the library with its values. Runs no app commands. |
| `gen-design-engine.solver-benchmark` | read | Runs the 2D SIMP optimiser on a test case and reports the work done (development tool). |

In SolidCraft: **Model › Plugins** lists them (and, where the version supports placement,
Manage and Inspect ribbon groups). Every entry's form uses only fields SolidCraft's dialog can
draw (a test checks this). From MCP or a script: `plugin.run` with the entry name.

### Stress check

- Faces: a face name as `model.faces` lists it (`F1:-x`, `F2:bore:side:0`), its number, or a
  side of the part (`-x`, `+y`, ...); several with commas. Supports can be holes (cylinder faces
  through the thickness) or flat side faces.
- The plate's plane is across its thinnest direction unless chosen. The force's part through
  the plate is ignored, with a warning: bending out of the plane is not covered.
- Checked against beam theory: a 100 x 10 x 5 mm aluminium cantilever under 100 N gives a tip
  deflection of 1.170 mm (theory 1.17) and a peak of 121 MPa (theory 120) in SolidCraft.
- Limits: linear elastic, small deflections, rigid supports, the load spread evenly; flat parts
  of one thickness (it warns when the outline changes through the thickness or the part is
  thick); peaks at sharp inside corners grow with finer grids (it warns; round the corner, or
  judge by the 99th percentile it also reports). A design check, not a certified analysis.

### Plate optimiser

- Input: width, height, thickness, material, load case, force, material share (0.1 to 0.9), grid
  cells across the width (20 to 120; the whole grid is capped at 8000 cells so a run finishes
  within SolidCraft's 30 s limit), design iterations and the smallest member (in cells).
- The shape comes from 2D SIMP optimisation (density filter, optimality criteria, multigrid
  solver). Cells above 0.5 density are material; only material joined to the fixed edge is kept;
  diagonal-only contacts are closed; the staircase boundary is smoothed (real corners are kept) and
  drawn as closed sketch lines, with holes as inner loops, and extruded.
- It is checked in SolidCraft: the extruded solid has the same volume as the optimiser's own area,
  and one undo removes sketch, body and material assignment.
- **Stiffness only.** Strength, stress concentrations, buckling and manufacturability are not
  checked. Run `stress-check` on the new body before relying on it. The deflection figure is a
  rough estimate (a force on one grid node).
- The force changes the deflection estimate, not the shape: one linear material and one load.

Allowable stress = yield strength / safety factor, or tensile strength / safety factor for
materials without a yield value (the printed plastics). **This is not a stress analysis**: the
report says how much stress a material may take, not how much the part gets.

## Material library (`materials/library.json`)

| Id | Material | Values |
|---|---|---|
| `aluminium-6061-t6` | Aluminium 6061-T6 | specification minimums for strength (yield 240, tensile 260 MPa), ASM typical for the rest |
| `aluminium-7075-t6` | Aluminium 7075-T6 | specification minimums for strength (yield 460, tensile 540 MPa), ASM typical for the rest |
| `steel-s235jr` | Steel S235JR | EN 10025-2 minimums, Eurocode E/nu/density |
| `steel-s355jr` | Steel S355JR | EN 10025-2 minimums, Eurocode E/nu/density |
| `stainless-304` | Stainless steel 304 | ASTM A240 minimum strengths, typical physical values |
| `stainless-316l` | Stainless steel 316L | ASTM A240 minimum strengths, typical physical values |
| `titanium-ti6al4v` | Titanium Ti-6Al-4V | ASTM B348 Grade 5 minimum strengths, typical physical values |
| `pla-printed` | PLA (3D printed, generic) | Prusament PLA data sheet (printed specimens) |
| `petg-printed` | PETG (3D printed, generic) | Prusament PETG data sheet (printed specimens) |
| `abs-printed` | ABS (3D printed, generic) | Bambu Lab ABS data sheet (annealed specimens) |
| `pa12-sls` | Nylon PA12 (SLS printed, generic) | EOS PA 2200 data sheet (laser sintered) |

Every property names its source in the card; SolidCraft shows them as "source review required".
**Checked on 19 Nov 2026 (0.4.0):** the metal values against public data sheets (ASM data via
reseller copies, mill data sheets, Eurocode summaries; the standards themselves are paywalled), and
the four plastics replaced by manufacturer data sheets of printed specimens. One correction came
out of it (7075 thermal expansion 23.6, not 23.4). Things to know:

- 6061-T6 and 7075-T6 use guaranteed **specification minimums** for yield and tensile strength
  (since 0.4.1; 0.4.0 had the ASM typical values 276 / 310 and 503 / 572 MPa, which the Aluminum
  Association marks not for design). 6061-T6: 240 / 260 MPa, the extrusion and bar minimum (sheet
  and plate guarantee 240 / 290). 7075-T6: 460 / 540 MPa, plate 6.3-12.5 mm unclad (ASTM
  B209M-14 table 3, Aluminum Association table 11.1, EN 485-2 / 755-2). Thicker plate, clad sheet
  and other forms differ; check the supplier's certificate. No shear strength is listed, because no
  specification guarantees one.
- S355JR tensile is the EN 10025-2 product minimum (470 MPa); Eurocode 3 designs with 490 MPa.
  Yield drops for sections thicker than 16 mm.
- 304 and 316L follow the ASTM/ASM side (density 8000, E 193 GPa); EN 1.4301/1.4404 data sheets
  use slightly different physical values. For 316L, Poisson ratio, conductivity and expansion are
  widely quoted typical values that were not on the ASM sheet.
- Printed plastics are **one brand each** and the Poisson ratios are generic estimates. The ABS
  sheet's specimens were annealed (stronger than ordinary parts), so its in-plane modulus
  (3.8 GPa) is not used; the library has the across-layers 2.96 GPa. Across the layers (Z) the
  datasheets give lower strengths (PLA 38, PETG 30, ABS 32, PA12 42 MPa); the library holds the
  in-plane values, so a part loaded across its layers is weaker than the card says.

To add or change a material: edit `materials/library.json`, run `./build.sh` (it regenerates
`plugin.json`, runs every check and packages the zip), then `tests/cli-test.sh`.

A design material with the same name as a library material but other values is never overwritten:
`apply-material` refuses and says so. The report leaves out strengths when a body's material has
a library name but a different density.

## Build, test, install

```bash
rustup target add wasm32-unknown-unknown
cargo install wasm-tools --locked
./build.sh              # tests, clippy -D warnings, wasm build, component, dist/*.zip + .sha256
tests/cli-test.sh       # end to end in the installed SolidCraft (headless solidcraft-cli)
```

Install: **Preferences › Plugins › Install plugin…** and pick `dist/gen-design-engine-<version>.zip`,
or `plugin.install {"path": ".../gen-design-engine-<version>.zip", "replace": true}`.
The plugin needs no command groups, folders, clock or network.

## Command-line tools (outside SolidCraft)

These run natively (no plugin limits on time or memory) and are built with
`cargo build --release --bins`. They are experimental.

| Tool | What it does |
|---|---|
| `stress3d` | 3D voxel stress check of an STL part: majority-vote voxelisation, preconditioned CG solve, von Mises heat maps (PNG) and a report (`report.json`) with the safety factor against a library material. |
| `stress3d_printed` | The same for a **3D-printed** part: solid walls and top/bottom skins around a lighter infill core, with layers along z. Infill stiffness is `infill^exp` (exp 1 = upper bound, exp 2 = lower bound). |
| `beam_opt` | Lightest plate-truss for a long cantilever beam: SIMP topology optimisation with a min-member filter, then strut and sideways (lateral-torsional) buckling checks. |

Usage is in the header comment of each file in `src/bin/`. Example:

```bash
cargo run --release --bin stress3d -- part.stl --material aluminium-6061-t6 \
  --fix zmin --pin z,50,30,4 --case pull:0,0,-200 --voxel 1.5 --out out/
```

The optional **GPU solver** (wgpu, f64 multigrid-preconditioned CG) speeds up `beam_opt`
(`gpu=1`) and needs a GPU with 64-bit float support:
`cargo build --release --features gpu --bins`. Without the feature, or without a suitable GPU,
the CPU solver is used. The plugin's WebAssembly build never uses it.

A design from `beam_opt` is a starting point: the plate-truss weight is limited by strut and
lateral-torsional buckling, not by stiffness or strength, and none of this replaces a proper
structural calculation.

## Hardening notes

- Inputs are hostile: unknown fields, non-objects, bad numbers, empty/duplicate/oversized body
  names are refused with `invalid-input` before any edit. No unwrap/expect/panic/indexing in
  non-test code (clippy-denied); `unsafe` only in the generated WIT bindings.
- Host errors pass through unchanged, so SolidCraft rolls the transaction back.
- Cancellation is checked before each edit.
- Limits: 2e9 fuel, 64 MiB, 15 s: inside SolidCraft's default ceilings, so no extra approval.
- The component imports only `craft:plugin` interfaces (no WASI); `build.sh` refuses anything else.

## SolidCraft 0.42 notes

- `material.list` and `material.inspect` only read, but 0.42 does not mark them read-only, so a
  read entry may not call them. The report and the stress check therefore take strengths from the
  library (matched by name, checked by density).
- The plugin dialog draws numbers, whole numbers, checkboxes, dropdowns and text only; a list
  field (such as `bodies`) switches Run off. Hence the `body` text field.
- 0.42's `plugin.install` takes only `path`, `replace` and `dry_run`.

## History

0.1.0 (not working): the manifest granted a non-existent command and asked for more than the
host allows; `plugin.wasm` was SolidCraft's sample box component, not this code; `optimize` and
`analyze` returned invented numbers. It stays in the git history.

## License

Dual-licensed, like SolidCraft itself: **MIT OR Apache-2.0**, at your option. See
[`LICENSE-MIT`](LICENSE-MIT) and [`LICENSE-APACHE`](LICENSE-APACHE). Copyright (c) 2026 LudaLunizz.
The software is provided "as is", without warranty of any kind; see both licences.
See [SECURITY.md](SECURITY.md) to report a vulnerability and [CONTRIBUTING.md](CONTRIBUTING.md) to contribute.
