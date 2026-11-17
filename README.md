# Generative Design Engine (SolidCraft plugin)

A SolidCraft command plugin (`craft:plugin@0.1.0`, capability `commands@0.1`).
Version 0.3: an engineering **material library** that becomes SolidCraft's own physical
material cards, a per-body **material report**, and a **stress check** for flat parts. A 2D
topology optimiser exists as a benchmark; generative design itself comes later.

## Entries

| Entry | Effect | What it does |
|---|---|---|
| `gen-design-engine.apply-material` | edit | Creates the chosen library material in the design (once, with its full engineering card) and assigns it to a body (empty = every body; scripts may pass a `bodies` list). One undo step. |
| `gen-design-engine.material-report` | read | Per body: material, volume, mass, stiffness, yield/tensile strength and the allowable stress at a safety factor. |
| `gen-design-engine.stress-check` | read | Plane-stress analysis of a flat body: fixed faces, load faces and a force in N; gives peak von Mises stress and where, the safety factor against the material's strength (pass/fail), the largest displacement, and a text heat map. |
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

Allowable stress = yield strength / safety factor, or tensile strength / safety factor for
materials without a yield value (the printed plastics). **This is not a stress analysis**: the
report says how much stress a material may take, not how much the part gets.

## Material library (`materials/library.json`)

| Id | Material | Values |
|---|---|---|
| `aluminium-6061-t6` | Aluminium 6061-T6 | typical (ASM handbook-style) |
| `aluminium-7075-t6` | Aluminium 7075-T6 | typical |
| `steel-s235jr` | Steel S235JR | EN 10025-2 minimums, Eurocode E/nu/density |
| `steel-s355jr` | Steel S355JR | EN 10025-2 minimums, Eurocode E/nu/density |
| `stainless-304` | Stainless steel 304 | ASTM A240 minimum strengths, typical physical values |
| `stainless-316l` | Stainless steel 316L | ASTM A240 minimum strengths, typical physical values |
| `titanium-ti6al4v` | Titanium Ti-6Al-4V | ASTM B348 Grade 5 minimum strengths, typical physical values |
| `pla-printed` | PLA (3D printed, generic) | illustrative |
| `petg-printed` | PETG (3D printed, generic) | illustrative |
| `abs-printed` | ABS (3D printed, generic) | illustrative |
| `pa12-sls` | Nylon PA12 (SLS printed, generic) | illustrative |

Every property names its source in the card; SolidCraft shows them as "source review required".
The values were compiled from standard handbook and standard tables, **not yet checked against
the documents themselves**; the printed-plastic values are placeholders to be replaced by your
filament's datasheet or your own test prints. Printed parts are much weaker across the layers.

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
`analyze` returned invented numbers. Kept in git history and in
`winExtra/0000 Claude deprecated/SolidCraft/.generative-design-engine.bak-before-rewrite`.

MIT licence.
