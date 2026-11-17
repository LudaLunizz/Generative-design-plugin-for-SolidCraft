# Generative Design Engine (SolidCraft plugin)

A SolidCraft command plugin (`craft:plugin@0.1.0`, capability `commands@0.1`).
Version 0.2 is the foundation: an engineering **material library** that becomes SolidCraft's
own physical material cards, and a per-body **material report**. Generative design (topology
optimisation) and stress analysis come later, once a real solver is behind them.

## Entries

| Entry | Effect | What it does |
|---|---|---|
| `gen-design-engine.apply-material` | edit | Creates the chosen library material in the design (once, with its full engineering card) and assigns it to the bodies (empty list = every body). One undo step. |
| `gen-design-engine.material-report` | read | Per body: material, volume, mass, stiffness, yield/tensile strength and the allowable stress at a safety factor. |
| `gen-design-engine.materials` | read | Lists the library with its values. Runs no app commands. |

In SolidCraft: **Model › Plugins** lists all three (and, where the version supports placement,
Manage and Inspect ribbon groups). From MCP or a script: `plugin.run` with the entry name.

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
  read entry may not call them. The report therefore takes strengths from the library.
- 0.42's `plugin.install` takes only `path`, `replace` and `dry_run`.

## History

0.1.0 (not working): the manifest granted a non-existent command and asked for more than the
host allows; `plugin.wasm` was SolidCraft's sample box component, not this code; `optimize` and
`analyze` returned invented numbers. Kept in git history and in
`winExtra/0000 Claude deprecated/SolidCraft/.generative-design-engine.bak-before-rewrite`.

MIT licence.
