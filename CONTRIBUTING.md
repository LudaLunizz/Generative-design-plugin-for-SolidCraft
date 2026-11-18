# Contributing

Issues and pull requests are welcome. The [README](README.md) explains what the
plugin does and how to install it; this file is the setup and the exact commands.

## Setup

| Path | What |
|---|---|
| `src/` | the plugin (WebAssembly component), the solvers and the command-line tools (`src/bin/`) |
| `materials/library.json` | the material library; `plugin.json` is generated from it |
| `tests/` | `manifest.rs` (the manifest and its dialog fields) and `cli-test.sh` (end to end in an installed SolidCraft) |
| `tools/gen_manifest.py` | writes `plugin.json` from the library |

You need Rust (stable, 1.85 or newer), the `wasm32-unknown-unknown` target and
`wasm-tools`:

```sh
rustup target add wasm32-unknown-unknown
cargo install wasm-tools --locked
```

## Before you open a pull request

These are what CI runs (`.github/workflows/ci.yml`):

```sh
./build.sh                         # tests, clippy -D warnings, wasm build, component, dist/*.zip
cargo check --all-targets --features gpu    # the optional GPU solver still builds
```

`tests/cli-test.sh` runs the packaged plugin in an installed SolidCraft
(`solidcraft-cli`); CI cannot, so say in the pull request whether you ran it and
with which SolidCraft version.

## Conventions

- **No `unwrap`, `expect`, `panic` or indexing in non-test code** (clippy denies
  them); inputs are treated as hostile.
- **The version is in `Cargo.toml` and `tools/gen_manifest.py`**; change both.
- **A material value names its source** in `materials/library.json`. Do not add
  values from memory; cite the standard or datasheet.
- **GitHub Actions are pinned to a full commit SHA** with the version in a
  comment; Dependabot proposes updates.

## Issues and pull requests

Open an issue first for anything beyond a typo, saying what you want and why.
Work on a branch named `<issue number>-short-name`, and start the pull request
description with `Closes #<number>`. Wait for CI to be green; the maintainer
merges with a rebase, so the history stays linear. Security problems go to
**Security → Report a vulnerability**, not to a public issue (see
[SECURITY.md](SECURITY.md)).

Developed with help from Claude Code; every change is reviewed, tested and
released by the maintainer.
