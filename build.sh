#!/usr/bin/env bash
# Build the plugin component and package it as dist/gen-design-engine-<version>.zip.
# Needs: rustup target wasm32-unknown-unknown, wasm-tools (cargo install wasm-tools --locked).
set -euo pipefail
cd "$(dirname "$0")"
export PATH="$HOME/.cargo/bin:$PATH"

version=$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -n 1)
python3 -I tools/gen_manifest.py
cargo test --quiet
cargo clippy --all-targets --quiet -- -D warnings
cargo build --release --target wasm32-unknown-unknown --quiet

core=target/wasm32-unknown-unknown/release/generative_design_engine.wasm
stage=target/package/gen-design-engine
rm -rf target/package && mkdir -p "$stage" dist
wasm-tools component new "$core" -o "$stage/plugin.wasm"
wasm-tools validate --features component-model "$stage/plugin.wasm"
# SolidCraft links no WASI: refuse a component that imports anything outside craft:plugin.
if wasm-tools component wit "$stage/plugin.wasm" | grep -E '^\s*import ' | grep -v 'craft:plugin/'; then
  echo "component imports something outside craft:plugin" >&2
  exit 1
fi
cp plugin.json README.md LICENSE "$stage/"
zip_path="dist/gen-design-engine-$version.zip"
rm -f "$zip_path"
(cd target/package && zip -qr -X "../../$zip_path" gen-design-engine)
sha256sum "$zip_path" | tee "$zip_path.sha256"
