# Building the Plugin

## Option A: Use pre-built (when available)
Download `plugin.wasm` from releases.

## Option B: Build locally

```bash
# Install tools
cargo install wit-bindgen-cli
cargo install wasm-tools

# Generate bindings
wit-bindgen generate --world command-plugin --target wasm32 \
  plugin.wit --out-dir src/generated

# Update src/lib.rs to use generated bindings
# Then build
cargo build --release --target wasm32-unknown-unknown
wasm-tools component new \
  --interface craft:plugin@0.1.0 \
  target/wasm32-unknown-unknown/release/generative_design_engine.wasm \
  -o plugin.wasm
```

## Current Status
- ✅ Rust implementation complete (src/lib.rs)
- ✅ Plugin manifest valid (plugin.json)
- ✅ WIT interface copied (plugin.wit)
- ⏳ WebAssembly component requires wit-bindgen

The Rust core library compiles successfully. The plugin.wasm component will be generated when wit-bindgen bindings are available.
