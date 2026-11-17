# Generative Design Engine Plugin

An AI-powered generative design plugin for SolidCraft that implements topology optimization and stress analysis, similar to Fusion 360's capabilities.

## 📁 Package Contents

| File | Description |
|------|-------------|
| `plugin.json` | SolidCraft plugin manifest |
| `plugin.wat` | WebAssembly component source (WAT text format) |
| `plugin.wasm` | Compiled WebAssembly component |
| `src/` | TypeScript/JavaScript reference implementation |
| `src/lib.rs` | Rust core logic |
| `BUILD.md` | Build instructions |
| `SKILL.md` | Full plugin documentation |

## 🚀 Quick Start

### Installation

1. Open SolidCraft
2. Go to **Preferences > Plugins**
3. Click **Install plugin** and select this ZIP file

**OR** manually place the folder:
```bash
# Linux
mkdir -p ~/.local/share/solidcraft/plugins/gen-design-engine
cp -r *.json *.wasm src/ ~/.local/share/solidcraft/plugins/gen-design-engine/

# macOS
mkdir -p ~/Library/Application Support/SolidCraft/plugins/gen-design-engine
cp -r *.json *.wasm src/ ~/Library/Application Support/SolidCraft/plugins/gen-design-engine/

# Windows
# %LOCALAPPDATA%\SolidCraft\plugins\gen-design-engine\
```

## 🔧 Plugin Entries

### `optimize` (edit effect)
Runs topology optimization with constraints and parameters.
```json
{
  "command": "gen-design-engine.optimize",
  "params": {
    "constraints": {"material": "Aluminum 6061-T6", "yield_strength": 276, "safety_factor": 2.0},
    "parameters": {"target_volume": 5000, "min_wall": 2},
    "optimization": {"method": "lsho", "iterations": 200}
  }
}
```

### `analyze` (read effect)
Analyzes stress distribution on a body.
```json
{
  "command": "gen-design-engine.analyze",
  "params": {"body": "Bracket", "stress_limits": {"yield_strength": 250}}
}
```

## 📚 Documentation

See `SKILL.md` for complete API reference, material properties, and Fusion 360 compatibility notes.

## 🏗️ Architecture

The plugin follows the SolidCraft Plugin RFC 0002:
- WebAssembly sandboxed component
- Invokes native SolidCraft commands through the provided `run()` interface
- Produces validated native features that save and rebuild without the plugin
- One undo step per successful call

## 🔗 Resources

- [SolidCraft Plugins Documentation](docs/plugins.md)
- [Plugin API WIT](plugin.wit)
- [Fusion 360 Generative Design Reference](https://help.autodesk.com/view/fusion360/ENU/?guid=GUID-...)
