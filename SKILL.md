# Generative Design Engine Plugin

An AI-powered generative design plugin for SolidCraft that implements topology optimization and stress analysis, similar to Fusion 360's generative design capabilities.

## Features

- **Topology Optimization**: Uses Solid Isoparametric (SIsoP) and Linear Static Half-Space Optimization (LS-HO) methods
- **Stress Analysis**: Real-time von Mises stress computation and safety margin evaluation
- **Material Library**: Pre-configured materials with properties (Aluminum 6061-T6, Titanium, Steel, etc.)
- **Constraint System**: Preset constraints for volume, mass, wall thickness, and safety factors
- **Lattice Generation**: Optional lattice infill structures for lightweight design

## Installation

### From ZIP Package

1. Build the plugin (see Build section)
2. Open SolidCraft
3. Go to **Preferences > Plugins**
4. Click **Install plugin** and select the ZIP file

### Manual Installation

Place the plugin folder in SolidCraft's plugins directory:
- **Linux**: `~/.local/share/solidcraft/plugins/gen-design-engine/`
- **macOS**: `~/Library/Application Support/SolidCraft/plugins/gen-design-engine/`
- **Windows**: `%LOCALAPPDATA%\SolidCraft\plugins\gen-design-engine\`

## Usage

### Optimize Entry Point

```json
{
  "command": "gen-design-engine.optimize",
  "params": {
    "constraints": {
      "material": "Aluminum 6061-T6",
      "yield_strength": 276,
      "safety_factor": 2.0
    },
    "parameters": {
      "target_volume": 5000,
      "max_mass": 50,
      "min_wall": 2
    },
    "optimization": {
      "method": "lsho",
      "iterations": 200,
      "filter_radius": 1.5
    }
  }
}
```

### Analyze Entry Point

```json
{
  "command": "gen-design-engine.analyze",
  "params": {
    "body": "Bracket",
    "stress_limits": {
      "yield_strength": 276
    }
  }
}
```

## API Reference

### `optimize` Entry (edit effect)

**Input Schema:**
- `constraints.object`: Design constraints
  - `material` (string): Material name
  - `density` (number, default: 2.7): Material density g/cm³
  - `yield_strength` (number, default: 276): Yield strength MPa
  - `safety_factor` (number, default: 2.0): Target safety factor
  - `max_stress` (number, default: 125): Maximum allowable stress MPa
- `parameters.object`: Design parameters
  - `target_volume` (number, default: 1000): Target volume mm³
  - `max_mass` (number, optional): Maximum mass grams
  - `min_bedding`/`min_wall` (number, default: 2): Minimum wall thickness mm
- `optimization.object`: Algorithm settings
  - `method` (enum: lsho|lsd|simulated_annealing|gradient, default: lsho)
  - `iterations` (integer, default: 100, range: 10-5000)
  - `filter_radius` (number, default: 1.5, range: 0.1-10): Minimum feature size
  - `convergence_threshold` (number, default: 0.01)
  - `smoothing_factor` (number, default: 0.5, range: 0-1)
- `constraints.loads.array`: Load conditions (optional)
- `output.object`: Output settings (optional)
- `tolerances` (number, default: 0.5): Convergence tolerance

**Returns:**
```json
{
  "success": true,
  "optimized_body_name": "Optimized_Bracket",
  "volume_reduction_percent": 35.5,
  "stress_distribution": {"max": 88.2, "avg": 44.1},
  "safety_margins": {"overall": 2.0, "critical_points": 1.6},
  "suggested_features": ["Add fillet at edge A", "Verify hole positions"],
  "warnings": [],
  "execution_log": ["Step 1 complete", "Step 2 complete"]
}
```

### `analyze` Entry (read effect)

**Input Schema:**
- `body` (string, required): Body name to analyze
- `stress_limits.object`: Material limits (optional)
  - `yield_strength` (number, default: 250)
  - `ultimate_strength` (number, default: 400)
  - `fatigue_strength` (number, default: 150)
- `mesh_sampling` (number, default: 0.5): Resolution factor

**Returns:**
```json
{
  "body": "Bracket",
  "max_stress": 125.5,
  "min_stress": 10.0,
  "average_stress": 44.2,
  "von_mises_stress": 78.3,
  "principal_stresses": [85.0, 65.0, 45.0],
  "safety_factor": 2.1,
  "critical_locations": [...],
  "stress_concentrations": [...]
}
```

## Building

### Requirements

- Rust 1.70+ with wasm32-unknown-unknown target
- `wasm-tools` for componentization

```bash
# Install dependencies
rustup target add wasm32-unknown-unknown
cargo install wasm-tools

# Build
./build.sh

# Or manually:
cargo build --release --target wasm32-unknown-unknown
wasm-tools component new \
    --interface craft:plugin@0.1.0 \
    --export=entry.execute=wasm32[func(func(i32,i32,i32,i32) -> i32)] \
    target/wasm32-unknown-unknown/release/generative_design_engine.wasm \
    -o plugin.wasm

# Package
zip -r gen-design-engine.zip plugin.json plugin.wasm src/
```

## SolidCraft Commands Used

The plugin internally uses these SolidCraft commands:
- `document.inspect` - Query design bodies
- `model.mesh` - Read body mesh for analysis
- `solid.cylinder`, `solid.extrude`, `solid.fillet` - Create optimized geometry
- `solid.shell` - Apply shell operations
- `parameters.change` - Update parametric constraints

## Development

### Project Structure

```
generative-design-engine/
├── Cargo.toml          # Rust package manifest
├── build.sh            # Build script
├── plugin.json         # SolidCraft plugin manifest
├── plugin.wit          # WIT interface definition
├── src/
│   └── lib.rs          # Plugin implementation
└── SKILL.md            # This documentation
```

### Testing

```bash
# Run in SolidCraft
solidcraft-cli mcp

# Or use a test script
echo '[{"command": "plugin.inspect", "params": {"path": "."}}]' > test.json
solidcraft-cli run test.json
```

## License

MIT OR Apache-2.0