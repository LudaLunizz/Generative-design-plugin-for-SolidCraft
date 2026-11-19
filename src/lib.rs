//! Generative Design Engine: a SolidCraft command plugin (`craft:plugin@0.1.0`, `commands@0.1`).
//!
//! Version 0.2 ships the engineering material library and a per-body material report. The
//! generative design (topology optimization) and stress analysis entries of 0.1 returned made-up
//! numbers and are withdrawn until a real solver backs them.

pub mod beam;
pub mod engine;
#[cfg(feature = "gpu")]
pub mod gpu;
pub mod library;
pub mod mg2d;
pub mod mg3d;
pub mod outline;
pub mod stress2d;
pub mod stress3d;
pub mod topo2d;

#[cfg(target_arch = "wasm32")]
mod component;
