//! The WebAssembly component: binds the engine to SolidCraft's `craft:plugin@0.1.0` WIT.
//! Generated binding code needs `unsafe` internally; nothing in this file adds any.
#![allow(unsafe_code)]

wit_bindgen::generate!({
    world: "command-plugin",
    path: "wit",
});

use crate::engine::{self, Failure, Level};
use craft::plugin::types::{CommandFailure, FieldFailure, PluginError};
use craft::plugin::{commands, control, log};
use serde_json::Value;

struct Component;

struct WitHost;

impl engine::Host for WitHost {
    type Error = PluginError;

    fn run(&mut self, command: &str, args: &Value) -> Result<Value, Failure<PluginError>> {
        let reply = commands::run(command, &args.to_string()).map_err(Failure::Host)?;
        serde_json::from_str(&reply).map_err(|e| Failure::Internal(format!("{command} gave a reply that is not JSON: {e}")))
    }

    fn log(&mut self, level: Level, message: &str) {
        let level = match level {
            Level::Debug => log::Level::Debug,
            Level::Info => log::Level::Info,
            Level::Warn => log::Level::Warn,
            Level::Error => log::Level::Error,
        };
        log::log(level, message);
    }

    fn progress(&mut self, fraction: f32, message: &str) {
        control::progress(fraction.clamp(0.0, 1.0), message);
    }

    fn cancelled(&mut self) -> bool {
        control::cancelled()
    }
}

impl exports::craft::plugin::entry::Guest for Component {
    fn execute(entry: String, params: String) -> Result<String, PluginError> {
        engine::execute(&mut WitHost, &entry, &params).map_err(|failure| match failure {
            Failure::Host(error) => error,
            Failure::Cancelled => PluginError::Cancelled,
            Failure::Input { field, reason } => PluginError::InvalidInput(FieldFailure { field, reason }),
            Failure::Internal(reason) => PluginError::CommandFailed(CommandFailure { command: format!("gen-design-engine.{entry}"), reason }),
        })
    }
}

export!(Component);
