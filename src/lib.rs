//! A Rust translation of CD-HIT 4.8.1.
//!
//! Module layout mirrors the C++ sources one-to-one; see `NAMING.md` at the repository
//! root for the name mapping and the small number of deliberate structural changes.

pub mod cdhit_common;
pub mod cdhit_utility;

// One module per C++ front-end, each holding that program's `main` as
// `pub fn main(argv: &[String])`. They are dispatched by `src/bin/cdhit.rs`, and are public
// so an embedding crate can call them directly instead of spawning a process.
pub mod cdhit;
pub mod cdhit_2d;
pub mod cdhit_454;
pub mod cdhit_div;
pub mod cdhit_est;
pub mod cdhit_est_2d;

/// A program's entry point: the translation of one `cdhit-*.c++` `main`, taking the
/// argument vector the standalone executable would have received.
pub type ProgramMain = fn(&[String]);

/// The six programs, by their upstream names, with the short aliases the single binary
/// also accepts.
pub const PROGRAMS: &[(&str, &str, ProgramMain)] = &[
    ("cd-hit", "prot", cdhit::main),
    ("cd-hit-est", "est", cdhit_est::main),
    ("cd-hit-2d", "2d", cdhit_2d::main),
    ("cd-hit-est-2d", "est-2d", cdhit_est_2d::main),
    ("cd-hit-454", "454", cdhit_454::main),
    ("cd-hit-div", "div", cdhit_div::main),
];

/// Resolve a program name or alias to its entry point.
pub fn program(name: &str) -> Option<ProgramMain> {
    PROGRAMS
        .iter()
        .find(|(full, alias, _)| *full == name || *alias == name)
        .map(|(_, _, f)| *f)
}
