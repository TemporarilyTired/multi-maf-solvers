//! Implementation of a preprocessing step, which takes a solver (function) and calls the solver
//! function on each reduced cluster, combining the result into a solution for the original
//! instance.
//!
//! Includes reduction rules from [Shi et al., 2018] (see ./../alg_dbs_lsi/mod.rs), and includes
//! cluster reduction.

mod arena_tree_ext;
mod instance_ext;
mod reduce_and_solve;

pub use reduce_and_solve::reduce_and_solve;
