//! Implementation for BINARY input trees of `Alg-MAF`, the depth-bounded search algorithm for computing the exact maximum agreement forest of multiple rooted (originally multifurcating) phylogenetic trees from:
//!
//! F. Shi, J. Chen, Q. Feng, J. Wang.
//! "A parameterized algorithm for the Maximum Agreement Forest problem on multiple rooted multifurcating trees".
//! Journal of Computer and System Sciences, 2018.
//! DOI: https://doi.org/10.1016/j.tcs.2013.12.025

pub mod arena_tree_lsi_ext;
mod branch_result;
pub mod instance_lsi_ext;
mod logging;
mod lsi_instance;
mod lsi_solver;
mod r1_workspace;

pub use lsi_solver::solve;
