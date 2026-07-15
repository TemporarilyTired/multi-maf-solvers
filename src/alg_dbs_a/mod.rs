//! Implementation of Rt-MAF, the depth bounded search algorithm for computing the exact maximum agreement forest of multiple rooted binary phylogenetic trees from:
//!
//! F. Shi, J. Wang, J. Chen, Q. Feng, J. Guo.
//! "Algorithms for parameterized maximum agreement forest problem on multiple trees".
//! Theoretical Computer Science, 2014.
//! DOI: https://doi.org/10.1016/j.tcs.2013.12.025

mod logging;
mod solver;
mod stack_arena_tree;
mod state;

pub use solver::solve;
pub mod variants;
