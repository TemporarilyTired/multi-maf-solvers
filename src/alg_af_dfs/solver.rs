use super::state::init_bu_state;

use crate::{
    maf_instance::{arena_tree::ArenaTree, instance::Instance},
    preprocessing::reduce_and_solve,
};

// Returns a MAF for the instance
// Perform global preprocessing, calling solve_fn on every produced cluster;
pub fn solve(instance: Instance) -> ArenaTree {
    reduce_and_solve(instance, &solve_fn)
}

// Returns a MAF for the instance;
// Finds this MAF using the AF-DFS algorith
fn solve_fn(instance: Instance) -> ArenaTree {
    let mut bu_state = init_bu_state(instance);
    bu_state.solve()
}
