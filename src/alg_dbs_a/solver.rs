use std::collections::HashSet;

use crate::maf_instance::arena_tree::ArenaTree;
use crate::maf_instance::arena_vertex::Idx;
use crate::maf_instance::instance::Instance;

use super::stack_arena_tree::ArenaTreeMaximalExt;
use super::state::{
    Branch::{self, *},
    Forest::*,
    Frame, State,
};
use super::variants::Variant;

// states are uniquely identified by the sorted set of root indices of f1
type StateHash = Vec<Idx>;

// import log macro; and only import other logging things if enabled
use super::logging::log;
#[cfg(feature = "logging")]
use super::logging::{LOGS, Logs};

// Use at most 7GB of memory for the caching to keep the full memory usage under 8GB
const MAX_MEM_GB: usize = 7;
const MAX_MEM_B: usize = MAX_MEM_GB * 1024 * 1024 * 1024;

pub fn solve<V: Variant>(instance: Instance) -> ArenaTree {
    #[cfg(feature = "logging")]
    println!("# running maximal af solver (variant: {})", V::NAME);

    if V::ENABLE_GLOBAL_PREPROCESSING {
        use crate::preprocessing::reduce_and_solve;
        reduce_and_solve(instance, &solve_fn::<V>)
    } else {
        solve_fn::<V>(instance)
    }
}

/// Computes a maximum (i.e., minimum-order) agreement forest for the given instance.
///
/// This method uses the depth bounded search algorithm in [Shi et al., 2014] as the subroutine until an optimal solution is found. Depending on the Variant, this method either performs binary search on the order of the optimal solution or iteratively checks k=1 to k=n until a solution is found.
/// The constant (compile time) booleans in the Variant type (chosen using command line arguments) also determine the use of other
/// features:
///
/// - modified branching rules when only two forests remain in the instance,
/// - lower bound calculation to prune execution branches,
/// - caching of infeasible states.
///
/// # Parameters
///
/// # `instance` - Input instance consisting of the forests to compute the MAF for
///
/// # Returns
///
/// A minimum-order agreement forests represented as an `ArenaTree`.
fn solve_fn<V: Variant>(mut instance: Instance) -> ArenaTree {
    if instance.is_completed() {
        return instance.extract_af();
    }

    // instantiate cache hashset if caching is enabled,
    // otherwise compiler should recognize visited is always None and compile out all caching
    let mut visited = if V::ENABLE_CACHING {
        // hashset storing the hashes of visited (pruned) states
        // assume hashset capacity <= 3 * hashset length; assume vector capacity <= 2 * vector length
        let hash_size_b = size_of::<StateHash>() * 3 + size_of::<Idx>() * instance.num_leaves * 2;
        let max_visited_states = MAX_MEM_B / hash_size_b;
        Some(HashSet::<StateHash>::with_capacity(max_visited_states))
    } else {
        None
    };
    let mut last_ran_k = if V::ENABLE_CACHING {
        Some(instance.num_leaves)
    } else {
        None
    };

    // find trivial upper bound
    let mut sol = dbs::<V>(instance.clone(), instance.num_leaves, visited.as_mut())
        .expect("there must be a valid AF of size #leaves");
    let mut ub = sol.ord();

    // find simple lower bound (or use lb=1)
    let mut lb = if V::ENABLE_LB_CALCULATION {
        let mut temp_state = State::new_from_instance(instance.clone());
        apply_reductions(&mut temp_state);
        calculate_simple_lb(&mut temp_state, ub)
    } else {
        1
    };

    // Depending of Variant, either:
    // - perform binary search on the order of the MAF, or
    // - just iterate k=1..n
    while lb < ub {
        let k = if V::ENABLE_BS_ON_SOLUTION {
            (lb + ub) / 2
        } else {
            lb
        };

        // reset all counters
        log!(|logs: &mut Logs| logs.reset_logs_for_new_k());

        // if we inspect a higher k, clear the visited states because some may have become viable
        if let Some(some_visited) = &mut visited {
            if last_ran_k.is_some_and(|last| last < k) {
                some_visited.clear();
            }
            last_ran_k = Some(k);
        }
        let result = dbs::<V>(instance.clone(), k, visited.as_mut());

        // print the call count logs
        log!(|logs: &Logs| logs.print_logs_after_k(k));
        log!(|logs: &mut Logs| logs.total_loops += logs.loops);
        log!(|logs: &mut Logs| logs.reset_logs_for_new_k());

        if let Some(maf) = result {
            ub = maf.ord();
            sol = maf;
        } else {
            lb = k + 1;
        }
    }
    sol
}

/// Finds an agreement forest of the instance of order at most `k` if it exists.
///
/// This method implements the depth bounded search algorithm, Rt-MAF, in [Shi et al., 2014].
///
/// The method is implemented without recursion, instead using a stack containing the (performed) branches and another stack containing information necessary to undo all performed operation to reset the state of the instance to before the branch was applied.
///
/// Exact functionality depends on Variant (see solve()):
///
/// - optionally uses given set of visited (infeasible) states to speed up the search and expands this set with newly found infeasible states,
/// - optionally calculates a lower bound each iteration to attempt to prune the branch.
///
/// # Parameters
///
/// # `instance` - Input instance.
/// # `k` - Maximum allowed order of the resulting agreement forest.
/// # `visited` - Optional cache of previously proven infeasible states; expected to be infeasible
/// for the current paremeter `k`.
///
/// # Returns
///
/// # `Some(solution)` if an agreement forest of order at most `k` exists.
/// # `None` if no such forest exists.
fn dbs<V: Variant>(
    instance: Instance,
    k: usize,
    mut visited: Option<&mut HashSet<StateHash>>,
) -> Option<ArenaTree> {
    // assume hashset capacity <= 3 * hashset length; assume vector capacity <= 2 * vector length
    let hash_size_b = size_of::<StateHash>() * 3 + size_of::<Idx>() * k * 2;
    let max_visited_states = MAX_MEM_B / hash_size_b;

    let mut state = State::new_from_instance(instance);
    log!(|logs: &mut Logs| logs.max_cache_size = max_visited_states as u64);

    // apply reductions before first branch
    apply_reductions(&mut state);

    // check if we already exceeded the target order
    if state.ord() > k {
        return None;
    }
    // then check if any branches are actually needed
    if state.is_solved() {
        return Some(state.extract_solution());
    }

    // otherwise: calculate initial branching options
    let mut stack: Vec<Frame> = vec![Frame {
        checkpoint: 0,
        branches: compute_branches::<V>(&state, k),
        next_branch: 0,
    }];

    while let Some(frame) = stack.last_mut() {
        log!(|logs: &mut Logs| logs.loops += 1);

        if frame.next_branch >= frame.branches.len() {
            // all branches were unsuccesful: the current state is infeasible
            // add its hash to the visited states
            if let Some(some_visited) = &mut visited {
                if some_visited.len() >= max_visited_states {
                    some_visited.clear();
                    log!(|logs: &mut Logs| logs.n_cache_flushes += 1);
                }
                some_visited.insert(state.get_hash());
            }

            state.rollback(frame.checkpoint);
            stack.pop();
            continue;
        }

        let branch = frame.branches[frame.next_branch].clone();
        frame.next_branch += 1;

        let checkpoint = state.undos.len();
        apply_branch(&mut state, &branch);

        apply_reductions(&mut state);

        if let Some(some_visited) = &mut visited {
            // if caching is enabled: check if state was visited before
            if some_visited.contains(&state.get_hash()) {
                state.rollback(checkpoint);
                continue;
            }
        }

        // calculate O(ntk) lower bound for this branch
        let ord = state.ord();
        let lower_bound: usize = if V::ENABLE_LB_CALCULATION && ord + 4 <= k {
            calculate_simple_lb(&mut state, k)
        } else {
            ord
        };

        // prune branch
        if lower_bound > k {
            state.rollback(checkpoint);
            continue;
        }

        if state.is_solved() {
            return Some(state.extract_solution());
        }

        let branches = compute_branches::<V>(&state, k);

        stack.push(Frame {
            checkpoint,
            branches,
            next_branch: 0,
        });
    }
    None
}

/// Computes branching options for the current state.
///
/// Branching options are as described in line 6 and line 8 in Rt-MAF in [Shi et al., 2014]; all
/// other operations that are described as branches in Rt-MAF are instead (assumed to be) performed exhaustively as
/// reductions before this function call.
///
/// Exact functionality depends on Variant (see solve()):
///
/// - modified branching rules when only two forest remain in the instance,
///
/// # Parameters
///
/// # `state` - Current state of the search, containing the input forests and undo stack.
/// # `k` - Maximum allowed order of teh resulting agreement forest.
///
/// # Returns
///
/// A list of candidate branches to explore from the current state (possibly empty if state is
/// infeasible).
pub fn compute_branches<V: Variant>(state: &State, k: usize) -> Vec<Branch> {
    let mut branches: Vec<Branch> = vec![];
    // find a sibling pair
    let Some(((label1, leaf1_f2), (label2, leaf2_f2))) = state.f2().find_a_sibling_pair() else {
        // current forest is completed
        return vec![];
    };
    debug_assert_ne!(label1, label2);

    // we need at least one more cut at this point
    let f1 = state.f1();
    if f1.ord() >= k {
        return vec![];
    }

    // check if the labels are in the same components in f1
    // if so: add branch for cutting all pendant subtrees
    let leaf1_f1 = f1.locate_label(label1);
    let leaf2_f1 = f1.locate_label(label2);
    // we need some lca and to connect them we need path_len - 3 cuts
    // because every node on the path except for the 2 leaves and the lca will have a child cut
    if let Some((lca_f1, path_len)) = f1.get_lca_and_path_len(leaf1_f1, leaf2_f1)
        && f1.ord() + path_len - 3 <= k
    {
        branches.push(Branch::CutPathBetween {
            leaf_a: leaf1_f1,
            leaf_b: leaf2_f1,
            lca: lca_f1,
        });

        // if we don't have any other trees after this
        // we can find ANY agreement forest instead of
        // iterating all maximal ones
        if V::ENABLE_BETTER_BRANCHING_ON_LAST_TREE && path_len == 4 && state.is_on_last_tree() {
            return branches;
        }
    }
    // add branches for cutting a label in both trees
    branches.push(CutLeaf {
        leaf_f1: leaf1_f1,
        leaf_f2: leaf1_f2,
    });
    branches.push(CutLeaf {
        leaf_f1: leaf2_f1,
        leaf_f2: leaf2_f2,
    });
    branches
}

/// Applies a branching decision to the current state.
///
/// Supported branch types are:
///
/// - `CutLeaf`: remove a leaf from both active forests.
/// - `CutPathBetween`: cut all pendant subtrees along the path connecting two leaves in F1.
///
/// # Parameters
///
/// # `state` - Current state of the search, containing the input forests and undo stack.
/// # `branch` - Branch to apply.
///
/// # Notes
///
/// This function does not perform reductions automatically. Call `apply_reductions()` afterwards
/// to restore expected invariants.
pub fn apply_branch(state: &mut State, branch: &Branch) {
    match branch {
        CutLeaf { leaf_f1, leaf_f2 } => {
            log!(|logs: &mut Logs| logs.cut_leaf_branches += 1);
            state.instance.forests[state.f1_idx].cut_branch_and_contract_with_undo(
                *leaf_f1,
                F1,
                &mut state.undos,
            );
            state.instance.forests[state.f2_idx].cut_branch_and_contract_with_undo(
                *leaf_f2,
                F2,
                &mut state.undos,
            );
        }
        CutPathBetween {
            leaf_a,
            leaf_b,
            lca,
        } => {
            log!(|logs: &mut Logs| logs.cut_path_branches += 1);
            let f1 = &mut state.instance.forests[state.f1_idx];
            for leaf in [*leaf_a, *leaf_b] {
                while let Some(parent) = f1.get(leaf).parent
                    && parent != *lca
                {
                    let [left, right] = f1.get(parent).children_unchecked();
                    if left != leaf {
                        f1.cut_branch_and_contract_with_undo(left, F1, &mut state.undos);
                    } else {
                        debug_assert!(right != leaf);
                        f1.cut_branch_and_contract_with_undo(right, F1, &mut state.undos);
                    }
                }
            }
        }
    }
}

/// Exhaustively applies all reduction rules to the current state, moving on to new forests
/// whenever possible.
///
///
/// These reduction are in fact exhaustive variants of the routines described in line 3, line 4 and
/// line 7 in Rt-MAF in [Shi et al., 2014]
///
/// The implementation relies on the fact that during the Rt-MAF algorithm, there are no single
/// vertex trees in F2 that are not SVT in F1 as well.
///
/// # Parameters
///
/// # `state` - Current state of the search, containing the input forests and undo stack.
///
/// # Invariants
///
/// After completion:
///
/// - F1 and F2 contain no common cherries,
/// - every label is a single vertex tree in F1 if and only if it is one in F2,
/// - the active forest is either not completed or it is the final forest.
pub fn apply_reductions(state: &mut State) {
    // exhaustively apply reduction after cutting an edge to maintain invariants
    fn reductions(state: &mut State) {
        log!(|logs: &mut Logs| logs.full_reductions += 1);
        state.sync_single_vertex_trees();
        state.exhaustively_merge_agreeing_subtrees();
    }

    reductions(state);

    // load next tree if current one is completed and keep applying reductions
    while state.current_forest_is_completed() && !state.is_on_last_tree() {
        state.load_next_tree();

        reductions(state);
    }
    #[cfg(feature = "assert_validity")]
    state.assert_fully_reduced();
}

/// Calculates a simple lower bound on the order of any agreement forest reachable from the current
/// state.
///
/// The implementation is a modified version of the 3-approximation algorithm `Apx-MAF` from:
///
/// J. Chen, F. Shi, J. Wang.
/// "Approximating Maximum Agreement Forest on Multiple Binary Trees".
/// Algorithmica, 2016.
/// DOI: https://doi.org/10.1007/s00453-015-0087-6
///
/// # Parameters
///
/// # `state` - Fully reduced search state.
/// # `k` - Current search bound.
///
/// # Returns
///
/// A lower bound on the order of a maximum agreement forest.
///
/// # Panics
///
/// Possibly panics if the state is not fully reduced before invocation.
///
/// # Complexity
///
/// Approximately `O(n t k)`
///
/// # Notes
///
/// Modifies the states, but restores it to a state equivalent (or equal) to the one before invocation.
pub fn calculate_simple_lb(state: &mut State, k: usize) -> usize {
    #[cfg(feature = "assert_validity")]
    state.assert_fully_reduced();

    let checkpoint = state.undos.len();
    let mut non_forced_cuts = 0;
    if state.ord() > k {
        return k + 1;
    }

    while let Some(((label1, leaf1_f2), (label2, leaf2_f2))) = state.f2().find_a_sibling_pair() {
        let is_last_tree: bool = state.is_on_last_tree();

        let ord = state.ord();
        let current_lb = ord - non_forced_cuts;
        if current_lb > k {
            // dont finish approximation if we know lb > k
            break;
        }

        if current_lb + (state.instance.num_leaves - ord) / 2 < k {
            // Don't finish approximation if we EXPECT final_lb < k:
            // If <=1/2 of the possible future cuts in the rest of the calculation being non_forced_cuts
            // is not enough to make lb>k: give up and return the current lb
            break;
            // This will occasionally fail to report lb>k when it could have,
            // but in general saves time in hopeless cases.
        }

        let [f1, f2] = state.instance.forests.get_disjoint_mut([state.f1_idx, state.f2_idx]).expect("Expected f1_idx and f2_idx to always be distinct and in bound when grabbing mutable references");
        let leaf1_f1 = f1.locate_label(label1);
        let leaf2_f1 = f1.locate_label(label2);

        let perform_leaf_cuts: bool;
        let mut num_cuts_on_f1_performed: usize = 0;
        // check if the labels are in the same components in f1
        match f1.get_lca_and_path_len(leaf1_f1, leaf2_f1) {
            Some((_, 0..=3)) => {
                unreachable!(
                    "path between two leaves cannot be shorter than 4 nodes, because we expected sibling pairs to be merged"
                );
            }
            None => {
                perform_leaf_cuts = true;
            }
            Some((lca, path_len)) => {
                // if path is length 4, there is some MAF of the 2 current trees
                // where we cut the single pendant subtree
                perform_leaf_cuts = !(path_len == 4 && is_last_tree);
                // cut one pendant subtree along the path
                'cut_one_branch_on_path: for leaf in [leaf1_f1, leaf2_f1] {
                    if let Some(parent) = f1.get(leaf).parent
                        && parent != lca
                    {
                        let [left, right] = f1.get(parent).children_unchecked();
                        if left != leaf {
                            f1.cut_branch_and_contract_with_undo(left, F1, &mut state.undos);
                        } else {
                            debug_assert!(right != leaf);
                            f1.cut_branch_and_contract_with_undo(right, F1, &mut state.undos);
                        }
                        num_cuts_on_f1_performed += 1;
                        break 'cut_one_branch_on_path;
                    }
                }
            }
        };

        if perform_leaf_cuts {
            f1.cut_branch_and_contract_with_undo(leaf1_f1, F1, &mut state.undos);
            num_cuts_on_f1_performed += 1;
            // NOTE: cutting leaf1 can make leaf2 a single vertex tree (possible when path_len=4): hence the checks
            if f1.get(leaf2_f1).parent.is_some() {
                num_cuts_on_f1_performed += 1;
                f1.cut_branch_and_contract_with_undo(leaf2_f1, F1, &mut state.undos);
            }
            f2.cut_branch_and_contract_with_undo(leaf1_f2, F2, &mut state.undos);
            if f2.get(leaf2_f2).parent.is_some() {
                f2.cut_branch_and_contract_with_undo(leaf2_f2, F2, &mut state.undos);
            }
        }
        non_forced_cuts += num_cuts_on_f1_performed - 1; // all but one of the cuts was forced

        // exhaustively apply reduction after cutting an edge to maintain invariants
        // apply_reductions(state);
        apply_reductions(state);
    }

    let lb = state.ord() - non_forced_cuts;
    state.rollback(checkpoint);
    lb
}

impl State {
    /// Computes a canonical hash representation of the current search state.
    ///
    /// The hash is defined as the sorted set of root indices (as type Idx) in forest F1.
    ///
    /// # Returns
    ///
    /// A canonical state identifier suitable for use in the visited-state cache.
    ///
    /// # Notes
    ///
    /// The uniqueness of this hash relies on the implementation of this maximal_af_dbs algorithm. In
    /// particular, on the the indices of nodes never changing and F1 always corresponding to the same input
    /// forest.
    pub fn get_hash(&self) -> StateHash {
        let mut roots: Vec<_> = self.f1().roots.iter().copied().collect();
        roots.sort_unstable();
        roots.shrink_to_fit();
        roots
    }
}
