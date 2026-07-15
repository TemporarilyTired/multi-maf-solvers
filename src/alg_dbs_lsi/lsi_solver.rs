#[cfg(not(feature = "LSI_DISABLE_LB_CALC"))]
use crate::alg_approx::calculate_simple_lb_on_general_instance;

#[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
use crate::maf_instance::arena_vertex::Label;

use crate::maf_instance::{arena_tree::ArenaTree, arena_vertex::NodeData::*, instance::Instance};

use super::arena_tree_lsi_ext::ArenaTreeLsiExt;
use super::branch_result::BranchResult::*;
use super::instance_lsi_ext::InstanceLsiDbsExt;
use super::lsi_instance::LsiInstance;
use crate::common::validity::assert_validity;

// import log macro; and only import other logging things if enabled
use super::logging::log;
#[cfg(feature = "logging")]
use super::logging::{LOGS, Logs};

#[cfg(not(feature = "LSI_DISABLE_GLOBAL_PREPROCESSING"))]
pub fn solve(instance: Instance) -> ArenaTree {
    use crate::preprocessing::reduce_and_solve;
    reduce_and_solve(instance, &solve_fn)
}
#[cfg(feature = "LSI_DISABLE_GLOBAL_PREPROCESSING")]
pub fn solve(instance: Instance) -> ArenaTree {
    solve_fn(instance)
}

/// Computes a maximum (i.e., minimum-order) agreement forest for the given instance.
///
/// This method uses the depth-bounded search algorithm in [Shi et al., 2018] (see ./mod.rs) as the subroutine until an optimal solution is found.
/// Performs binary search on the order of the optimal solution.
///
/// # Parameters
///
/// - `instance` - Input instance consisting of the forests to compute the MAF for
///
/// # Returns
///
/// A minimum-order agreement forests represented as an `ArenaTree`.
#[cfg(feature = "LSI_DISABLE_BINARY_SEARCH_ON_SOL")]
fn solve_fn(instance: Instance) -> ArenaTree {
    let lsi_instance = LsiInstance::from_sub_instance(instance);
    assert_validity!(lsi_instance);

    for k in 1..=lsi_instance.main.num_leaves {
        // reset all counters
        log!(|logs: &mut Logs| logs.reset_logs_for_new_k());

        // find an AF of ord() = k means cutting k-1 edges
        let result = dbs(lsi_instance.clone(), k);

        // print the call count logs and update totals
        log!(|logs: &Logs| logs.print_logs_after_k(k));
        log!(|logs: &mut Logs| logs.total_calls += logs.maf_calls);

        if let Some(maf) = result {
            debug_assert!(
                maf.ord() == k,
                "size of MAF must be {} (but is {}), since dbs(k') returned None for all k'<k",
                k,
                maf.ord()
            );
            log!(|logs: &Logs| logs.print_logs_full_run());
            return maf;
        }
    }
    unreachable!(
        "Solver has not found a MAF after considering k = num_leaves cuts, which should be impossible",
    )
}

/// Computes a maximum (i.e., minimum-order) agreement forest for the given instance.
///
/// This method uses the depth-bounded search algorithm in [Shi et al., 2018] (see ./mod.rs) as the subroutine until an optimal solution is found.
/// Iteratively checks k=1 to k=n until a solution is found.
///
/// # Parameters
///
/// - `instance` - Input instance consisting of the forests to compute the MAF for
///
/// # Returns
///
/// A minimum-order agreement forests represented as an `ArenaTree`.
#[cfg(not(feature = "LSI_DISABLE_BINARY_SEARCH_ON_SOL"))]
fn solve_fn(instance: Instance) -> ArenaTree {
    #[cfg(feature = "LSI_DISABLE_LB_CALC")]
    use crate::alg_approx::calculate_simple_lb_on_general_instance;

    let lsi_instance = LsiInstance::from_sub_instance(instance);
    let upper_bound_af = dbs(lsi_instance.clone(), lsi_instance.get_num_leaves())
        .expect("to find an AF of size <= #leaves");
    let mut ub = upper_bound_af.ord();
    let mut lb = calculate_simple_lb_on_general_instance(&lsi_instance.main, ub, false);
    let mut sol = upper_bound_af;

    while lb < ub {
        let mid = (lb + ub) / 2;

        // reset all counters
        log!(|logs: &mut Logs| logs.reset_logs_for_new_k());

        // find an AF of ord() = k means cutting k-1 edges
        let result = dbs(lsi_instance.clone(), mid);

        // print the call count logs and update totals
        log!(|logs: &Logs| logs.print_logs_after_k(mid));
        log!(|logs: &mut Logs| logs.total_calls += logs.maf_calls);

        if let Some(maf) = result {
            ub = maf.ord();
            sol = maf;
        } else {
            lb = mid + 1;
        }
    }
    log!(|logs: &Logs| logs.print_logs_full_run());
    sol
}

/// Finds an agreement forest of the instance with order at most max_ord if it exists.
///
/// # Parameters
///
/// - `instance` - Input instance.
/// - `max_ord` - Maximum allowed order of the resulting agreement forest.
///
/// # Returns
///
/// - `Some(solution)` if an agreement forest of order at most `max_ord` exists.
/// - `None` if no such forest exists.
#[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
pub fn solve_sub_instance(instance: Instance, max_ord: usize) -> Option<ArenaTree> {
    assert_validity!(instance);
    log!(|logs: &mut Logs| logs.subinstances_created += 1);

    let lsi_instance = LsiInstance::from_sub_instance(instance);
    for k in 0..=max_ord {
        let result = dbs(lsi_instance.clone(), k);

        if let Some(maf) = result {
            log!(|logs: &mut Logs| logs.subinstances_solved += 1);
            debug_assert_eq!(maf.ord(), k);
            return Some(maf);
        }
    }
    None
}

/// Finds an agreement forest of the instance with order at most max_ord if it exists.
///
/// This method is an implementation of the `Alg-MAF` algorithm from [Shi et al., 2018],
/// but only for binary input trees.
///
/// Some compile time features also determine the internal behaviour of the depth-bounded search:
///
/// - `LSI_DISABLE_LB_CALC`: disables lower bound calculation to prune execution branches,
/// - `LSI_DISABLE_CLUSTER_REDUCTION`: disables all cluster reduction,
/// - `LSI_DISABLE_MODIFIED_BRANCHING`: disables the significant modifications to the original branching order,
///
/// # Parameters
///
/// - `instance` - Input instance.
/// - `max_ord` - Maximum allowed order of the resulting agreement forest.
///
/// # Returns
///
/// - `Some(solution)` if an agreement forest of order at most `max_ord` exists.
/// - `None` if no such forest exists.
///
/// # Notes
///
/// The parameter max_ord is equal to 'k + ord(instance)' where k is the parameter
/// the original algorithm description uses for the 'budget'.
/// To reduce the complexity of bookkeeping this implementation
/// uses max_ord instead of keeping track of a remaining budget
pub fn dbs(mut instance: LsiInstance, max_ord: usize) -> Option<ArenaTree> {
    // increment call count logs
    log!(|logs: &mut Logs| logs.maf_calls += 1);

    assert_validity!(instance);

    // INFO: alg-maf line 1
    let starting_ord = instance.ord();
    if starting_ord > max_ord {
        return None;
    }

    // INFO: alg-maf line 2
    if !instance.satisfies_lsi_mut() {
        return br_lsi(instance, max_ord);
    }
    // INFO: alg-maf line 3
    else if !instance.f1_subgraph_of_all() {
        // exhaustively perform reduction rule 2
        // and remove SVTs
        let (performed_reductions, n_svt_removed) = instance.fully_r2_reduce_and_remove_svt();
        let mut sol = solve_r2_reduced_lsi_instance(instance, max_ord.checked_sub(n_svt_removed)?)?;

        // undo any performed reductions
        if !performed_reductions.is_empty() {
            for reduction in performed_reductions.into_iter().rev() {
                sol.undo_reduction(reduction);
            }
        }
        return Some(sol);
    }

    // INFO: alg-maf line 4
    Some(instance.extract_af())
}

/// Finds an agreement forest of the instance with order at most max_ord if it exists.
///
/// This method is the last part of the depth-bounded search procedure, solving only r2-reduced
/// instances that satisfy the LSI property
///
/// Some compile time features also determine the internal behaviour of the depth-bounded search:
///
/// - `LSI_DISABLE_LB_CALC`: disables lower bound calculation to prune execution branches,
/// - `LSI_DISABLE_CLUSTER_REDUCTION`: disables all cluster reduction,
/// - `LSI_DISABLE_MODIFIED_BRANCHING`: disables the significant modifications to the original branching order,
///
/// # Parameters
///
/// - `instance` - An r2-reduced instance satisfying the LSI property.
/// - `max_ord` - Maximum allowed order of the resulting agreement forest.
///
/// # Returns
///
/// - `Some(solution)` if an agreement forest of order at most `max_ord` exists.
/// - `None` if no such forest exists.
pub fn solve_r2_reduced_lsi_instance(
    mut instance: LsiInstance,
    max_ord: usize,
) -> Option<ArenaTree> {
    // if we have at least one separate lsi cluster: solve them separately
    #[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
    if !instance.lsi_clusters.is_empty() {
        return lsi_cluster_reduction(instance, max_ord);
    }

    debug_assert!(instance.lsi_clusters.is_empty());

    // Calculate O(ntk) lower bound for this branch (if enabled)
    #[cfg(not(feature = "LSI_DISABLE_LB_CALC"))]
    {
        let lb = calculate_simple_lb_on_general_instance(&instance.main, max_ord, true);
        if lb > max_ord {
            return None;
        }
    }

    let ord = instance.ord();
    #[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
    if ord + 4 < max_ord {
        // if the instance has a useful common cluster, split it into clusters and solve separately
        // if the addition of the dummy leaf does not increase the MAF in both parts: some component
        // can span the cut edge of the subtree
        // So: test if opt(above)+opt(below) <= k, if so: return Some(_)
        // otherwise: test if opt(above)+opt(below) > k+1, if so: return None
        // At this point opt(above)+opt(below) == k+1,
        // If opt(above with dummy) == opt(above)  (we can just check for   opt(above w dummy) <= opt(above))
        // and opt(below with dummy) == opt(below)  (we can just check for   opt(below w dummy) <= opt(below))
        // Then: return Some(opt(above) + opt(below) - 1)
        // otherwise: return None

        if let Some((
            above,
            below,
            (above_w_dummy, used_dummy_above),
            (below_w_dummy, used_dummy_below),
        )) = instance.find_clusters_w_dummy()
        {
            return solve_split_with_dummy(
                above,
                below,
                (above_w_dummy, used_dummy_above),
                (below_w_dummy, used_dummy_below),
                max_ord,
            );
        }
    }

    // NOTE: Added early return check, saving 2 instance clones and a search.
    // If ord == max_ord (i.e., budget = 0), then BR 2.1 is the same as BR 2.2.2 or BR 2.2.1,
    // so we just apply it on any cherry.
    // This is optimal because the branches where a cherry is made a SVT will have ord > max_ord.
    if ord == max_ord {
        let (_, label_a, label_b) = instance.main.forests[0]
            .iterate_cherries()
            .next()
            .expect("there must be a cherry in every tree");
        let lcas_and_paths_len = instance
            .main
            .get_lca_and_path_len_in_all_forests(label_a, label_b);

        // remove all subtrees incident to the path between label_a and label_b in all forests
        instance
            .main
            .cut_path_between_labels_in_all_forests(lcas_and_paths_len, label_a, label_b);

        return dbs(instance, max_ord);
    }

    // NOTE: If modified branching is enabled:
    // - put reduction rule 2.2.1 before any branching rules to try to save branches,
    // - otherwise always branch on the largest path possible
    #[cfg(not(feature = "LSI_DISABLE_MODIFIED_BRANCHING"))]
    {
        // INFO: Reduction Rule 2.2.1
        if instance.main.try_apply_rr2_2_1() {
            log!(|logs: &mut Logs| logs.rr_2_2_1_applied += 1);
            return dbs(instance, max_ord);
        };
        // NOTE: Branching similar to Branching Rule 2.1, but on the cherry with the largest path
        // between them in any forest
        match instance.main.apply_branching_on_largest_path(max_ord) {
            Branch3(branches) => {
                log!(|logs: &mut Logs| logs.br_2_1_applied += 1);
                return branches.into_iter().find_map(|inst| dbs(inst, max_ord));
            }
            Branch2(branches) => {
                log!(|logs: &mut Logs| logs.br_2_1_applied += 1);
                return branches.into_iter().find_map(|inst| dbs(inst, max_ord));
            }
            NotApplicable(_inst) => unreachable!("the should alway be a cherry to apply this on"),
        }
    }

    // NOTE: If modified branching is disabled:
    // apply first applicable rule from: BR 2.1, RR 2.2.1, BR 2.2.2
    #[cfg(feature = "LSI_DISABLE_MODIFIED_BRANCHING")]
    {
        // INFO: Branching Rule 2.1
        match instance.main.try_apply_br2_1(max_ord) {
            NotApplicable(inst) => {
                instance.main = inst;
            }
            Branch3(branches) => {
                log!(|logs: &mut Logs| logs.br_2_1_applied += 1);
                return branches.into_iter().find_map(|inst| dbs(inst, max_ord));
            }
            Branch2(branches) => {
                log!(|logs: &mut Logs| logs.br_2_1_applied += 1);
                return branches.into_iter().find_map(|inst| dbs(inst, max_ord));
            }
        }

        // INFO: Reduction Rule 2.2.1
        if instance.main.try_apply_rr2_2_1() {
            log!(|logs: &mut Logs| logs.rr_2_2_1_applied += 1);
            return dbs(instance, max_ord);
        };

        // we can now apply BR 2.2.2 on any cherry, and find two disagreeing forests for that label pair
        let (_, label_a, label_b) = instance.main.forests[0]
            .iterate_cherries()
            .next()
            .expect("there must be a cherry in every tree");
        let (disagreeing_forest_i, disagreeing_forest_j) = instance
            .main
            .find_disagreeing_forests(label_a, label_b)
            .expect("there must be an agreement, otherwise RR 2.2.1 would have been applied");

        // INFO: Branching Rule 2.2.2
        match instance.main.apply_br2_2_2(label_a, label_b) {
            Branch3([inst_label_a_cut, inst_label_b_cut, inst_paths_cut]) => {
                log!(|logs: &mut Logs| logs.br_2_2_2_applied += 1);
                if let Some(sol) = dbs(inst_label_a_cut, max_ord) {
                    return Some(sol);
                }
                if let Some(sol) = dbs(inst_label_b_cut, max_ord) {
                    return Some(sol);
                }

                // perform special case of BR-LSI, use the disagreeing forests as suggestions to
                // start with in the BR-PROCESS within BR-LSI
                br_lsi_with_suggestion(
                    inst_paths_cut,
                    max_ord,
                    disagreeing_forest_i,
                    disagreeing_forest_j,
                )
            }
            _ => {
                unreachable!()
            }
        }
    }
}

/// Finds an agreement forest of the instance with order at most max_ord if it exists.
///
/// This method is an implementation of the `BR-LSI` process from [Shi et al., 2018]
///
/// # Parameters
///
/// - `instance` - Input instance.
/// - `max_ord` - Maximum allowed order of the resulting agreement forest.
///
/// # Returns
///
/// - `Some(solution)` if an agreement forest of order at most `max_ord` exists.
/// - `None` if no such forest exists.
pub fn br_lsi(instance: LsiInstance, max_ord: usize) -> Option<ArenaTree> {
    log!(|logs: &mut Logs| logs.br_lsi_calls += 1);

    // INFO: br-lsi line 1
    if instance.ord() > max_ord {
        return None;
    }
    br_process(instance, max_ord)
}

/// Finds an agreement forest of the instance with order at most max_ord if it exists.
///
/// This method is an implementation of the `BR-PROCESS` from [Shi et al., 2018]
///
/// Some compile time features also determine the internal behaviour of the depth-bounded search:
///
/// - `LSI_DISABLE_CLUSTER_REDUCTION`: disables all cluster reduction.
///
/// # Parameters
///
/// - `instance` - Input instance.
/// - `max_ord` - Maximum allowed order of the resulting agreement forest.
///
/// # Returns
///
/// - `Some(solution)` if an agreement forest of order at most `max_ord` exists.
/// - `None` if no such forest exists.
///
/// # Notes
///
/// The BR-PROCESS is modified as described in the descriptions of BR-LSI and Alg-MAF from [Shi et al., 2018].
pub fn br_process(mut instance: LsiInstance, max_ord: usize) -> Option<ArenaTree> {
    log!(|logs: &mut Logs| logs.br_process_calls += 1);

    // INFO: BR-PROCESS line 1
    instance.main.fully_r1_reduce();
    assert_validity!(instance);

    // INFO: BR-LSI line 2a
    let new_ord = instance.ord();
    if new_ord > max_ord {
        return None;
    }

    // INFO: BR-PROCESS line 2
    if !instance.satisfies_lsi_mut() {
        // If cluster reduction is enabled: try it before defaulting to BR2-PROCESS
        #[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
        {
            // NOTE: eager cluster reduction (when some but not all clusters satisfy LSI)

            // the r2 and svt reductions avoid performing cluster reduction when clusters are isomorphic
            let (performed_reductions, n_svt_removed) = instance.fully_r2_reduce_and_remove_svt();
            let new_max_ord = max_ord.checked_sub(n_svt_removed)?;

            // if there is at least one LSI component and at least two in total:
            // solve the LSI components separately from each other and the rest
            let mut sol = if !instance.lsi_clusters.is_empty() {
                lsi_cluster_reduction(instance, new_max_ord)?
            } else {
                // otherwise enter the BR2-PROCESS
                let (i, j) = instance.main.find_max_ord_pair_without_lsi();
                br2_process(instance, i, j, new_max_ord, true)?
            };

            for reduction in performed_reductions.into_iter().rev() {
                sol.undo_reduction(reduction);
            }
            Some(sol)
        }

        // enter the BR2-PROCESS
        #[cfg(feature = "LSI_DISABLE_CLUSTER_REDUCTION")]
        {
            let (i, j) = instance.main.find_max_ord_pair_without_lsi();
            br2_process(instance, i, j, max_ord, true)
        }
    } else {
        // INFO: BR-PROCESS line 3, modified by BR-LSI line 2b and Alg-MAF line 2
        dbs(instance, max_ord)
    }
}

/// Finds an agreement forest of the instance with order at most max_ord if it exists.
///
/// This method is an implementation of the `2BR-PROCESS` from [Shi et al., 2018]
///
/// # Parameters
///
/// - `instance` - Input instance.
/// - `i` - Index of one of the two forests to apply the process on.
/// - `j` - Index of one of the two forests to apply the process on.
/// - `max_ord` - Maximum allowed order of the resulting agreement forest.
/// - `is_r1_reduced` - Whether the instance is already r1-reduced before invocation.
///
/// # Returns
///
/// - `Some(solution)` if an agreement forest of order at most `max_ord` exists.
/// - `None` if no such forest exists.
///
/// # Notes
///
/// The `is_r1_reduced` parameter is introduced in this implementation to prevent
/// r1-reduction from being called when the instance is known to be fully r1-reduced.
pub fn br2_process(
    mut instance: LsiInstance,
    i: usize,
    j: usize,
    max_ord: usize,
    is_r1_reduced: bool,
) -> Option<ArenaTree> {
    log!(|logs: &mut Logs| logs.br2_process_calls += 1);

    let [f_i, f_j] = instance
        .main
        .forests
        .get_disjoint_mut([i, j])
        .expect("unique indices to forests i and j");

    // If we are coming from BR-PROCESS, the pair is already r1-reduced
    if !is_r1_reduced {
        // INFO: 2BR_PROCESS line 1
        ArenaTree::r1_reduce_pair(f_i, f_j);
    }

    // NOTE: Extra pruning introduced here because otherwise an instance might recursively call `br2_process`
    // many times before noticing that the `max_ord` was exceeded.
    let new_ord = f_i.ord().max(f_j.ord());
    if new_ord > max_ord {
        log!(|logs: &mut Logs| logs.early_returns_br2_process_after_r1 += 1);
        return None;
    }

    assert_validity!(f_i);
    assert_validity!(f_j);

    // INFO: 2BR_PROCESS line 2.x
    if let Some(v) = f_i.find_br1_target_with(f_j) {
        let Internal { left, right } = f_i.get(v).data else {
            unreachable!("target for branching rule 1 must be an internal node");
        };
        let mut new_instance = instance.clone();
        new_instance.main.forests[i].cut_branch(left);
        if let Some(res) = br2_process(new_instance, i, j, max_ord, false) {
            return Some(res);
        }

        instance.main.forests[i].cut_branch(right);
        br2_process(instance, i, j, max_ord, false)
    } else if let Some(v) = f_j.find_br1_target_with(f_i) {
        let Internal { left, right } = f_j.get(v).data else {
            unreachable!("target for branching rule 1 must be an internal node");
        };
        let mut new_instance = instance.clone();
        new_instance.main.forests[j].cut_branch(left);
        if let Some(res) = br2_process(new_instance, i, j, max_ord, false) {
            return Some(res);
        }

        instance.main.forests[j].cut_branch(right);
        br2_process(instance, i, j, max_ord, false)
    } else {
        debug_assert!(&instance.main.forests[i].has_same_label_sets(&instance.main.forests[j]));

        // INFO: 2BR_PROCESS line 3
        br_process(instance, max_ord)
    }
}

/// Finds an agreement forest of the instance with order at most max_ord if it exists.
///
/// This method is an implementation of `BR-LSI` with the modification
/// from branching rule 2.2.2, branch 3, from [Shi et al., 2018].
///
/// # Parameters
///
/// - `instance` - Input instance.
/// - `max_ord` - Maximum allowed order of the resulting agreement forest.
/// - `suggested_forest_i` - Index of one of the two forests suggested to try first in the br2_process calls.
/// - `suggested_forest_j` - Index of one of the two forests suggested to try first in the br2_process calls.
///
/// # Returns
///
/// - `Some(solution)` if an agreement forest of order at most `max_ord` exists.
/// - `None` if no such forest exists.
#[cfg(feature = "LSI_DISABLE_MODIFIED_BRANCHING")]
pub fn br_lsi_with_suggestion(
    mut instance: LsiInstance,
    max_ord: usize,
    suggested_forest_i: usize,
    suggested_forest_j: usize,
) -> Option<ArenaTree> {
    log!(|logs: &mut Logs| logs.br_process_calls += 1);
    log!(|logs: &mut Logs| logs.br_2_2_2_special_case += 1);

    // INFO: BR-LSI line 1
    let ord = instance.ord();
    if ord > max_ord {
        return None;
    }

    // INFO: BR-PROCESS line 1
    instance.main.fully_r1_reduce();
    assert_validity!(instance);

    // INFO: BR-LSI line 2a
    let new_ord = instance.ord();
    if new_ord > max_ord {
        return None;
    }

    // INFO: BR-PROCESS line 2
    if !instance.satisfies_lsi_mut() {
        // If cluster reduction is enabled: try it before defaulting to BR2-PROCESS
        #[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
        {
            // NOTE: eager cluster reduction (when some but not all clusters satisfy LSI)

            // the r2 and svt reductions avoid performing cluster reduction when clusters are isomorphic
            let (performed_reductions, n_svt_removed) = instance.fully_r2_reduce_and_remove_svt();
            let new_max_ord = max_ord.checked_sub(n_svt_removed)?;

            // if there is at least one LSI component and at least two in total:
            // solve the LSI components separately from each other and the rest
            let mut sol = if !instance.lsi_clusters.is_empty() {
                lsi_cluster_reduction(instance, new_max_ord)?
            } else {
                // Otherwise enter the BR2-PROCESS
                // Use the suggested forests for the subroutine if they have maximum order
                let (i, j) = if instance.main.forests[suggested_forest_i]
                    .ord()
                    .max(instance.main.forests[suggested_forest_j].ord())
                    == instance.main.ord()
                {
                    (suggested_forest_i, suggested_forest_j)
                } else {
                    instance.main.find_max_ord_pair_without_lsi()
                };
                br2_process(instance, i, j, new_max_ord, true)?
            };

            for reduction in performed_reductions.into_iter().rev() {
                sol.undo_reduction(reduction);
            }
            Some(sol)
        }

        #[cfg(feature = "LSI_DISABLE_CLUSTER_REDUCTION")]
        {
            // Enter the BR2-PROCESS
            // Use the suggested forests for the subroutine if they have maximum order
            let (i, j) = if instance.main.forests[suggested_forest_i]
                .ord()
                .max(instance.main.forests[suggested_forest_j].ord())
                == instance.main.ord()
            {
                (suggested_forest_i, suggested_forest_j)
            } else {
                instance.main.find_max_ord_pair_without_lsi()
            };
            br2_process(instance, i, j, max_ord, true)
        }
    } else {
        // INFO: BR-PROCESS line 3, modified by BR-LSI line 2b and Alg-MAF line 2
        dbs(instance, max_ord)
    }
}

/// Finds an agreement forest of the instance with order at most max_ord if it exists.
/// Performs LSI cluster reduction, which solves the separate clusters before recombinig them into an
/// agreement forest for the main problem.
///
/// # Parameters
///
/// - `{main, lsi_clusters}` - The input instance where main is any cluster and lsi_clusters in a
/// collection of the clusters extracted from the instance.
/// - `max_ord` - Maximum allowed order of the resulting agreement forest.
///
/// # Returns
///
/// - `Some(solution)` if an agreement forest of order at most `max_ord` exists.
/// - `None` if no such forest exists.
#[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
fn lsi_cluster_reduction(
    LsiInstance { main, lsi_clusters }: LsiInstance,
    max_ord: usize,
) -> Option<ArenaTree> {
    #[cfg(feature = "LSI_DISABLE_LB_CALC")]
    use crate::alg_approx::calculate_simple_lb_on_general_instance;

    log!(|logs: &mut Logs| logs.lsi_cluster_reduc += 1);

    // Calculate lower bounds for all other clusters
    let lbs_of_clusters: Vec<usize> = lsi_clusters
        .iter()
        .map(|cluster| calculate_simple_lb_on_general_instance(cluster, max_ord, false))
        .collect();
    let mut remaining_lb_total = lbs_of_clusters.iter().sum();
    let remaining_budget = max_ord.checked_sub(remaining_lb_total)?;
    // Use the lower bound to decrease the remaining max_ord for the main cluster
    let mut sol = solve_sub_instance(main, remaining_budget)?;

    let mut n_lsi_clusters_to_do = lsi_clusters.len();

    debug_assert_ne!(n_lsi_clusters_to_do, 0);
    for (i, lsi_cluster) in lsi_clusters.into_iter().enumerate() {
        if n_lsi_clusters_to_do >= 2 {
            // Use the lower bounds for the upcoming and the solutions for the previous clusters
            // To lower the budget (i.e., search depth)
            remaining_lb_total -= lbs_of_clusters[i];
            let remaining_budget = max_ord.checked_sub(sol.ord() + remaining_lb_total)?;

            let cluster_sol = solve_sub_instance(lsi_cluster, remaining_budget)?;
            sol = sol.join_with(cluster_sol);
            n_lsi_clusters_to_do -= 1;
        } else {
            let dbs_sol = dbs(
                LsiInstance::from_sub_instance(lsi_cluster),
                max_ord.checked_sub(sol.ord())?,
            )?;
            sol = sol.join_with(dbs_sol);
        }
    }

    Some(sol)
}

/// Finds an agreement forest of the instance with order at most max_ord if it exists.
/// Performs subtree cluster reduction, on a detached subtree with a common label-set,
/// and solves the subtree cluster separately from the main instance,
/// before recombinig them into an agreement forest for the main problem.
///
/// See the post https://skelk.sdf-eu.org/cluster/ by Kelk for a
/// pedagogical implementation and explanation of this cluster redcution.
///
/// # Parameters
///
/// - `above` - The main input instance, where the `below` cluster is removed.
/// - `below` - The extracted subtree cluster.
/// - `(above_with_dummy, used_dummy_above)` - The main input instance, but with a dummy
/// (used_dummy_above) added.
/// - `(below_with_dummy, used_dummy_below)` - The subtree cluster, but with a dummy
/// (used_dummy_below) added.
/// - `max_ord` - Maximum allowed order of the resulting agreement forest.
///
/// # Returns
///
/// - `Some(solution)` if an agreement forest of order at most `max_ord` exists.
/// - `None` if no such forest exists.
#[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
pub fn solve_split_with_dummy(
    above: LsiInstance,
    below: LsiInstance,
    (above_with_dummy, used_dummy_above): (LsiInstance, Label),
    (below_with_dummy, used_dummy_below): (LsiInstance, Label),
    max_ord: usize,
) -> Option<ArenaTree> {
    #[cfg(feature = "LSI_DISABLE_LB_CALC")]
    use crate::alg_approx::calculate_simple_lb_on_general_instance;

    debug_assert!(above.lsi_clusters.is_empty());
    debug_assert!(below.lsi_clusters.is_empty());

    log!(|logs: &mut Logs| logs.subtree_cluster_reduc += 1);

    let lb_a = calculate_simple_lb_on_general_instance(&above.main, max_ord, false);

    let mut sol_b = dbs(below.clone(), (max_ord + 1).checked_sub(lb_a)?)?;
    let mut ub_b = sol_b.ord();
    let budget_a = (max_ord + 1).saturating_sub(ub_b);
    let mut maybe_sol_a = dbs(above.clone(), budget_a);

    // repeatedly decrease the upper bound on b
    // and try to find a solution to a with total ord <= k+1
    while maybe_sol_a.is_none() {
        sol_b = dbs(below.clone(), ub_b.checked_sub(1)?)?;
        ub_b = sol_b.ord();
        maybe_sol_a = dbs(above.clone(), (max_ord + 1).saturating_sub(ub_b));
    }
    // we have found a solution of ord <= k+1

    let sol_a = maybe_sol_a.unwrap();
    let sol_ord_a = sol_a.ord();
    let sol_ord_b = ub_b as usize;

    debug_assert_eq!(sol_ord_b, sol_b.ord());

    // if we are lucky: we have a solution of <= k already
    if sol_ord_a + sol_ord_b <= max_ord {
        return Some(sol_a.join_with(sol_b));
    }

    // we found a solution of size exactly k+1: try to decrease it to k
    if let Some(smaller_sol_b) = dbs(below.clone(), sol_ord_b - 1) {
        return Some(sol_a.join_with(smaller_sol_b));
    }
    if let Some(smaller_sol_a) = dbs(above.clone(), sol_ord_a - 1) {
        return Some(smaller_sol_a.join_with(sol_b));
    }

    // decreasing to k failed: see if clusters with dummies have the same optimal
    // if so: join them into a solution of size k
    let sol_a_w_dummy = dbs(above_with_dummy, sol_ord_a)?;
    let sol_b_w_dummy = dbs(below_with_dummy, sol_ord_b)?;

    Some(sol_a_w_dummy.join_at_dummy(sol_b_w_dummy, used_dummy_above, used_dummy_below))
}
