use crate::alg_dbs_lsi::arena_tree_lsi_ext::ArenaTreeLsiExt;
use crate::alg_dbs_lsi::instance_lsi_ext::InstanceLsiDbsExt;
use crate::maf_instance::instance::Instance;

/// Calculate a simple lower bound on the number of cuts needed to consolidate the currently
/// loaded forests (which is a lower bound for the entire branch).
/// Assumes instance is fully reduced before invocation
///
/// # Panics
///
/// Possibly panics if forest is not fully reduced before call
pub fn calculate_simple_lb_on_general_instance(
    instance: &Instance,
    k: usize,
    terminate_early: bool,
) -> usize {
    let mut f2_idx: usize = 1;

    let mut two_forest_instance = Instance {
        num_leaves: instance.num_leaves,
        forests: vec![
            instance.forests[0].clone(),
            instance.forests[f2_idx].clone(),
        ],
    };

    let mut non_forced_cuts: usize = 0;

    let old_ord = two_forest_instance.forests[0].ord();
    if old_ord > k {
        return k + 1;
    }

    let (mut performed_reductions, mut n_removed_svts) =
        two_forest_instance.fully_r2_reduce_and_remove_svt();

    loop {
        let [f1, f2] = two_forest_instance
            .forests
            .get_disjoint_mut([0, 1])
            .unwrap();

        let ord = f1.ord();
        let current_lb = ord + n_removed_svts - non_forced_cuts;
        if current_lb > k {
            // dont finish approximation if we know lb > k
            break;
        }

        if terminate_early {
            let approximated_lb_increase = (instance.num_leaves - ord) / 2;
            if current_lb + approximated_lb_increase < k {
                // Don't finish approximation if we EXPECT final_lb < k:
                // If <=1/2 of the possible future cuts in the rest of the calculation being non_forced_cuts
                // is not enough to make lb>k: give up and return the current lb
                break;
                // This will occasionally fail to report lb>k when it could have,
                // but in general saves time in hopeless cases.
            }
        }

        let is_last_tree: bool = f2_idx == instance.forests.len() - 1;
        let Some((_, label1, label2)) = f2.iterate_cherries().next() else {
            for performed_reduction in performed_reductions.into_iter().rev() {
                f1.undo_reduction(performed_reduction);
                two_forest_instance.num_leaves += 1;
            }
            n_removed_svts = 0;

            f2_idx += 1;
            if f2_idx >= instance.forests.len() {
                break;
            }

            two_forest_instance.forests[1] = instance.forests[f2_idx].clone();
            (performed_reductions, n_removed_svts) =
                two_forest_instance.fully_r2_reduce_and_remove_svt();
            continue;
        };

        let leaf1_f1 = f1.locate_label(label1);
        let leaf2_f1 = f1.locate_label(label2);

        // check if the labels are in the same components in F1
        // if so: cut any pendant subtree between them
        let perform_leaf_cuts: bool;
        let mut num_cuts_on_f1_performed: usize = 0;
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
                            f1.cut_branch(left);
                        } else {
                            debug_assert!(right != leaf);
                            f1.cut_branch(right);
                        }
                        num_cuts_on_f1_performed += 1;
                        break 'cut_one_branch_on_path;
                    }
                }
            }
        };

        // if required: make the labels single vertex trees in both F1 and F2
        if perform_leaf_cuts {
            f1.make_svt(label1);
            num_cuts_on_f1_performed += 1;
            // NOTE: cutting leaf1 can make leaf2 a single vertex tree (possible when path_len=4): hence the checks
            if f1.get(f1.locate_label(label2)).parent.is_some() {
                num_cuts_on_f1_performed += 1;
                f1.make_svt(label2);
            }

            f2.make_svt(label1);
            if f2.get(f2.locate_label(label2)).parent.is_some() {
                f2.make_svt(label2);
            }
        }
        non_forced_cuts += num_cuts_on_f1_performed - 1;

        // exhaustively apply reduction after cutting an edge to maintain invariants
        let (new_reductions, new_removed_svts) =
            two_forest_instance.fully_r2_reduce_and_remove_svt();
        performed_reductions.extend(new_reductions);
        n_removed_svts += new_removed_svts;
    }

    // we have an upper bound:  ord + n_removed_svts
    // and a lower bound:       ord + n_removed_svts - non_forced_cuts
    two_forest_instance.forests[0].ord() + n_removed_svts - non_forced_cuts
}
