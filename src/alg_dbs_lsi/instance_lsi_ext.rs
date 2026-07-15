use std::collections::HashSet;

use crate::common::validity::assert_validity;
use crate::maf_instance::{
    arena_tree::ArenaTree,
    arena_vertex::{Idx, Label, NodeData::*, Status::*},
    instance::Instance,
    performed_reduction::PerformedReduction::{self, *},
    tree_traversal::TreeTraversal,
};

use super::arena_tree_lsi_ext::ArenaTreeLsiExt;
use super::branch_result::BranchResult::{self, *};
use super::lsi_instance::LsiInstance;
use super::r1_workspace::R1Workspace;

pub trait InstanceLsiDbsExt {
    #[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
    fn extract_lsi_cluster_with_label(&mut self, label: Label) -> Instance;

    fn fully_r1_reduce(&mut self) -> bool;

    /// Exhaustively performs reductions:
    /// - merges common subtrees between all forests
    /// - syncs single vertex trees between the forests and removes them
    ///
    /// returns the removed svts and the merged labels
    fn fully_r2_reduce_and_remove_svt(&mut self) -> (Vec<PerformedReduction>, usize);

    fn find_max_ord_pair_without_lsi(&self) -> (usize, usize);

    fn try_apply_br2_1(self, max_ord: usize) -> BranchResult;

    /// Takes two labels that form a cherry in one of the forests
    /// and finds the indices of two forests for which the labels set of the
    /// lca of the labels is different (if they exist, otherwise None)
    fn find_disagreeing_forests(&self, label_a: Label, label_b: Label) -> Option<(usize, usize)>;

    /// Applies reduction rule 2.2.1 on cherry if applicable
    /// Assumes instance satisfies LSI
    /// If not applied: returns a tuple of forests for which the label set of the lca
    /// does not match
    fn try_apply_rr2_2_1(&mut self) -> bool;

    /// Applies branching rule 2.2.2 on cherry (a,b)
    /// Assumes instance satisfies LSI
    fn apply_br2_2_2(self, label_a: Label, label_b: Label) -> BranchResult;

    #[cfg(not(feature = "LSI_DISABLE_MODIFIED_BRANCHING"))]
    fn apply_branching_on_largest_path(self, max_ord: usize) -> BranchResult;

    /// Cuts off the leaf with the specified label in all forests
    ///NOTE: Panics if any leaf already is a single vertex tree
    /// This can only happen if instance did not satisfy LSI
    fn cut_label_in_all_forests(&mut self, label: Label);

    /// Cuts off all subtrees pendant to the path between two labels
    /// in all forests
    /// Takes in a vector with precomputed LCAs and path lengths in each forest
    fn cut_path_between_labels_in_all_forests(
        &mut self,
        lcas_and_paths_len: Vec<(Idx, usize)>,
        label_a: Label,
        label_b: Label,
    );

    /// Takes an instance satisfying the LSI and two labels
    /// and finds the length of the path between them in each forest
    /// and finds the lca in each forest
    ///
    /// # Panics if labels are not in the same component in some forest
    fn get_lca_and_path_len_in_all_forests(
        &self,
        label_a: Label,
        label_b: Label,
    ) -> Vec<(Idx, usize)>;

    /// Tests if f1 is a subgraph of all other forests.
    /// Note: if instance satisfies LSI then it is an AF
    fn f1_subgraph_of_all(&self) -> bool;

    fn satisfies_lsi(&self) -> bool;
}

impl InstanceLsiDbsExt for Instance {
    #[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
    fn extract_lsi_cluster_with_label(&mut self, label: Label) -> Instance {
        let mut lsi_cluster_forests = vec![];
        for forest in self.forests.iter_mut() {
            let lsi_comp = forest.find_root_of(forest.locate_label(label));
            let cluster = forest.extract_lsi_comp(lsi_comp);
            lsi_cluster_forests.push(cluster);
        }
        let n_cluster_labels = lsi_cluster_forests[0].leaf_map.len();
        self.num_leaves -= n_cluster_labels;
        Instance {
            num_leaves: n_cluster_labels,
            forests: lsi_cluster_forests,
        }
    }

    fn fully_r1_reduce(&mut self) -> bool {
        let n_forests = self.forests.len();
        let n_nodes = self
            .forests
            .iter()
            .map(|f| f.arena.len())
            .max()
            .unwrap_or(0);
        let max_label = self.forests[0].max_label_value() as usize;
        let mut done: Vec<Vec<bool>> = vec![vec![false; n_forests]; n_forests];
        let workspace: &mut R1Workspace = &mut R1Workspace::new(max_label, n_nodes);

        let mut has_performed_cut = false;

        loop {
            let mut any_changed = false;

            for i in 1..n_forests {
                for j in 0..i {
                    if done[i][j] {
                        continue;
                    }

                    let [f_i, f_j] = self.forests.get_disjoint_mut([i, j]).unwrap();
                    let (i_changed, j_changed) =
                        ArenaTree::r1_reduce_pair_w_workspace(f_i, f_j, workspace);

                    // mark pair as fully reduced
                    done[i][j] = true;

                    // if a forest changes: invalidate all 'done' pairs involving it
                    if i_changed {
                        any_changed = true;
                        for k in 0..n_forests {
                            if k != j {
                                let (a, b) = if i > k { (i, k) } else { (k, i) };
                                done[a][b] = false;
                            }
                        }
                        debug_assert!(done[i][j]);
                    }
                    if j_changed {
                        any_changed = true;
                        for k in 0..n_forests {
                            if k != i {
                                let (a, b) = if j > k { (j, k) } else { (k, j) };
                                done[a][b] = false;
                            }
                        }
                        debug_assert!(done[i][j]);
                    }
                }
            }

            if !any_changed {
                break;
            }
            has_performed_cut |= any_changed;
        }
        has_performed_cut
    }

    /// Exhaustively performs reductions:
    /// - merges common subtrees between all forests
    /// - syncs single vertex trees between the forests and removes them
    ///
    /// returns the removed svts and the merged labels
    fn fully_r2_reduce_and_remove_svt(&mut self) -> (Vec<PerformedReduction>, usize) {
        let mut performed_reductions: Vec<PerformedReduction> = vec![];
        let mut n_removed_svts: usize = 0;

        let mut to_check: Vec<Label> = self.forests[0].leaf_map.keys().copied().collect();
        let mut to_check_set: HashSet<Label> = to_check.iter().copied().collect();

        'outer: while let Some(label) = to_check.pop() {
            to_check_set.remove(&label);

            // check if label is part of a common cherry
            // or if it is a svt that is not yet synced
            let f1 = &self.forests[0];
            let Some(leaf_idx) = f1.try_locate_label(label) else {
                continue 'outer;
            };

            let v = f1.get(leaf_idx);
            debug_assert_eq!(v.status, Present);

            'check_common_cherry: {
                let Some(parent_idx) = v.parent else {
                    break 'check_common_cherry;
                };

                let parent = f1.get(parent_idx);

                let Internal { left, right } = parent.data else {
                    break 'check_common_cherry;
                };

                let Leaf { label: label_left } = f1.get(left).data else {
                    break 'check_common_cherry;
                };
                let Leaf { label: label_right } = f1.get(right).data else {
                    break 'check_common_cherry;
                };

                // check if cherry (left, right) is common to all forests
                for forest in &self.forests[1..] {
                    let leaf1 = forest.locate_label(label_left);
                    let leaf2 = forest.locate_label(label_right);
                    if forest.get(leaf1).parent.is_none()
                        || forest.get(leaf1).parent != forest.get(leaf2).parent
                    {
                        break 'check_common_cherry;
                    }
                }

                // (arbitrarily) assign the new label to be the one with the lowest number
                let new_label = label_left.min(label_right);

                // apply merge in each forest
                self.merge_common_sibling(label_left, label_right, new_label);
                performed_reductions.push(LabelsMerged {
                    original1: label_left,
                    original2: label_right,
                    new_label,
                });

                // NOTE: all cut-opts containing the old-label and new-label should be
                // outdated cut-opts where the other is a s.v.t

                // queue up the new-label for (possibly) syncing svt's
                // or merging the (possibly) newly formed common cherry
                if to_check_set.insert(new_label) {
                    to_check.push(new_label);
                }
                continue 'outer;
            }

            'check_svt: {
                let mut is_svt_in_any = false;
                let mut is_svt_in_all = true;
                for f in self.forests.iter_mut() {
                    let label_idx = f.locate_label(label);
                    let is_svt_in_this_forests = f.get(label_idx).parent.is_none();

                    is_svt_in_all &= is_svt_in_this_forests;
                    is_svt_in_any |= is_svt_in_this_forests;
                }
                match (is_svt_in_any, is_svt_in_all) {
                    (_, true) => {
                        // label is already an svt in every forest
                    }
                    (false, false) => {
                        // label not an svt in any forest
                        break 'check_svt;
                    }
                    (true, false) => {
                        // label is an svt in some forests, but not all:
                        //  - make it an svt in all forests
                        let affected_sibling_labels = self.cut_svt_return_siblings_labels(label);

                        for sibling_label in affected_sibling_labels {
                            // queue up the all affected sibling labels to check them again
                            if to_check_set.insert(sibling_label) {
                                to_check.push(sibling_label);
                            }
                        }
                    }
                }

                // label is an svt in every forest: remove it from the instance
                self.remove_svt(label);
                n_removed_svts += 1;
                performed_reductions.push(SvtRemoved { label });
            }
        }

        (performed_reductions, n_removed_svts)
    }

    fn find_max_ord_pair_without_lsi(&self) -> (usize, usize) {
        let (i, max_ord_forest) = self
            .forests
            .iter()
            .enumerate()
            .max_by_key(|&(_, forest)| forest.ord())
            .expect("non-empty instance");
        for j in (0..i).chain(i + 1..) {
            if !max_ord_forest.has_same_label_sets(&self.forests[j]) {
                return (i, j);
            }
        }

        unreachable!("expected instance not to satisfy LSI property");
    }

    fn try_apply_br2_1(mut self, max_ord: usize) -> BranchResult {
        assert_validity!(self);
        let target = {
            self.forests[0]
                .iterate_cherries()
                .find_map(|(_, label_l, label_r)| {
                    let lcas_and_paths_len =
                        self.get_lca_and_path_len_in_all_forests(label_l, label_r);
                    if lcas_and_paths_len
                        .iter()
                        .any(|&(_, path_len)| path_len >= 5)
                    {
                        // some path length is >= 5  (thus at least two pendant subtrees to be cut)
                        // so conditions for case 2.1 are met for this cherry
                        return Some((label_l, label_r, lcas_and_paths_len));
                    }
                    None
                })
        };

        let Some((label_l, label_r, lcas_and_paths_len)) = target else {
            return NotApplicable(self);
        };

        for f in self.forests.iter_mut() {
            f.roots.reserve(1);
        }

        // apply branches
        // branch [1]: remove edge incident to label_l in all forests
        let mut branch1_instance = self.clone();
        branch1_instance.cut_label_in_all_forests(label_l);

        let budget = max_ord as i32 - self.ord() as i32;
        // check if all forests will stay within budget; if not: skip branch [3]
        if lcas_and_paths_len
            .iter()
            .all(|&(_, path_len)| path_len as i32 - 3 <= budget)
        {
            // branch [2]: remove edge incident to label_r in all forests
            let mut branch2_instance = self.clone();
            branch2_instance.cut_label_in_all_forests(label_r);

            // branch [3]: remove all subtrees incident to the path between label_l and label_r in all forests
            self.cut_path_between_labels_in_all_forests(lcas_and_paths_len, label_l, label_r);

            assert_validity!(branch1_instance);
            assert_validity!(branch2_instance);
            assert_validity!(self);
            return Branch3([
                LsiInstance::from_sub_instance(branch1_instance),
                LsiInstance::from_sub_instance(branch2_instance),
                LsiInstance::from_sub_instance(self),
            ]);
        }

        // branch [2]: remove edge incident to label_r in all forests
        self.cut_label_in_all_forests(label_r);
        assert_validity!(branch1_instance);
        assert_validity!(self);

        Branch2([
            LsiInstance::from_sub_instance(branch1_instance),
            LsiInstance::from_sub_instance(self),
        ])
    }

    #[cfg(not(feature = "LSI_DISABLE_MODIFIED_BRANCHING"))]
    fn apply_branching_on_largest_path(mut self, max_ord: usize) -> BranchResult {
        assert_validity!(self);
        let target = self.forests[0]
            .iterate_cherries()
            .map(|(_, label_l, label_r)| {
                let lcas_and_paths_len = self.get_lca_and_path_len_in_all_forests(label_l, label_r);
                (label_l, label_r, lcas_and_paths_len)
            })
            .max_by_key(|(_, _, lcas_and_paths_len)| {
                lcas_and_paths_len
                    .iter()
                    .map(|(_, path_len)| *path_len)
                    .max()
                    .unwrap_or_default()
            });

        let Some((label_l, label_r, lcas_and_paths_len)) = target else {
            unreachable!("expect branching rule to be applicable: expect at least one cherry")
        };

        for f in self.forests.iter_mut() {
            f.roots.reserve(1);
        }

        // apply branches
        // branch [1]: remove edge incident to label_l in all forests
        let mut branch1_instance = self.clone();
        branch1_instance.cut_label_in_all_forests(label_l);

        let budget = max_ord as i32 - self.ord() as i32;
        // check if all forests will stay within budget; if not: skip branch [3]
        if lcas_and_paths_len
            .iter()
            .all(|&(_, path_len)| path_len as i32 - 3 <= budget)
        {
            // branch [2]: remove edge incident to label_r in all forests
            let mut branch2_instance = self.clone();
            branch2_instance.cut_label_in_all_forests(label_r);

            // branch [3]: remove all subtrees incident to the path between label_l and label_r in all forests
            self.cut_path_between_labels_in_all_forests(lcas_and_paths_len, label_l, label_r);

            assert_validity!(branch1_instance);
            assert_validity!(branch2_instance);
            assert_validity!(self);
            return Branch3([
                LsiInstance::from_sub_instance(branch1_instance),
                LsiInstance::from_sub_instance(branch2_instance),
                LsiInstance::from_sub_instance(self),
            ]);
        }

        // branch [2]: remove edge incident to label_r in all forests
        self.cut_label_in_all_forests(label_r);
        assert_validity!(branch1_instance);
        assert_validity!(self);

        Branch2([
            LsiInstance::from_sub_instance(branch1_instance),
            LsiInstance::from_sub_instance(self),
        ])
    }

    /// Takes two labels that form a cherry in one of the forests
    /// and finds the indices of two forests for which the labels set of the
    /// lca of the labels is different (if they exist, otherwise None)
    fn find_disagreeing_forests(&self, label_a: Label, label_b: Label) -> Option<(usize, usize)> {
        if self.forests.len() == 2 {
            return None;
        }
        let lcas_and_paths_len = self.get_lca_and_path_len_in_all_forests(label_a, label_b);

        let mut iter = lcas_and_paths_len
            .iter()
            .zip(&self.forests)
            .enumerate()
            .filter(|(_, ((_, path_len), _))| *path_len > 3);

        // actually in a fully r2-reduced instance there must always be some forest with path_len > 3
        let (ref_idx, ((ref_lca, _), ref_forest)) =
            iter.next().expect("expected fully r2-reduced instance");

        let reference_labels: HashSet<Label> = ref_forest.dfs_from(*ref_lca).labels().collect();

        for (i, ((lca, path_len), forest)) in iter {
            if *path_len <= 3 {
                continue;
            }
            let mut count = 0;

            for label in forest.dfs_from(*lca).labels() {
                if !reference_labels.contains(&label) {
                    return Some((ref_idx, i));
                }
                count += 1;
            }
            if count != reference_labels.len() {
                return Some((ref_idx, i));
            }
        }
        None
    }

    /// Applies reduction rule 2.2.1 on cherry if applicable
    /// Assumes instance satisfies LSI
    /// If not applied: returns a tuple of forests for which the label set of the lca
    /// does not match
    // pub fn try_apply_rr2_2_1(&mut self) -> Option<(usize, usize)> {
    fn try_apply_rr2_2_1(&mut self) -> bool {
        fn rr2_2_1_applicable(
            instance: &Instance,
            lcas_and_paths_len: &[(Idx, usize)],
            label_set: &mut [u8],
            generation: u8,
        ) -> bool {
            let mut iter = lcas_and_paths_len
                .iter()
                .zip(&instance.forests)
                .filter(|((_, path_len), _)| *path_len > 3);

            // in a fully r2-reduced instance there must always be some forest with path_len > 3
            let (&(ref_lca, ref_path_len), ref_forest) =
                iter.next().expect("expected fully r2-reduced instance");
            if ref_path_len >= 5 {
                // the rr 2.2.1 rule is only optimal with one subtree on the path between the
                // labels, i.e, a path length of at most 4 in total
                return false;
            }

            // store label set using generation index
            let mut ref_count = 0;
            for Label(label) in ref_forest.dfs_from(ref_lca).labels() {
                label_set[label as usize] = generation;
                ref_count += 1;
            }

            for (&(lca, path_len), forest) in iter {
                match path_len {
                    ..=2 => unreachable!(),
                    3 => continue,
                    4 => (),
                    5.. => return false,
                }
                let mut count = 0;

                for Label(label) in forest.dfs_from(lca).labels() {
                    if label_set[label as usize] != generation {
                        return false;
                    }
                    count += 1;
                }
                if count != ref_count {
                    return false;
                }
            }
            true
        }

        let mut label_set_workspace: Vec<u8> =
            vec![0; self.forests[0].max_label_value() as usize + 1];
        let mut generation = 1;
        let target = {
            self.forests[0]
                .iterate_cherries()
                .find_map(|(_, label_l, label_r)| {
                    let lcas_and_paths_len =
                        self.get_lca_and_path_len_in_all_forests(label_l, label_r);
                    generation = if generation == u8::MAX {
                        label_set_workspace.fill(0);
                        1
                    } else {
                        generation + 1
                    };
                    if rr2_2_1_applicable(
                        self,
                        &lcas_and_paths_len,
                        &mut label_set_workspace,
                        generation,
                    ) {
                        return Some((label_l, label_r, lcas_and_paths_len));
                    }
                    None
                })
        };

        let Some((label_l, label_r, lcas_and_paths_len)) = target else {
            return false;
        };

        // apply reduction: remove the (at most 1) subtrees incident to the path between label_a and label_b in all forests
        self.cut_path_between_labels_in_all_forests(lcas_and_paths_len, label_l, label_r);
        true
    }

    /// Applies branching rule 2.2.2 on cherry (a,b)
    /// Assumes instance satisfies LSI
    fn apply_br2_2_2(mut self, label_a: Label, label_b: Label) -> BranchResult {
        let lcas_and_paths_len = self.get_lca_and_path_len_in_all_forests(label_a, label_b);

        // apply branches
        // branch [1]: remove edge incident to label_l in all forests
        let mut branch1_instance = self.clone();
        branch1_instance.cut_label_in_all_forests(label_a);

        // branch [2]: remove edge incident to label_r in all forests
        let mut branch2_instance = self.clone();
        branch2_instance.cut_label_in_all_forests(label_b);

        // branch [3]: remove all subtrees incident to the path between label_l and label_r in all forests
        self.cut_path_between_labels_in_all_forests(lcas_and_paths_len, label_a, label_b);

        Branch3([
            LsiInstance::from_sub_instance(branch1_instance),
            LsiInstance::from_sub_instance(branch2_instance),
            LsiInstance::from_sub_instance(self),
        ])
    }

    /// Cuts off the leaf with the specified label in all forests
    ///NOTE: Panics if any leaf already is a single vertex tree
    /// This can only happen if instance did not satisfy LSI
    #[inline]
    fn cut_label_in_all_forests(&mut self, label: Label) {
        for f in self.forests.iter_mut() {
            let label_idx = f.locate_label(label);
            f.cut_branch(label_idx);
        }
    }

    /// Cuts off all subtrees pendant to the path between two labels
    /// in all forests
    /// Takes in a vector with precomputed lcas and path lengths in each forest
    #[inline]
    fn cut_path_between_labels_in_all_forests(
        &mut self,
        lcas_and_paths_len: Vec<(Idx, usize)>,
        label_a: Label,
        label_b: Label,
    ) {
        debug_assert_eq!(lcas_and_paths_len.len(), self.forests.len());
        for (f, (lca, _)) in self.forests.iter_mut().zip(lcas_and_paths_len) {
            for leaf in [f.locate_label(label_a), f.locate_label(label_b)] {
                while let Some(parent) = f.get(leaf).parent
                    && parent != lca
                {
                    let [left, right] = f.get(parent).children_unchecked();
                    if left != leaf {
                        f.cut_branch(left);
                    } else {
                        debug_assert!(right != leaf);
                        f.cut_branch(right);
                    }
                }
            }
        }
    }

    /// Takes an instance satisfying the LSI and two labels
    /// and finds the length of the path between them in each forest
    /// and finds the lca in each forest
    ///
    /// # Panics if labels are not in the same component in some forest
    fn get_lca_and_path_len_in_all_forests(
        &self,
        label_a: Label,
        label_b: Label,
    ) -> Vec<(Idx, usize)> {
        let mut paths: Vec<(Idx, usize)> = Vec::with_capacity(self.forests.len());
        for f in self.forests.iter() {
            let a = f.locate_label(label_a);
            let b = f.locate_label(label_b);

            let Some(lca_and_path_len) = f.get_lca_and_path_len(a, b) else {
                unreachable!("expected an instance satisfying LSI");
            };

            paths.push(lca_and_path_len);
        }
        paths
    }

    /// Tests if f1 is a subgraph of all other forests.
    /// Note: if instance satisfies LSI then it is an AF
    fn f1_subgraph_of_all(&self) -> bool {
        for other in self.forests[1..].iter() {
            if !self.forests[0].is_subgraph_of(other) {
                return false;
            }
        }
        true
    }

    fn satisfies_lsi(&self) -> bool {
        // if the number of components does not match, it certainly does not satisfy LSI
        let f1 = &self.forests[0];
        let n_comps = f1.roots.len();

        if self.forests[1..].iter().any(|f| f.roots.len() != n_comps) {
            return false;
        }

        // store component index of each label in f1
        let mut comp_idx_of_label: Vec<usize> =
            vec![usize::MAX; self.forests[0].max_label_value() as usize + 1];
        for (comp_idx, &root) in f1.roots.iter().enumerate() {
            // denote the component idx of all labels in this component
            for Label(label) in f1.dfs_from(root).labels() {
                comp_idx_of_label[label as usize] = comp_idx;
            }
        }

        // check if component index of each label in f1 is the same
        for other in self.forests[1..].iter() {
            for &root in other.roots.iter() {
                let mut labels_other = other.dfs_from(root).labels();
                let comp_idx = comp_idx_of_label[labels_other
                    .next()
                    .expect("each component must have a label")
                    .0 as usize];
                if labels_other.any(|Label(label)| comp_idx_of_label[label as usize] != comp_idx) {
                    return false;
                }
            }
        }
        true
    }
}
