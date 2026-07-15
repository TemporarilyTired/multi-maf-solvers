use itertools::kmerge;
use std::collections::{HashMap, HashSet};

use super::arena_tree_ext::ReductionsArenaTreeExt;
use crate::{
    alg_dbs_lsi::{arena_tree_lsi_ext::ArenaTreeLsiExt, instance_lsi_ext::InstanceLsiDbsExt},
    maf_instance::{
        arena_tree::ArenaTree,
        arena_vertex::{Idx, Label, NodeData::*, Status::*},
        instance::Instance,
        performed_reduction::PerformedReduction::{self, *},
        tree_traversal::TreeTraversal,
    },
};

pub trait ReductionsInstanceExt {
    fn fully_reduce(&mut self) -> Vec<PerformedReduction>;
    fn rr_2_2_1_reduce(&mut self) -> bool;
    fn reduce_svt_and_merge(&mut self) -> (Vec<PerformedReduction>, usize);

    fn split_into_clusters(self) -> Vec<Instance>;
    fn copy_cluster(&self, cluster_labels: &HashSet<Label>) -> Instance;

    #[allow(clippy::type_complexity)]
    fn find_clusters_w_dummy(
        &self,
    ) -> Option<(Instance, Instance, (Instance, Label), (Instance, Label))>;
}

impl ReductionsInstanceExt for Instance {
    /// Exhaustively performs reductions:
    /// - merge common subtrees between all forests
    /// - sync single vertex trees between the forests
    /// - perform reduction rule 2.2.1
    fn fully_reduce(&mut self) -> Vec<PerformedReduction> {
        let (mut performed_reductions, _) = self.reduce_svt_and_merge();

        while self.rr_2_2_1_reduce() || self.fully_r1_reduce() {
            let (new_reductions, _) = self.reduce_svt_and_merge();
            performed_reductions.extend(new_reductions);
        }
        performed_reductions
    }

    /// Exhaustively performs reductions:
    /// - merges common subtrees between all forests
    /// - syncs single vertex trees between the forests and removes them
    ///
    /// returns the removed svts and the merged labels
    fn reduce_svt_and_merge(&mut self) -> (Vec<PerformedReduction>, usize) {
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

    /// Find cherries s.t. in each forest
    /// - there are 0 or 1 pendant subtrees between it, and
    /// - the label set of these pendant subtrees is equal
    ///
    /// These pendant subtrees are cut
    /// NOTE: does NOT perform this reduction exhaustively
    ///
    /// Returns true if at least one reduction was applied
    fn rr_2_2_1_reduce(&mut self) -> bool {
        let mut applied_a_reduction = false;

        let cherries: Vec<(Label, Label)> = self
            .iterate_all_cherries()
            .map(|(_, a, b)| (a, b))
            .collect();

        'cherries: for (a, b) in cherries {
            // check if labels still exist (after previous reduction iterations)
            let f1 = &self.forests[0];
            if f1.try_locate_label(a).is_none() || f1.try_locate_label(b).is_none() {
                continue 'cherries;
            }

            // check if label is a cherry or uncle-nephew in every forest with the same label set
            let mut label_set_opt: Option<HashSet<Label>> = None;
            for forest in &self.forests {
                let leaf_a = forest.locate_label(a);
                let leaf_b = forest.locate_label(b);

                let Some(parent_a) = forest.get(leaf_a).parent else {
                    continue 'cherries;
                };
                let Some(parent_b) = forest.get(leaf_b).parent else {
                    continue 'cherries;
                };
                let sibling_a = forest.find_sibling(leaf_a).unwrap();
                let sibling_b = forest.find_sibling(leaf_b).unwrap();

                let pendant_subtree = if parent_a == parent_b {
                    // this is a cherry, no further checks needed
                    continue;
                } else if sibling_a == parent_b {
                    // a is the uncle of b: ((b,_), a)
                    sibling_b
                } else if sibling_b == parent_a {
                    // b is the uncle of a: ((a,_), b)
                    sibling_a
                } else {
                    // there is more than 1 pendant subtree between a and b
                    // (or they are in different components)
                    continue 'cherries;
                };

                let new_label_set = forest.dfs_from(pendant_subtree).labels().collect();
                if let Some(label_set) = &label_set_opt {
                    if new_label_set != *label_set {
                        continue 'cherries;
                    }
                } else {
                    label_set_opt = Some(new_label_set);
                }
            }

            // perform the cut in each forest (if not a cherry already)
            for forest in self.forests.iter_mut() {
                let leaf_a = forest.locate_label(a);
                let leaf_b = forest.locate_label(b);

                let parent_a = forest.get(leaf_a).parent.unwrap();
                let parent_b = forest.get(leaf_b).parent.unwrap();
                let sibling_a = forest.find_sibling(leaf_a).unwrap();
                let sibling_b = forest.find_sibling(leaf_b).unwrap();

                let pendant_subtree = if parent_a == parent_b {
                    // this is a cherry, no further checks needed
                    continue;
                } else if sibling_a == parent_b {
                    // a is the uncle of b: ((b,_), a)
                    sibling_b
                } else if sibling_b == parent_a {
                    // b is the uncle of a: ((a,_), b)
                    sibling_a
                } else {
                    // there is more than 1 pendant subtree between a and b
                    // (or they are in different components)
                    unreachable!("this cut should not have been performed");
                };

                forest.cut_branch(pendant_subtree);
            }
            applied_a_reduction = true;
        }
        applied_a_reduction
    }

    fn split_into_clusters(self) -> Vec<Instance> {
        let mut expanded_labels = HashSet::<Label>::new();
        let mut clusters: Vec<HashSet<Label>> = vec![];

        // loop through all labels, skipping the ones already processed
        for label in self.forests[0]
            .leaf_map
            .keys()
            .cloned()
            .collect::<Vec<Label>>()
        {
            if expanded_labels.contains(&label) {
                continue;
            }

            // collect all labels that are connected via some route
            // (possibly alternating between forests)
            let mut current_cluster = HashSet::<Label>::new();

            let mut to_expand = vec![label];
            let mut to_expand_set = HashSet::<Label>::from_iter(to_expand.clone());
            while let Some(neighbor) = to_expand.pop() {
                to_expand_set.remove(&neighbor);
                debug_assert_eq!(
                    HashSet::<Label>::from_iter(to_expand.clone()),
                    to_expand_set
                );

                if !expanded_labels.insert(neighbor) {
                    continue;
                }
                current_cluster.insert(neighbor);

                for f in self.forests.iter() {
                    let leaf = f.locate_label(neighbor);
                    let comp = f.find_root_of(leaf);
                    for neighbor2 in f.dfs_from(comp).labels() {
                        if !to_expand.contains(&neighbor2) && to_expand_set.insert(neighbor2) {
                            to_expand.push(neighbor2);
                        }
                    }
                }
            }
            clusters.push(current_cluster);
        }

        clusters
            .into_iter()
            .map(|cluster| self.copy_cluster(&cluster))
            .collect()
    }

    fn copy_cluster(&self, cluster_labels: &HashSet<Label>) -> Instance {
        Instance {
            num_leaves: cluster_labels.len(),
            forests: self
                .forests
                .iter()
                .map(|f| f.copy_cluster(cluster_labels))
                .collect(),
        }
    }

    fn find_clusters_w_dummy(
        &self,
    ) -> Option<(Instance, Instance, (Instance, Label), (Instance, Label))> {
        // use super::arena_tree_lsi_ext::ArenaTreeLsiExt;
        const MIN_CLUSTER_SIZE_BELOW: usize = 3;
        const MIN_CLUSTER_SIZE_ABOVE: usize = 3;

        let n = self.forests[0].leaf_map.len();
        let max_cluster_size = n.checked_sub(MIN_CLUSTER_SIZE_ABOVE)?;

        if n < MIN_CLUSTER_SIZE_BELOW + MIN_CLUSTER_SIZE_ABOVE {
            return None;
        }

        let goal = n / 2;

        // build map storing the label set of each subtree as
        // a sorted vec
        let mut clusters: Vec<HashMap<Vec<Label>, Idx>> = vec![];

        for f in self.forests.iter() {
            let mut label_sets: HashMap<Idx, Vec<Label>> = HashMap::new();
            for &root in f.roots.iter() {
                for (idx, node) in f.dfs_postorder(root) {
                    // NOTE: we skip roots here, because they are handled without dummies
                    if node.parent.is_none() {
                        continue;
                    }
                    let label_set = match node.data {
                        Internal { left, right } => {
                            kmerge([&label_sets[&left], &label_sets[&right]])
                                .copied()
                                .collect()
                        }
                        Leaf { label } => vec![label],
                    };

                    label_sets.insert(idx, label_set);
                }
            }

            // invert the map
            let node_of_label_set = label_sets.drain().map(|(k, v)| (v, k)).collect();
            clusters.push(node_of_label_set);
        }

        let Some((f1_clusters, other_clusters)) = clusters.split_first() else {
            unreachable!("can only try to find clusters on two or more trees");
        };

        let mut best_cluster: Option<&Vec<Label>> = None;

        for cluster in f1_clusters.keys() {
            let cluster_size = cluster.len();
            if cluster_size < MIN_CLUSTER_SIZE_BELOW || cluster_size > max_cluster_size {
                // this cluster cuts off too little nodes to be useful
                continue;
            }

            if best_cluster
                .is_some_and(|best| best.len().abs_diff(goal) < cluster_size.abs_diff(goal))
            {
                // this cluster is farther from the goal size than the current best cluster
                continue;
            }

            if other_clusters
                .iter()
                .any(|other| !other.contains_key(cluster))
            {
                // this is not a common cluster
                continue;
            }

            // this is a common cluster, better balanced than the current best
            best_cluster = Some(cluster);
        }

        let cluster = best_cluster?;

        // let dummy_label_above = Label(self.main.forests[0].max_label_value() + 1);
        let dummy_label_above = *cluster.iter().min().expect("at least one label in cluster");
        let dummy_label_below = Label(0);

        let mut cluster_above: Vec<ArenaTree> = vec![];
        let mut cluster_above_w_dummy: Vec<ArenaTree> = vec![];
        let mut cluster_below: Vec<ArenaTree> = vec![];
        let mut cluster_below_w_dummy: Vec<ArenaTree> = vec![];
        for (forest, clusters_f) in self.forests.iter().zip(clusters.iter()) {
            let cluster_subtree = clusters_f[cluster];
            let (forest_above, tree_below, forest_above_w_dummy, tree_below_w_dummy) =
                forest.clone().split_at_subtree_with_dummy(
                    cluster_subtree,
                    (dummy_label_above, dummy_label_below),
                );

            cluster_above.push(forest_above);
            cluster_above_w_dummy.push(forest_above_w_dummy);
            cluster_below.push(tree_below);
            cluster_below_w_dummy.push(tree_below_w_dummy);
        }

        let lsi_instance_above = Instance {
            forests: cluster_above,
            num_leaves: n - cluster.len(),
        };
        let lsi_instance_above_w_dummy = Instance {
            forests: cluster_above_w_dummy,
            num_leaves: n - cluster.len() + 1,
        };
        let lsi_instance_below = Instance {
            forests: cluster_below,
            num_leaves: cluster.len(),
        };
        let lsi_instance_below_w_dummy = Instance {
            forests: cluster_below_w_dummy,
            num_leaves: cluster.len() + 1,
        };

        Some((
            lsi_instance_above,
            lsi_instance_below,
            (lsi_instance_above_w_dummy, dummy_label_above),
            (lsi_instance_below_w_dummy, dummy_label_below),
        ))
    }
}
