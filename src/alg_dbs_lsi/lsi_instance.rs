#[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
use crate::maf_instance::arena_vertex::Label;

use super::instance_lsi_ext::InstanceLsiDbsExt;
use crate::maf_instance::{
    arena_tree::ArenaTree, instance::Instance, performed_reduction::PerformedReduction,
};

#[derive(Debug, Clone)]
pub struct LsiInstance {
    pub main: Instance,
    pub lsi_clusters: Vec<Instance>,
}

impl LsiInstance {
    pub fn get_num_leaves(&self) -> usize {
        self.main.num_leaves
            + self
                .lsi_clusters
                .iter()
                .map(|inst| inst.num_leaves)
                .sum::<usize>()
    }

    #[inline]
    pub fn ord(&self) -> usize {
        let lsi_clusters_ord: usize = self.lsi_clusters.iter().map(|inst| inst.ord()).sum();
        self.main.ord() + lsi_clusters_ord
    }

    pub fn extract_af(&mut self) -> ArenaTree {
        let mut sol = self.main.extract_af();
        for lsi_cluster in self.lsi_clusters.iter_mut() {
            sol = sol.join_with(lsi_cluster.extract_af());
        }
        sol
    }

    pub fn f1_subgraph_of_all(&self) -> bool {
        self.main.f1_subgraph_of_all() && self.lsi_clusters.iter().all(|c| c.f1_subgraph_of_all())
    }

    /// Exhaustively performs reductions:
    /// - merges common subtrees between all forests
    /// - syncs single vertex trees between the forests and removes them
    ///
    /// returns the removed svts and the merged labels
    pub fn fully_r2_reduce_and_remove_svt(&mut self) -> (Vec<PerformedReduction>, usize) {
        let (mut performed_reductions, mut n_removed_svts) =
            self.main.fully_r2_reduce_and_remove_svt();

        for lsi_cluster in self.lsi_clusters.iter_mut() {
            let (new_performed_reductions, new_n_removed_svts) =
                lsi_cluster.fully_r2_reduce_and_remove_svt();
            performed_reductions.extend(new_performed_reductions);
            n_removed_svts += new_n_removed_svts;
        }

        self.lsi_clusters.retain(|f| f.num_leaves > 0);

        (performed_reductions, n_removed_svts)
    }

    pub fn from_sub_instance(instance: Instance) -> LsiInstance {
        LsiInstance {
            main: instance,
            lsi_clusters: vec![],
        }
    }

    /// Tests if the instance satisfies the LSI property.
    /// Extracts all components that do satisfy it into the lsi_clusters list
    /// to save later recalculation
    #[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
    pub fn satisfies_lsi_mut(&mut self) -> bool {
        use crate::common::generation_set::GenerationSet;
        use crate::maf_instance::arena_vertex::Idx;
        use crate::maf_instance::tree_traversal::TreeTraversal;

        let f1 = &self.main.forests[0];

        // iterate through all the roots of f1 that may not have LSI
        let mut is_fully_lsi = true;
        let mut newly_lsi_roots: Vec<Idx> = Vec::with_capacity(f1.roots.len());
        let mut label_set: GenerationSet = GenerationSet::new(f1.max_label_value() as usize + 1);
        'components: for &root in f1.roots.iter() {
            let mut labels = f1.dfs_from(root).labels().peekable();
            let first_label = *labels.peek().expect("at least one label in each component");
            label_set.advance();
            let mut f1_count = 0;
            for Label(label) in labels {
                label_set.insert(label as usize);
                f1_count += 1;
            }
            for f in self.main.forests.iter().skip(1) {
                let mut other_count = 0;
                let leaf_in_other = f.locate_label(first_label);
                let root_other = f.find_root_of(leaf_in_other);
                let labels_other = f.dfs_from(root_other).labels();
                for Label(label) in labels_other {
                    other_count += 1;
                    if !label_set.contains(label as usize) {
                        is_fully_lsi = false;
                        continue 'components;
                    }
                }
                if other_count != f1_count {
                    is_fully_lsi = false;
                    continue 'components;
                }
            }
            newly_lsi_roots.push(root);
        }

        // extract all lsi-clusters from the main forest
        let first_label_of_newly_lsi_roots: Vec<Label> = newly_lsi_roots
            .iter()
            .map(|root_f1| {
                self.main.forests[0]
                    .dfs_from(*root_f1)
                    .labels()
                    .next()
                    .unwrap()
            })
            .collect();

        for label in first_label_of_newly_lsi_roots
            .into_iter()
            .skip(if is_fully_lsi { 1 } else { 0 })
        {
            let cluster = self.main.extract_lsi_cluster_with_label(label);
            self.lsi_clusters.push(cluster);
        }
        is_fully_lsi
    }

    /// Tests if the instance satisfies the LSI property.
    #[cfg(feature = "LSI_DISABLE_CLUSTER_REDUCTION")]
    pub fn satisfies_lsi_mut(&mut self) -> bool {
        debug_assert!(self.lsi_clusters.is_empty());

        self.main.satisfies_lsi()
    }

    #[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
    pub fn find_clusters_w_dummy(
        &self,
    ) -> Option<(
        LsiInstance,
        LsiInstance,
        (LsiInstance, Label),
        (LsiInstance, Label),
    )> {
        use super::arena_tree_lsi_ext::ArenaTreeLsiExt;
        use crate::maf_instance::arena_vertex::Idx;
        use crate::maf_instance::arena_vertex::NodeData::*;
        use itertools::kmerge;
        use std::collections::HashMap;

        const MIN_CLUSTER_SIZE_BELOW: usize = 5;
        const MIN_CLUSTER_SIZE_ABOVE: usize = 5;

        let ord = self.main.ord();
        let n = self.main.forests[0].leaf_map.len();

        let max_cluster_size = n.checked_sub(MIN_CLUSTER_SIZE_ABOVE)?;
        let corrected_max_cluster_size = max_cluster_size.checked_sub(ord - 1)?;

        if n < MIN_CLUSTER_SIZE_BELOW + corrected_max_cluster_size {
            return None;
        }

        let goal = (n.checked_sub(ord - 1)? as f64 * 0.5).floor() as usize;
        // let goal = (n.checked_sub(ord - 1)? as f64 * 0.75).floor() as usize;

        // build map storing the label set of each subtree as
        // a sorted vec
        let mut clusters: Vec<HashMap<Vec<Label>, Idx>> = vec![];

        for f in self.main.forests.iter() {
            let mut label_sets: HashMap<Idx, Vec<Label>> = HashMap::new();
            for &root in f.roots.iter() {
                for (idx, node) in f.dfs_postorder(root) {
                    // NOTE: we skip roots here, because those are cases for LSI cluster reduction
                    // instead
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
            if cluster_size < MIN_CLUSTER_SIZE_BELOW || cluster_size > corrected_max_cluster_size {
                // this cluster cuts off too little nodes to be useful
                // clusters of size 3 are also useless, because
                // reduction rule 2.2.1 is always applicable on clusters of size 3
                // resulting in zero branches in the clusters
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
        for (forest, clusters_f) in self.main.forests.iter().zip(clusters.iter()) {
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

        let lsi_instance_above = LsiInstance::from_sub_instance(Instance {
            forests: cluster_above,
            num_leaves: n - cluster.len(),
        });
        let lsi_instance_above_w_dummy = LsiInstance::from_sub_instance(Instance {
            forests: cluster_above_w_dummy,
            num_leaves: n - cluster.len() + 1,
        });
        let lsi_instance_below = LsiInstance::from_sub_instance(Instance {
            forests: cluster_below,
            num_leaves: cluster.len(),
        });
        let lsi_instance_below_w_dummy = LsiInstance::from_sub_instance(Instance {
            forests: cluster_below_w_dummy,
            num_leaves: cluster.len() + 1,
        });

        Some((
            lsi_instance_above,
            lsi_instance_below,
            (lsi_instance_above_w_dummy, dummy_label_above),
            (lsi_instance_below_w_dummy, dummy_label_below),
        ))
    }

    #[cfg(feature = "assert_validity")]
    pub fn assert_validity(&self) {
        use crate::maf_instance::arena_vertex::Label;
        use std::collections::HashSet;

        self.main.assert_validity();

        for lsi_cluster in self.lsi_clusters.iter() {
            lsi_cluster.assert_validity();
        }

        // verify that the labels in lsi_clusters do not overlap
        // eachother or those in the main instance
        let mut seen_labels: HashSet<Label> =
            self.main.forests[0].leaf_map.keys().copied().collect();
        for lsi_cluster in self.lsi_clusters.iter() {
            for &label in lsi_cluster.forests[0].leaf_map.keys() {
                assert!(seen_labels.insert(label));
            }
        }
    }
}
