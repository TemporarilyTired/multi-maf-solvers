use itertools::Itertools;
use std::collections::HashMap;

use crate::maf_instance::{
    arena_vertex::{Idx, Label, NodeData::*},
    instance::Instance,
    tree_traversal::TreeTraversal,
};

pub trait AfDbsInstanceExt {
    fn calculate_lcas(&self) -> Vec<Vec<Vec<usize>>>;

    fn find_incompatible_triples(
        &self,
        labels: &[Label],
        lcas: &[Vec<Vec<usize>>],
        comp_groups: &[usize],
    ) -> Vec<u64>;

    fn find_non_crossing_paths_in_any(
        &self,
        labels: &[Label],
        lcas: &[Vec<Vec<usize>>],
        comp_groups: &[usize],
    ) -> Vec<u64>;

    fn calculate_comp_groups(&self) -> Vec<usize>;
}

impl AfDbsInstanceExt for Instance {
    /// calculate the 2d map of the lca of two nodes for each forest
    fn calculate_lcas(&self) -> Vec<Vec<Vec<usize>>> {
        // NOTE: this could be way more efficient; but even this takes negligible time
        // since it is only used once per instance/cluster
        debug_assert!(!self.forests.is_empty());

        let n_forests = self.forests.len();
        let n_nodes = self
            .forests
            .iter()
            .map(|f| f.arena.len())
            .max()
            .unwrap_or_default();
        let mut lcas: Vec<Vec<Vec<usize>>> =
            vec![vec![vec![usize::MAX; n_nodes]; n_nodes]; n_forests];
        for (lcas_f, f) in lcas.iter_mut().zip_eq(self.forests.iter()) {
            for &root in f.roots.iter() {
                for (node_idx, node) in f.dfs_from(root) {
                    let node_idx_u = node_idx as usize;
                    match node.data {
                        Leaf { .. } => lcas_f[node_idx_u][node_idx_u] = node_idx_u,
                        Internal { left, right } => {
                            let parent: u16 = node_idx;
                            let parent_u: usize = node_idx_u;
                            let left_indices: Vec<Idx> =
                                f.dfs_from(left).indices().chain([parent]).collect();
                            let right_indices: Vec<Idx> =
                                f.dfs_from(right).indices().chain([parent]).collect();
                            for left_idx in left_indices {
                                for &right_idx in right_indices.iter() {
                                    lcas_f[left_idx as usize][right_idx as usize] = parent_u;
                                    lcas_f[right_idx as usize][left_idx as usize] = parent_u;
                                }
                            }
                        }
                    }
                }
            }

            #[cfg(debug_assertions)]
            {
                for idx_a in f.iterate_all().indices() {
                    for idx_b in f.iterate_all().indices() {
                        let ab = lcas_f[idx_a as usize][idx_b as usize];
                        let ba = lcas_f[idx_b as usize][idx_a as usize];
                        debug_assert_eq!(ab, ba);
                        if ab == usize::MAX {
                            debug_assert_ne!(f.find_root_of(idx_a), f.find_root_of(idx_b));
                        } else {
                            debug_assert_eq!(f.find_root_of(idx_a), f.find_root_of(idx_b));
                            debug_assert!(f.ancestors_incl(idx_a).contains(&(ab as Idx)));
                            debug_assert!(f.ancestors_incl(idx_b).contains(&(ab as Idx)));
                        }
                    }
                }
            }
        }
        lcas
    }

    /// Calculate a comp-group identifier for each label s.t.
    /// the identifier is equal for labels a and b iff labels
    /// a and b are in the same component in every forest
    fn calculate_comp_groups(&self) -> Vec<usize> {
        let mut comp_groups: Vec<usize> =
            vec![usize::MAX; self.forests[0].max_label_value() as usize + 1];
        let mut comp_identifier: HashMap<Vec<Idx>, usize> = HashMap::new();
        for label in self.forests[0].iterate_all().labels() {
            let roots_of_label: Vec<Idx> = self
                .forests
                .iter()
                .map(|f| f.find_root_of(f.locate_label(label)))
                .collect();
            if let Some(&comp_group) = comp_identifier.get(&roots_of_label) {
                comp_groups[label.0 as usize] = comp_group;
            } else {
                let ident = comp_identifier.len();
                comp_identifier.insert(roots_of_label, ident);
                comp_groups[label.0 as usize] = ident;
            }
        }
        comp_groups
    }

    /// Find all pairs of pairs of labels that are disjoint in each input tree.
    /// Uses pre-calculated map of lcas of each pair of nodes in each tree
    /// Is the complement of 'find_crossing_paths_in_any', with the same complexity
    /// (when implemented the same way), but returns 50%-95% less pairs.
    /// WARN: result can be over 2BG of memory if trees have >300 labels
    fn find_non_crossing_paths_in_any(
        &self,
        labels: &[Label],
        lcas: &[Vec<Vec<usize>>],
        comp_groups: &[usize],
    ) -> Vec<u64> {
        let n_forests = self.forests.len();
        let max_label = labels.iter().max().unwrap_or(&Label(0)).0 as usize;
        let max_node_idx = self
            .forests
            .iter()
            .map(|f| f.arena.len())
            .max()
            .unwrap_or_default();

        // n^2 is a conservative guess on #non-crossing pairs
        let mut res: Vec<_> = Vec::with_capacity(labels.len() * labels.len());

        // NOTE: pre-calculate lca of each pair of labels in each tree
        let mut label_lcas: Vec<Vec<Vec<usize>>> =
            vec![vec![vec![usize::MAX; max_label + 1]; max_label + 1]; n_forests];
        for ((label_lcas_f, lcas_f), f) in label_lcas
            .iter_mut()
            .zip_eq(lcas.iter())
            .zip_eq(&self.forests)
        {
            for &a in labels.iter() {
                let leaf_a = f.locate_label(a) as usize;
                for &b in labels.iter() {
                    let leaf_b = f.locate_label(b) as usize;
                    label_lcas_f[a.0 as usize][b.0 as usize] = lcas_f[leaf_a][leaf_b];
                }
            }
        }

        // NOTE: pre-calculate depth of each node in each tree
        let mut depth: Vec<Vec<usize>> = vec![vec![usize::MAX; max_node_idx + 1]; n_forests];
        for (depth_f, f) in depth.iter_mut().zip_eq(&self.forests) {
            for &root in f.roots.iter() {
                for (node_idx, node) in f.dfs_from(root) {
                    depth_f[node_idx as usize] = match node.parent {
                        None => 0,
                        Some(parent) => {
                            debug_assert_ne!(depth_f[parent as usize], usize::MAX);
                            depth_f[parent as usize] + 1
                        }
                    };
                }
            }
        }

        // iterate over all pairs of labels (a,b),(c,d) where (a<b) (c<d) and (a<c)
        // and label a lies in the same component as label b in each forest (same with c and d)
        //
        // We can skip instances where a and b (or c and d) are separated in SOME forest
        // because we know that in any AF they must also be separated due to this.
        for &Label(d) in labels.iter() {
            let comp_group_d = comp_groups[d as usize];
            for &Label(c) in labels.iter() {
                if c >= d {
                    break;
                }
                let comp_group_c = comp_groups[c as usize];
                if comp_group_c != comp_group_d {
                    continue;
                }
                for &Label(b) in labels.iter() {
                    let comp_group_b = comp_groups[b as usize];
                    'label_a_loop: for &Label(a) in labels.iter() {
                        if a >= b || a >= c {
                            break;
                        }
                        let comp_group_a = comp_groups[a as usize];
                        if comp_group_a != comp_group_b {
                            continue;
                        }
                        for ((label_lcas_f, depth_f), f) in label_lcas
                            .iter()
                            .zip_eq(depth.iter())
                            .zip_eq(self.forests.iter())
                        {
                            // if pair a,b lies in a different component than pair c,d
                            // there is no need to check if their paths overlap
                            if f.find_root_of(f.locate_label(Label(a)))
                                != f.find_root_of(f.locate_label(Label(c)))
                            {
                                continue;
                            }

                            if f.paths_intersect(
                                a as usize,
                                b as usize,
                                c as usize,
                                d as usize,
                                label_lcas_f,
                                depth_f,
                            ) {
                                continue 'label_a_loop;
                            }
                        }
                        res.push(quad_key(Label(a), Label(b), Label(c), Label(d)));
                    }
                }
            }
        }
        res
    }

    /// Find all triples of labels that have a different embedding in at least one
    /// pair of input trees.
    /// Uses pre-calculated map of lcas of each pair of nodes in each tree
    fn find_incompatible_triples(
        &self,
        labels: &[Label],
        lcas: &[Vec<Vec<usize>>],
        comp_groups: &[usize],
    ) -> Vec<u64> {
        let n_forests = self.forests.len();
        let max_label = labels.iter().max().unwrap_or(&Label(0)).0 as usize;
        // guess that around n^2 incompatible triples will be found (at most n^3, but that is a lot of memory)
        let mut res: Vec<_> = Vec::with_capacity(labels.len() * labels.len());

        let mut label_lcas: Vec<Vec<Vec<usize>>> =
            vec![vec![vec![usize::MAX; max_label + 1]; max_label + 1]; n_forests];

        for (label_lcas_f, (lcas_f, f)) in label_lcas
            .iter_mut()
            .zip_eq(lcas.iter().zip_eq(&self.forests))
        {
            for &a in labels.iter() {
                let leaf_a = f.locate_label(a) as usize;
                for &b in labels.iter() {
                    let leaf_b = f.locate_label(b) as usize;
                    label_lcas_f[a.0 as usize][b.0 as usize] = lcas_f[leaf_a][leaf_b];
                }
            }
        }

        // iterate over all sorted triples of labels (a,b,c) that
        // are in the same component in every forest
        for &Label(c) in labels.iter() {
            let comp_group_c = comp_groups[c as usize];
            for &Label(b) in labels.iter() {
                if b >= c {
                    break;
                }
                let comp_group_b = comp_groups[b as usize];
                if comp_group_b != comp_group_c {
                    continue;
                }
                'label_a_loop: for &Label(a) in labels.iter() {
                    if a >= b {
                        break;
                    }
                    let comp_group_a = comp_groups[a as usize];
                    if comp_group_a != comp_group_b {
                        continue;
                    }

                    let lca_ab = label_lcas[0][a as usize][b as usize];
                    let lca_bc = label_lcas[0][b as usize][c as usize];
                    let lca_ac = label_lcas[0][a as usize][c as usize];

                    for label_lcas_other in label_lcas.iter().skip(1) {
                        let lca_ab_other = label_lcas_other[a as usize][b as usize];
                        let lca_bc_other = label_lcas_other[b as usize][c as usize];
                        let lca_ac_other = label_lcas_other[a as usize][c as usize];

                        if ((lca_ab == lca_ac) != (lca_ab_other == lca_ac_other))
                            || ((lca_ac == lca_bc) != (lca_ac_other == lca_bc_other))
                            || ((lca_ab == lca_bc) != (lca_ab_other == lca_bc_other))
                        {
                            // abc is an incompatible triple (in sorted order)
                            res.push(triple_key(Label(a), Label(b), Label(c)));
                            continue 'label_a_loop;
                        }
                    }
                }
            }
        }

        res
    }
}

#[inline(always)]
pub fn triple_key(a: Label, b: Label, c: Label) -> u64 {
    ((a.0 as u64) << 32) | ((b.0 as u64) << 16) | (c.0 as u64)
}

#[inline(always)]
pub fn quad_key(a: Label, b: Label, c: Label, d: Label) -> u64 {
    ((a.0 as u64) << 48) | ((b.0 as u64) << 32) | ((c.0 as u64) << 16) | (d.0 as u64)
}
