use std::collections::{HashMap, HashSet};

use crate::common::arena_tree_canonical_id::{CanonId, Canonicalizer};
use crate::common::validity::assert_validity;
use crate::maf_instance::arena_tree::ArenaTree;
use crate::maf_instance::arena_vertex::{Idx, Label, NodeData::*, Status::*};
use crate::maf_instance::tree_traversal::TreeTraversal;

use super::r1_workspace::R1Workspace;

pub trait ArenaTreeLsiExt {
    // -- functions for LSI algorithm --
    #[allow(deprecated)]
    fn find_br1_target_with(&self, other: &ArenaTree) -> Option<Idx>;

    /// exhaustively r1-reduce a pair of forests
    /// returns a tuple denoting for each forests if it changed
    fn r1_reduce_pair(a: &mut ArenaTree, b: &mut ArenaTree) -> (bool, bool);

    /// exhaustively r1-reduce a pair of forests
    /// returns a tuple denoting for each forests if it changed
    /// uses a workspace to avoid (re)allocation
    fn r1_reduce_pair_w_workspace(
        a: &mut ArenaTree,
        b: &mut ArenaTree,
        workspace: &mut R1Workspace,
    ) -> (bool, bool);

    fn find_r1_target_with(&self, other: &ArenaTree) -> Option<Idx>;

    fn find_r1_target_with_w_workspace(
        &self,
        other: &ArenaTree,
        ws: &mut R1Workspace,
    ) -> Option<Idx>;

    fn get_lca_and_path_len(&self, leaf_a: Idx, leaf_b: Idx) -> Option<(Idx, usize)>;

    /// Tests if self and other satisfy the Label-Set-Isomorphism
    /// (i.e., have the same set of sets of labels)
    fn has_same_label_sets(&self, other: &ArenaTree) -> bool;

    /// Tests if self is a subgraph of other.
    /// IMPORTANT: assumes that self and other satisfy LSI
    fn is_subgraph_of(&self, other: &ArenaTree) -> bool;

    #[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
    fn extract_lsi_comp(&mut self, lsi_comp: Idx) -> ArenaTree;

    fn split_at_subtree_with_dummy(
        self,
        subtree: Idx,
        dummy_labels: (Label, Label),
    ) -> (ArenaTree, ArenaTree, ArenaTree, ArenaTree);
}

impl ArenaTreeLsiExt for ArenaTree {
    // -- functions for LSI algorithm --
    fn find_br1_target_with(&self, other: &ArenaTree) -> Option<Idx> {
        for &root in self.roots.iter() {
            for (v, left, right) in self.dfs_from(root).internals() {
                if v == root {
                    continue;
                }
                let label_set_left: HashSet<Label> = self.dfs_from(left).labels().collect();
                let label_set_right: HashSet<Label> = self.dfs_from(right).labels().collect();

                // try both options for the children
                'try_children: for (label_set_a, label_set_b) in [
                    (&label_set_left, &label_set_right),
                    (&label_set_right, &label_set_left),
                ] {
                    for &label_a in label_set_a {
                        let leaf_idx_in_other = other.locate_label(label_a);
                        let root_other = other.find_root_of(leaf_idx_in_other);
                        if other
                            .dfs_from(root_other)
                            .labels()
                            .any(|label_other| label_set_b.contains(&label_other))
                        {
                            // some compoment has labels in subtree idx and in its sibling subtree
                            continue 'try_children;
                        }
                    }
                    return Some(v);
                }
            }
        }
        None
    }

    /// exhaustively r1-reduce a pair of forests
    /// returns a tuple denoting for each forests if it changed
    fn r1_reduce_pair(a: &mut ArenaTree, b: &mut ArenaTree) -> (bool, bool) {
        // keep reducing until it is not applicable in either direction
        let mut a_changed = false;
        let mut b_changed = false;
        loop {
            if let Some(cut) = a.find_r1_target_with(b) {
                a_changed = true;
                a.cut_branch(cut);
                // a.cut_branch(cut, false);
            } else if let Some(cut) = b.find_r1_target_with(a) {
                b_changed = true;
                b.cut_branch(cut);
                // b.cut_branch(cut, false);
            } else {
                break;
            }
        }
        (a_changed, b_changed)
    }

    /// exhaustively r1-reduce a pair of forests
    /// returns a tuple denoting for each forests if it changed
    /// uses a workspace to avoid (re)allocation
    fn r1_reduce_pair_w_workspace(
        a: &mut ArenaTree,
        b: &mut ArenaTree,
        workspace: &mut R1Workspace,
    ) -> (bool, bool) {
        // keep reducing until it is not applicable in either direction
        let mut a_changed = false;
        let mut b_changed = false;
        loop {
            if let Some(cut) = a.find_r1_target_with_w_workspace(b, workspace) {
                a_changed = true;
                a.cut_branch(cut);
                // a.cut_branch(cut, false);
            } else if let Some(cut) = b.find_r1_target_with_w_workspace(a, workspace) {
                b_changed = true;
                b.cut_branch(cut);
                // b.cut_branch(cut, false);
            } else {
                break;
            }
        }
        (a_changed, b_changed)
    }

    fn find_r1_target_with(&self, other: &ArenaTree) -> Option<Idx> {
        let mut root_other_of_label: Vec<Idx> = vec![0; other.max_label_value() as usize + 1];
        for &root_other in other.roots.iter() {
            for Label(label) in other.dfs_from(root_other).labels() {
                root_other_of_label[label as usize] = root_other;
            }
        }

        let max_root_other = *other.roots.iter().max().unwrap_or(&0) as usize;
        let mut component_checked: Vec<bool> = vec![false; max_root_other + 1];
        for &root in self.roots.iter() {
            let labels_component: HashSet<Label> = self.dfs_from(root).labels().collect();
            'vertices_in_comp: for v in self.dfs_from(root).indices() {
                if v == root {
                    continue;
                }
                let labels_subtree: HashSet<Label> = self.dfs_from(v).labels().collect();
                debug_assert!(labels_subtree.is_subset(&labels_component));

                component_checked.fill(false);
                for &Label(label) in &labels_subtree {
                    let root_other = root_other_of_label[label as usize];
                    if component_checked[root_other as usize] {
                        continue;
                    }
                    component_checked[root_other as usize] = true;
                    for label_other in other.dfs_from(root_other).labels() {
                        if labels_component.contains(&label_other)
                            && !labels_subtree.contains(&label_other)
                        {
                            continue 'vertices_in_comp;
                        }
                    }
                }
                // at this point all components in other overlapping with this subtree
                // have all their labels either not in this component or inside of this subtree
                debug_assert!(!self.has_same_label_sets(other));
                return Some(v);
            }
        }

        None
    }

    fn find_r1_target_with_w_workspace(
        &self,
        other: &ArenaTree,
        ws: &mut R1Workspace,
    ) -> Option<Idx> {
        for &root_other in other.roots.iter() {
            for label in other.dfs_from(root_other).labels() {
                ws.label_to_root[label.0 as usize] = root_other;
            }
        }

        for &root in self.roots.iter() {
            ws.component_labels.advance();
            for label in self.dfs_from(root).labels() {
                ws.component_labels.insert(label.0 as usize);
            }
            'vertices_in_comp: for v in self.dfs_from(root).indices() {
                if v == root {
                    continue;
                }
                ws.subtree_labels.advance();
                for label in self.dfs_from(v).labels() {
                    ws.subtree_labels.insert(label.0 as usize);
                }

                ws.checked_roots.advance();
                for label in self.dfs_from(v).labels() {
                    let root_other = ws.label_to_root[label.0 as usize];
                    debug_assert!((root_other as usize) < other.arena.len());
                    if ws.checked_roots.contains(root_other as usize) {
                        continue;
                    }
                    ws.checked_roots.insert(root_other as usize);

                    for Label(label_other) in other.dfs_from(root_other).labels() {
                        if ws.component_labels.contains(label_other as usize)
                            && !ws.subtree_labels.contains(label_other as usize)
                        {
                            continue 'vertices_in_comp;
                        }
                    }
                }
                // at this point all components in other overlapping with this subtree
                // have all their labels either not in this component or inside of this subtree
                debug_assert!(self.find_r1_target_with(other).is_some());
                debug_assert!(!self.has_same_label_sets(other));
                return Some(v);
            }
        }

        debug_assert!(self.find_r1_target_with(other).is_none());
        None
    }

    fn get_lca_and_path_len(&self, leaf_a: Idx, leaf_b: Idx) -> Option<(Idx, usize)> {
        let path_a = self.path_to_root(leaf_a);
        let path_b = self.path_to_root(leaf_b);

        let mut path_len: usize = path_a.len() + path_b.len();

        let path_a_rev = path_a.into_iter().rev();
        let path_b_rev = path_b.into_iter().rev();

        let mut lca = None;
        for (a, b) in path_a_rev.zip(path_b_rev) {
            if a != b {
                break;
            }
            path_len -= 2;
            lca = Some(a);
        }
        lca.map(|lca| (lca, path_len + 1))
    }

    /// Tests if self and other satisfy the Label-Set-Isomorphism
    /// (i.e., have the same set of sets of labels)
    fn has_same_label_sets(&self, other: &ArenaTree) -> bool {
        for &root in self.roots.iter() {
            let mut leaves = self.dfs_from(root).leaves().peekable();
            let (_, first_label) = leaves
                .peek()
                .expect("Expeted all components to contain at least one leaf");
            let first_leaf_other = other.locate_label(*first_label);
            let root_other = other.find_root_of(first_leaf_other);
            let leaves_other = other.dfs_from(root_other).leaves();

            let mut leaf_set: HashSet<Label> = HashSet::new();
            for (_, label) in leaves {
                leaf_set.insert(label);
            }
            for (_, label) in leaves_other {
                if !leaf_set.remove(&label) {
                    return false;
                }
            }
            if !leaf_set.is_empty() {
                return false;
            }
        }
        true
    }

    /// Tests if self is a subgraph of other.
    /// WARN: assumes that self and other satisfy LSI
    fn is_subgraph_of(&self, other: &ArenaTree) -> bool {
        let mut canon = Canonicalizer::new();
        let ids: Vec<CanonId> = self.canonical_ids_with(&mut canon);
        let ids_other: Vec<CanonId> = other.canonical_ids_with(&mut canon);

        let mut root_ids: Vec<CanonId> = self.roots.iter().map(|&idx| ids[idx as usize]).collect();
        // let mut root_ids: Vec<CanonId> = self.roots.keys().map(|&idx| ids[idx as usize]).collect();
        let mut root_ids_other: Vec<CanonId> = other
            .roots
            // .keys()
            .iter()
            .map(|&idx| ids_other[idx as usize])
            .collect();
        root_ids.sort_unstable();
        root_ids_other.sort_unstable();

        root_ids == root_ids_other
    }

    #[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
    fn extract_lsi_comp(&mut self, lsi_comp: Idx) -> ArenaTree {
        assert_validity!(self);

        let indices_of_subtree: HashSet<Idx> = self.dfs_from(lsi_comp).indices().collect();

        let mut translation_above: HashMap<Idx, Idx> = HashMap::new();
        let mut translation_below: HashMap<Idx, Idx> = HashMap::new();

        let mut above = vec![];
        let mut below = vec![];
        for (idx, node) in self.arena.iter().cloned().enumerate() {
            if node.status != Present {
                continue;
            }
            if indices_of_subtree.contains(&(idx as Idx)) {
                let new_idx = below.len() as Idx;
                below.push(node);
                translation_below.insert(idx as Idx, new_idx);
            } else {
                let new_idx = above.len() as Idx;
                above.push(node);
                translation_above.insert(idx as Idx, new_idx);
            }
        }

        let roots_above: HashSet<Idx> = above
            .iter()
            .enumerate()
            .filter_map(|(idx, node)| node.parent.is_none().then_some(idx as Idx))
            .collect();
        let roots_below: HashSet<Idx> = below
            .iter()
            .enumerate()
            .filter_map(|(idx, node)| node.parent.is_none().then_some(idx as Idx))
            .collect();

        let leaf_map_above: HashMap<Label, Idx> = above
            .iter()
            .enumerate()
            .filter_map(|(idx, node)| node.label().map(|label| (label, idx as Idx)))
            .collect();
        let leaf_map_below: HashMap<Label, Idx> = below
            .iter()
            .enumerate()
            .filter_map(|(idx, node)| node.label().map(|label| (label, idx as Idx)))
            .collect();

        for node in above.iter_mut() {
            node.parent = node.parent.map(|p| translation_above[&p]);
            node.data = match node.data.clone() {
                Internal { left, right } => Internal {
                    left: translation_above[&left],
                    right: translation_above[&right],
                },
                leaf => leaf,
            };
        }
        for node in below.iter_mut() {
            node.parent = node.parent.map(|p| translation_below[&p]);
            node.data = match node.data.clone() {
                Internal { left, right } => Internal {
                    left: translation_below[&left],
                    right: translation_below[&right],
                },
                leaf => leaf,
            };
        }

        let (c_above, c_below) = (
            ArenaTree {
                arena: above,
                roots: roots_above,
                leaf_map: leaf_map_above,
            },
            ArenaTree {
                arena: below,
                roots: roots_below,
                leaf_map: leaf_map_below,
            },
        );

        assert_validity!(c_above);
        assert_validity!(c_below);

        *self = c_above;
        c_below
    }

    fn split_at_subtree_with_dummy(
        mut self,
        subtree: Idx,
        (dummy_label_above, dummy_label_below): (Label, Label),
    ) -> (ArenaTree, ArenaTree, ArenaTree, ArenaTree) {
        assert_validity!(self);

        // store the sibling of the subtree, to later insert the dummy leaf as a sibling again
        let sibling_of_cut = self.find_sibling(subtree);

        if self.get(subtree).parent.is_some() {
            self.cut_branch(subtree);
        }

        let indices_of_subtree: HashSet<Idx> = self.dfs_from(subtree).indices().collect();

        let mut translation_above: HashMap<Idx, Idx> = HashMap::new();
        let mut translation_below: HashMap<Idx, Idx> = HashMap::new();

        let mut above = vec![];
        let mut below = vec![];
        for (idx, node) in self.arena.iter().cloned().enumerate() {
            if node.status != Present {
                continue;
            }
            if indices_of_subtree.contains(&(idx as Idx)) {
                let new_idx = below.len() as Idx;
                below.push(node);
                translation_below.insert(idx as Idx, new_idx);
            } else {
                let new_idx = above.len() as Idx;
                above.push(node);
                translation_above.insert(idx as Idx, new_idx);
            }
        }

        let roots_above: HashSet<Idx> = above
            .iter()
            .enumerate()
            .filter_map(|(idx, node)| node.parent.is_none().then_some(idx as Idx))
            .collect();
        let roots_below: HashSet<Idx> = below
            .iter()
            .enumerate()
            .filter_map(|(idx, node)| node.parent.is_none().then_some(idx as Idx))
            .collect();

        let leaf_map_above: HashMap<Label, Idx> = above
            .iter()
            .enumerate()
            .filter_map(|(idx, node)| node.label().map(|label| (label, idx as Idx)))
            .collect();
        let leaf_map_below: HashMap<Label, Idx> = below
            .iter()
            .enumerate()
            .filter_map(|(idx, node)| node.label().map(|label| (label, idx as Idx)))
            .collect();

        for node in above.iter_mut() {
            node.parent = node.parent.map(|p| translation_above[&p]);
            node.data = match node.data.clone() {
                Internal { left, right } => Internal {
                    left: translation_above[&left],
                    right: translation_above[&right],
                },
                leaf => leaf,
            };
        }
        for node in below.iter_mut() {
            node.parent = node.parent.map(|p| translation_below[&p]);
            node.data = match node.data.clone() {
                Internal { left, right } => Internal {
                    left: translation_below[&left],
                    right: translation_below[&right],
                },
                leaf => leaf,
            };
        }

        let (c_above, c_below) = (
            ArenaTree {
                arena: above,
                roots: roots_above,
                leaf_map: leaf_map_above,
            },
            ArenaTree {
                arena: below,
                roots: roots_below,
                leaf_map: leaf_map_below,
            },
        );

        // add the dummy leaves:
        // - for instance above:
        //  subdivide the edge above the sibling and give it dummy as other child

        // - for instance below:
        //  add a parent to the root (i.e., subtree) and give it dummy as other child
        let mut c_above_dummy = c_above.clone();
        c_above_dummy.add_dummy_leaf_as_sibling_of(
            dummy_label_above,
            sibling_of_cut.map(|sibling| translation_above[&sibling]),
        );
        let mut c_below_dummy = c_below.clone();
        c_below_dummy
            .add_dummy_leaf_as_sibling_of(dummy_label_below, Some(translation_below[&subtree]));

        assert_validity!(c_above);
        assert_validity!(c_below);
        assert_validity!(c_above_dummy);
        assert_validity!(c_below_dummy);

        (c_above, c_below, c_above_dummy, c_below_dummy)
    }
}
