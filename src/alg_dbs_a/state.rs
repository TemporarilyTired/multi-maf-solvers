use std::collections::HashSet;

use crate::maf_instance::{
    arena_tree::ArenaTree,
    arena_vertex::{Idx, Label},
    instance::Instance,
};

use super::stack_arena_tree::ArenaTreeMaximalExt;

/// Stack frame used by the depth-bounded search
pub struct Frame {
    pub checkpoint: usize,
    pub branches: Vec<Branch>,
    pub next_branch: usize,
}

/// Branching decisions explored by the search algorithm.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Branch {
    /// Cut one leaf in each active forest
    CutLeaf { leaf_f1: Idx, leaf_f2: Idx },
    /// Cut all pendant subtrees on the path between two leaves in F1
    CutPathBetween { leaf_a: Idx, leaf_b: Idx, lca: Idx },
}

/// Identifiers for the two currently active forests
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Forest {
    F1,
    F2,
}
use Forest::*;

/// Entries storing applied modifications to a forest, with information needed to perform backtracking.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Undo {
    CutAndContracted {
        forest: Forest,
        cut: Idx,
        contracted_parent: Idx,
    },
    MergedLeaves {
        forest: Forest, // actually always done in two forests simultaneously
        new_leaf: Idx,
        old_leaf_a: Idx,
        old_leaf_b: Idx,
    },
    UnmergedLeaves {
        forest: Forest, // actually always done in two forests simultaneously
        old_label: Label,
        old_leaf: Idx,
        new_leaf_a: Idx,
        new_leaf_b: Idx,
    },
    /// Advanced the active forest to the next
    NextTree,
}
use Undo::*;

/// Search state of the depth-bounded search exploration
///
/// Stores:
/// - the current instance (list of forests),
/// - the indices of the active forest pair,
/// - an undo stack used for backtracking.
pub struct State {
    pub instance: Instance,
    pub f1_idx: usize,
    pub f2_idx: usize,
    pub undos: Vec<Undo>,
}

impl State {
    pub fn new_from_instance(instance: Instance) -> State {
        State {
            instance,
            f1_idx: 0,
            f2_idx: 1,
            undos: vec![],
        }
    }

    /// Restore the state to a previously recorded checkpoint, where the given checkpoint is the
    /// length of the undo stack at that point
    pub fn rollback(&mut self, checkpoint: usize) {
        debug_assert!(self.undos.len() >= checkpoint);

        while self.undos.len() > checkpoint {
            let undo = self.undos.pop().unwrap();

            match undo {
                CutAndContracted {
                    forest,
                    cut,
                    contracted_parent,
                } => self
                    .forest_mut(forest)
                    .undo_cut_and_contract(cut, contracted_parent),
                MergedLeaves {
                    forest,
                    new_leaf,
                    old_leaf_a,
                    old_leaf_b,
                    ..
                } => self
                    .forest_mut(forest)
                    .unmerge_node(new_leaf, old_leaf_a, old_leaf_b),
                UnmergedLeaves {
                    forest,
                    old_label,
                    old_leaf,
                    new_leaf_a,
                    new_leaf_b,
                    ..
                } => self
                    .forest_mut(forest)
                    .apply_merge(old_leaf, new_leaf_a, new_leaf_b, old_label),
                NextTree => self.undo_next_tree(),
            }
        }
    }

    /// Exhaustively and recursively merge all common sibling leaves between f1 and f2.
    ///
    /// # Notes
    ///
    /// Useful when there are possibly many unknown unmerged sibling leaves.
    /// After performing some known cuts or merges, use 'recursively_try_merge_siblings' instead.
    pub fn exhaustively_merge_agreeing_subtrees(&mut self) {
        #[cfg(feature = "assert_validity")]
        self.assert_svts_synced();

        let siblings_f1 = self.f1().find_all_cherries();
        let siblings_f2 = self.f2().find_all_cherries();

        for ((leaf1_label, leaf2_label), &parent_f2) in siblings_f2.iter() {
            if let Some(&parent_f1) = siblings_f1.get(&(*leaf1_label, *leaf2_label)) {
                self.try_recursively_merge_siblings(parent_f1, parent_f2);
            }
        }

        #[cfg(feature = "assert_validity")]
        self.assert_fully_reduced();
    }

    fn try_recursively_merge_siblings(&mut self, parent_f1: Idx, parent_f2: Idx) {
        // check if parent_f1 is a cherry in f1
        let Some(((leaf_a_f1_label, leaf_a_f1), (leaf_b_f1_label, leaf_b_f1))) = self
            .f1()
            .get_indices_labels_if_sibling_pair(self.f1().get(parent_f1))
        else {
            return;
        };
        // check if parent_f2 is a cherry in f2
        let Some(((leaf_a_f2_label, leaf_a_f2), (leaf_b_f2_label, leaf_b_f2))) = self
            .f2()
            .get_indices_labels_if_sibling_pair(self.f2().get(parent_f2))
        else {
            return;
        };
        // check if the labels match
        if leaf_a_f1_label != leaf_a_f2_label || leaf_b_f1_label != leaf_b_f2_label {
            return;
        }

        // get references to grandparents
        let grandparent_f1_opt = self.f1().get(parent_f1).parent;
        let grandparent_f2_opt = self.f2().get(parent_f2).parent;

        let new_label = leaf_a_f1_label.min(leaf_b_f1_label);
        // actually merge the leaves
        // marking the leaves as merged and converting the parents to leaves with a new label
        // NOTE: the new label is a copy of one of the labels of the children
        self.f1_mut()
            .apply_merge(parent_f1, leaf_a_f1, leaf_b_f1, new_label);

        self.f2_mut()
            .apply_merge(parent_f2, leaf_a_f2, leaf_b_f2, new_label);

        self.undos.push(MergedLeaves {
            forest: F1,
            new_leaf: parent_f1,
            old_leaf_a: leaf_a_f1,
            old_leaf_b: leaf_b_f1,
        });
        self.undos.push(MergedLeaves {
            forest: F2,
            new_leaf: parent_f2,
            old_leaf_a: leaf_a_f2,
            old_leaf_b: leaf_b_f2,
        });

        // recursively check if grandparents can alo be merged
        // or
        // recursively handle possible newly generated single-vertex-trees
        match (grandparent_f1_opt, grandparent_f2_opt) {
            (Some(grandparent_f1), Some(grandparent_f2)) => {
                self.try_recursively_merge_siblings(grandparent_f1, grandparent_f2);
            }
            (None, Some(_)) => {
                // parent_f1 is a new single-vertex-tree in f1 and not in f2
                self.recursively_sync_single_vertex_tree_f1_to_f2(parent_f1, parent_f2);
            }
            (Some(_), None) => {
                // parent_f2 is a new single-vertex-tree in f2 and not in f1
                self.recursively_sync_single_vertex_tree_f2_to_f1(parent_f1, parent_f2);
            }
            _ => (),
        }
    }

    fn recursively_sync_single_vertex_tree_f1_to_f2(&mut self, _leaf_f1: Idx, leaf_f2: Idx) {
        // Takes a single-vertex-tree in f1 and a leaf with the same label with parent in f2
        // and cuts off the leaf in f2, then recursively merges new sibling leaves created by the cut

        // store the grandparent of f2 for after the cut and the sibling of the label
        let grandparent_f2_opt = self
            .f2()
            .get(leaf_f2)
            .parent
            .and_then(|parent_f2| self.f2().get(parent_f2).parent);
        let sibling_f2 = self
            .f2()
            .find_sibling(leaf_f2)
            .expect("a parent and thus a sibling");

        self.instance.forests[self.f2_idx].cut_branch_and_contract_with_undo(
            leaf_f2,
            F2,
            &mut self.undos,
        );

        if let Some(grandparent_f2) = grandparent_f2_opt {
            // if leaf_f2 had a grandparent, this might be a new cherry,
            // possibly prompting a sibling leaf merge with some unknown cherry in f1
            let f2 = self.f2();
            let Some(label_below_grandparent) = f2
                .get(grandparent_f2)
                .children()
                .next()
                .and_then(|child| f2.get(child).label())
            else {
                return;
            };
            let Some(parent_f1) = self
                .f1()
                .get(self.f1().locate_label(label_below_grandparent))
                .parent
            else {
                return;
            };
            self.try_recursively_merge_siblings(parent_f1, grandparent_f2);
        } else {
            // the sibling of svt_label will be a new root after the cut:
            // check if it is a leaf -> then it is a new SVT
            if let Some(sibling_label) = self.f2().get(sibling_f2).label() {
                self.recursively_sync_single_vertex_tree_f2_to_f1(
                    self.f1().locate_label(sibling_label),
                    sibling_f2,
                );
            }
        }
    }

    fn recursively_sync_single_vertex_tree_f2_to_f1(&mut self, leaf_f1: Idx, _leaf_f2: Idx) {
        // Takes a single-vertex-tree in f2 and a leaf with the same label with parent in f1
        // and cuts off the leaf in f1, then recursively merges new sibling leaves created by the cut

        // store the grandparent of f1 and the sibling for after the cut
        let grandparent_f1_opt = self
            .f1()
            .get(leaf_f1)
            .parent
            .and_then(|parent_f1| self.f1().get(parent_f1).parent);
        let sibling_f1 = self
            .f1()
            .find_sibling(leaf_f1)
            .expect("a parent and thus a sibling");

        self.instance.forests[self.f1_idx].cut_branch_and_contract_with_undo(
            leaf_f1,
            F1,
            &mut self.undos,
        );

        if let Some(grandparent_f1) = grandparent_f1_opt {
            // if leaf_f1 had a grandparent, this might be a new cherry,
            // possibly prompting a sibling leaf merge with some unknown cherry in f2
            let f1 = self.f1();
            let Some(label_below_grandparent) = f1
                .get(grandparent_f1)
                .children()
                .next()
                .and_then(|child| f1.get(child).label())
            else {
                return;
            };
            let Some(parent_f2) = self
                .f2()
                .get(self.f2().locate_label(label_below_grandparent))
                .parent
            else {
                return;
            };
            self.try_recursively_merge_siblings(grandparent_f1, parent_f2);
        } else {
            // the sibling of svt_label will be a new root after the cut:
            // check if it is a leaf -> then it is a new SVT
            if let Some(sibling_label) = self.f1().get(sibling_f1).label() {
                self.recursively_sync_single_vertex_tree_f1_to_f2(
                    sibling_f1,
                    self.f2().locate_label(sibling_label),
                );
            }
        }
    }

    /// Syncs all single vertex trees between F1 and F2
    ///
    /// # Notes
    ///
    /// Does NOT perform recursive sibling-leaf-merging afterwards.
    pub fn sync_single_vertex_trees(&mut self) {
        let [f1, f2] = self.instance.forests.get_disjoint_mut([self.f1_idx, self.f2_idx]).expect("Expected f1_idx and f2_idx to always be distinct and in bound when grabbing mutable references");

        let svts_f1: HashSet<Label> = f1
            .roots
            .iter()
            .filter_map(|&root| f1.get(root).label())
            .collect();
        let svts_f2: HashSet<Label> = f2
            .roots
            .iter()
            .filter_map(|&root| f2.get(root).label())
            .collect();
        let mut svts_to_sync: Vec<Label> =
            svts_f1.symmetric_difference(&svts_f2).copied().collect();

        while let Some(svt_label) = svts_to_sync.pop() {
            let f1_label_idx = f1.locate_label(svt_label);
            let f2_label_idx = f2.locate_label(svt_label);

            if let Some(f1_parent) = f1.get(f1_label_idx).parent {
                if f1.get(f1_parent).parent.is_none() {
                    // the sibling of svt_label will be a new root after the cut:
                    // check if it is a leaf -> then it is a new SVT
                    let sibling = f1
                        .find_sibling(f1_label_idx)
                        .expect("a sibling, because it has a parent");
                    if let Some(sibling_label) = f1.get(sibling).label() {
                        svts_to_sync.push(sibling_label);
                    }
                }
                f1.cut_branch_and_contract_with_undo(f1_label_idx, F1, &mut self.undos);
            } else if let Some(f2_parent) = f2.get(f2_label_idx).parent {
                if f2.get(f2_parent).parent.is_none() {
                    // the sibling of svt_label will be a new root after the cut:
                    // check if it is a leaf -> then it is a new SVT
                    let sibling = f2
                        .find_sibling(f2_label_idx)
                        .expect("a sibling, because it has a parent");
                    if let Some(sibling_label) = f2.get(sibling).label() {
                        svts_to_sync.push(sibling_label);
                    }
                }
                f2.cut_branch_and_contract_with_undo(f2_label_idx, F2, &mut self.undos);
            }
        }
        #[cfg(feature = "assert_validity")]
        self.assert_svts_synced();
    }

    #[inline]
    pub fn ord(&self) -> usize {
        // the order of the solution is always equal to ord(f1),
        // because f1 always has the most components
        self.f1().ord()
    }

    #[inline]
    fn forest_mut(&mut self, forest: Forest) -> &mut ArenaTree {
        match forest {
            F1 => self.f1_mut(),
            F2 => self.f2_mut(),
        }
    }

    #[inline]
    pub fn f1(&self) -> &ArenaTree {
        &self.instance.forests[self.f1_idx]
    }
    #[inline]
    pub fn f2(&self) -> &ArenaTree {
        &self.instance.forests[self.f2_idx]
    }

    #[inline]
    pub fn f1_mut(&mut self) -> &mut ArenaTree {
        &mut self.instance.forests[self.f1_idx]
    }

    #[inline]
    pub fn f2_mut(&mut self) -> &mut ArenaTree {
        &mut self.instance.forests[self.f2_idx]
    }

    pub fn current_forest_is_completed(&self) -> bool {
        self.f2().find_a_sibling_pair().is_none()
    }

    #[inline]
    pub fn is_on_last_tree(&self) -> bool {
        assert!(self.f2_idx < self.instance.forests.len());
        self.f2_idx == self.instance.forests.len() - 1
    }

    #[inline]
    pub fn is_solved(&self) -> bool {
        self.is_on_last_tree()
            && self.f1().find_a_sibling_pair().is_none()
            && self.f2().find_a_sibling_pair().is_none()
    }

    pub fn extract_solution(mut self) -> ArenaTree {
        self.unmerge_f1_with_undos();
        self.instance.forests.swap_remove(self.f1_idx)
    }

    pub fn load_next_tree(&mut self) {
        debug_assert!(self.f2_idx + 1 < self.instance.forests.len());
        self.unmerge_f1_with_undos();

        self.f2_idx += 1;

        self.undos.push(Undo::NextTree);
    }

    pub fn undo_next_tree(&mut self) {
        self.f2_idx -= 1;
        debug_assert!(self.f2_idx > self.f1_idx);
    }

    fn unmerge_f1_with_undos(&mut self) {
        let leaf_merges: Vec<Undo> = self
            .undos
            .iter()
            .rev()
            .take_while(|&undo| !matches!(undo, NextTree))
            .filter(|&undo| matches!(undo, MergedLeaves { forest: F1, .. }))
            .copied()
            .collect();
        self.instance.forests[self.f1_idx].unmerge_with_undo(leaf_merges, &mut self.undos);
    }

    #[cfg(feature = "assert_validity")]
    fn assert_svts_synced(&self) {
        use crate::maf_instance::tree_traversal::TreeTraversal;

        let f1 = self.f1();
        let f2 = self.f2();

        for label in f1.iterate_all().labels() {
            let leaf = f1.locate_label(label);
            let is_svt_f1 = f1.get(leaf).parent.is_none();

            let leaf2 = f2.locate_label(label);
            let is_svt_f2 = f2.get(leaf2).parent.is_none();
            assert_eq!(is_svt_f1, is_svt_f2);
        }
    }

    #[cfg(feature = "assert_validity")]
    fn assert_common_cherries_merged(&self) {
        let cherries_f2 = self.f2().find_all_cherries();

        for (_, a, b) in self.f1().iterate_cherries() {
            assert!(!cherries_f2.contains_key(&(a, b)));
            assert!(!cherries_f2.contains_key(&(b, a)));
        }
    }

    #[cfg(feature = "assert_validity")]
    pub fn assert_fully_reduced(&self) {
        self.assert_svts_synced();
        self.assert_common_cherries_merged();
    }
}
