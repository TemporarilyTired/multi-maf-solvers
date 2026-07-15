use std::collections::HashMap;

use super::state::{
    Forest,
    Undo::{self, *},
};
use crate::maf_instance::{
    arena_tree::ArenaTree,
    arena_vertex::{Idx, Label, Node, NodeData::*, Status::*},
};

/// Extension trait for ArenaTree defining functions necessary for this specific algorithm.
/// Mainly contains operations that record their changes to the undo stack.
pub trait ArenaTreeMaximalExt {
    fn unmerge_with_undo(&mut self, merged_leaves: Vec<Undo>, undos: &mut Vec<Undo>);
    fn get_lca_and_path_len(&self, leaf_a: Idx, leaf_b: Idx) -> Option<(Idx, usize)>;
    fn get_labels_if_sibling_pair(&self, parent: &Node) -> Option<(Label, Label)>;
    fn get_indices_labels_if_sibling_pair(
        &self,
        parent: &Node,
    ) -> Option<((Label, Idx), (Label, Idx))>;
    fn find_a_sibling_pair(&self) -> Option<((Label, Idx), (Label, Idx))>;
    fn find_all_cherries(&self) -> HashMap<(Label, Label), Idx>;
    fn unmerge_node(&mut self, parent_leaf: Idx, old_leaf_a: Idx, old_leaf_b: Idx);
    fn cut_branch_and_contract_with_undo(
        &mut self,
        cut: Idx,
        forest: Forest,
        undos: &mut Vec<Undo>,
    );
    fn undo_cut_and_contract(&mut self, cut: Idx, contracted_parent: Idx);
}

impl ArenaTreeMaximalExt for ArenaTree {
    /// Cuts of the branch and contracts the parent vertex, without recording it in the undo stack
    fn undo_cut_and_contract(&mut self, cut: Idx, contracted_parent: Idx) {
        debug_assert!(matches!(
            self.get(contracted_parent),
            Node {
                status: Contracted,
                data: Internal { left, right },
                ..
            } if *left == cut || *right == cut
        ));
        debug_assert!(matches!(
            self.get(cut),
            Node {
                status: Present,
                parent: None,
                ..
            }
        ));
        // -- first uncontract the parent node --
        let grandparent_opt = self.get(contracted_parent).parent;

        // find other child of the contracted parent
        let sibling = self
            .get(contracted_parent)
            .children()
            .find(|&v| v != cut)
            .expect(
                "One of the two children of the parent of the cut must be different from the cut",
            );

        // set contracted parent as child of grandparent instead of sibling
        if let Some(grand_parent) = grandparent_opt {
            self.get_mut(grand_parent).replace_child(sibling, contracted_parent)
                .expect("Undoing Cut and Contract failed: couldn't uncontract parent vertex because sibling is not a child of grand_parent");
        } else {
            // if contracted parent was a root before: make it root again and make sibling not a
            // root anymore
            debug_assert!(self.roots.contains(&sibling));
            self.roots.remove(&sibling);
            self.roots.insert(contracted_parent);
        }

        // uncontract parent
        self.get_mut(contracted_parent).status = Present;
        self.get_mut(sibling).parent = Some(contracted_parent);

        // -- then attach child to parent --
        self.get_mut(cut).parent = Some(contracted_parent);
        debug_assert!(self.roots.contains(&cut));
        self.roots.remove(&cut);
    }

    /// Cuts of the branch and contracts the parent vertex, recording it in the undo stack
    fn cut_branch_and_contract_with_undo(
        &mut self,
        cut: Idx,
        forest: Forest,
        undos: &mut Vec<Undo>,
    ) {
        // find parent
        let parent = self
            .get(cut)
            .parent
            .expect("Branch cut failed: tried to cut off a vertex without a parent");

        // find grand_parent (may be None)
        let grand_parent_opt = self.get(parent).parent;

        // find other child of the parent
        let sibling = self.get(parent).children().find(|&v| v != cut).expect(
            "One of the two children of the parent of the cut must be different from the cut",
        );

        // cut off the branch
        self.get_mut(cut).parent = None;
        self.roots.insert(cut);

        // contract parent vertex
        self.get_mut(parent).status = Contracted;
        self.get_mut(sibling).parent = grand_parent_opt;

        if let Some(grand_parent) = grand_parent_opt {
            self.get_mut(grand_parent).replace_child(parent, sibling)
                .expect("Branch cut failed (not cleanly): couldn't contract parent vertex because parent is not a child of grand_parent (tree representation is invalid)");
        } else {
            // if contracted vertex was a root: make their child a root instead
            debug_assert!(self.roots.contains(&parent));
            self.roots.remove(&parent);
            self.roots.insert(sibling);
        }

        undos.push(CutAndContracted {
            forest,
            cut,
            contracted_parent: parent,
        });
    }

    /// Unmerge a previously merged leaf pair
    fn unmerge_node(&mut self, parent_leaf: Idx, old_leaf_a: Idx, old_leaf_b: Idx) {
        debug_assert!(
            matches!(self.get(old_leaf_a),Node {status: MergedIntoParent, data: Leaf { .. }, parent: Some(parent)} if *parent == parent_leaf)
        );
        debug_assert!(
            matches!(self.get(old_leaf_b),Node {status: MergedIntoParent, data: Leaf { .. }, parent: Some(parent)} if *parent == parent_leaf)
        );
        debug_assert!(matches!(
            self.get(parent_leaf),
            Node {
                status: Present,
                data: Leaf { .. },
                ..
            }
        ));

        // update leaf map entries
        self.leaf_map
            .remove(&self.get(parent_leaf).label().expect("leaf to have a label"));
        self.leaf_map.insert(
            self.get(old_leaf_a).label().expect("leaf to have a label"),
            old_leaf_a,
        );
        self.leaf_map.insert(
            self.get(old_leaf_b).label().expect("leaf to have a label"),
            old_leaf_b,
        );

        // split parent node
        self.get_mut(old_leaf_a).status = Present;
        self.get_mut(old_leaf_b).status = Present;

        self.get_mut(parent_leaf).data = Internal {
            left: old_leaf_a,
            right: old_leaf_b,
        };
    }

    /// Unmerge the leaves of the forest, recording perfomed unmerges in the undo stack
    fn unmerge_with_undo(&mut self, merged_leaves: Vec<Undo>, undos: &mut Vec<Undo>) {
        for undo in merged_leaves {
            let MergedLeaves {
                forest,
                new_leaf,
                old_leaf_a,
                old_leaf_b,
            } = undo
            else {
                unreachable!()
            };
            let old_label = self.get(new_leaf).label().expect("leaf to have a label");
            undos.push(UnmergedLeaves {
                forest,
                old_label,
                old_leaf: new_leaf,
                new_leaf_a: old_leaf_a,
                new_leaf_b: old_leaf_b,
            });
            self.unmerge_node(new_leaf, old_leaf_a, old_leaf_b);
        }
    }

    /// Finds all siblings pairs in all components.
    /// Returns hashMap from labels of pair -> idx of common parent
    fn find_all_cherries(&self) -> HashMap<(Label, Label), Idx> {
        let mut siblings = HashMap::<(Label, Label), Idx>::new();

        for (idx, vertex) in self.arena.iter().enumerate() {
            self.get_labels_if_sibling_pair(vertex).inspect(|labels| {
                siblings.insert(*labels, Idx::try_from(idx).unwrap());
            });
        }
        siblings
    }

    #[inline]
    fn find_a_sibling_pair(&self) -> Option<((Label, Idx), (Label, Idx))> {
        self.arena
            .iter()
            .find_map(|vertex| self.get_indices_labels_if_sibling_pair(vertex))
    }

    #[inline]
    fn get_indices_labels_if_sibling_pair(
        &self,
        parent: &Node,
    ) -> Option<((Label, Idx), (Label, Idx))> {
        // If the given vertex has two leaves as children, return their labels sorted and their indices, otherwise None
        if parent.status == Present
            && let Internal { left, right } = parent.data
            && let Node {
                status: Present,
                data: Leaf { label: label_left },
                ..
            } = self.get(left)
            && let Node {
                status: Present,
                data: Leaf { label: label_right },
                ..
            } = self.get(right)
        {
            if label_left <= label_right {
                return Some(((*label_left, left), (*label_right, right)));
            } else {
                return Some(((*label_right, right), (*label_left, left)));
            }
        };
        None
    }

    /// If the given vertex has two leaves as children, return their labels sorted, otherwise None
    #[inline]
    fn get_labels_if_sibling_pair(&self, parent: &Node) -> Option<(Label, Label)> {
        if parent.status == Present
            && let Internal { left, right } = parent.data
            && let Node {
                status: Present,
                data: Leaf { label: label1 },
                ..
            } = self.get(left)
            && let Node {
                status: Present,
                data: Leaf { label: label2 },
                ..
            } = self.get(right)
        {
            let min_label = label1.min(label2);
            let max_label = label1.max(label2);
            return Some((*min_label, *max_label));
        };
        None
    }

    #[inline]
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
}
