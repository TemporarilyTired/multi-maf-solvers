use std::collections::{HashMap, HashSet};

use crate::maf_instance::{
    arena_tree::ArenaTree,
    arena_vertex::{Idx, Label, NodeData::*, Status::*},
    tree_traversal::TreeTraversal,
};

pub trait ReductionsArenaTreeExt {
    fn copy_cluster(&self, cluster_labels: &HashSet<Label>) -> ArenaTree;
}

impl ReductionsArenaTreeExt for ArenaTree {
    fn copy_cluster(&self, cluster_labels: &HashSet<Label>) -> ArenaTree {
        let roots: HashSet<Idx> = cluster_labels
            .iter()
            .map(|&label| self.find_root_of(self.locate_label(label)))
            .collect();
        let indices_of_subtree: HashSet<Idx> = roots
            .iter()
            .flat_map(|&root| self.dfs_from(root).indices())
            .collect();

        let mut translation: HashMap<Idx, Idx> = HashMap::new();

        let mut new_arena = vec![];
        for (idx, node) in self.arena.iter().cloned().enumerate() {
            if node.status != Present {
                continue;
            }
            if indices_of_subtree.contains(&(idx as Idx)) {
                let new_idx = new_arena.len() as Idx;
                new_arena.push(node);
                translation.insert(idx as Idx, new_idx);
            }
        }

        let new_roots: HashSet<Idx> = new_arena
            .iter()
            .enumerate()
            .filter_map(|(idx, node)| node.parent.is_none().then_some(idx as Idx))
            .collect();

        let new_leaf_map: HashMap<Label, Idx> = new_arena
            .iter()
            .enumerate()
            .filter_map(|(idx, node)| node.label().map(|label| (label, idx as Idx)))
            .collect();

        for node in new_arena.iter_mut() {
            node.parent = node.parent.map(|p| translation[&p]);
            node.data = match node.data.clone() {
                Internal { left, right } => Internal {
                    left: translation[&left],
                    right: translation[&right],
                },
                leaf => leaf,
            };
        }

        ArenaTree {
            arena: new_arena,
            roots: new_roots,
            leaf_map: new_leaf_map,
        }
    }
}
