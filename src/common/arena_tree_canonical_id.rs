use std::collections::HashMap;

use crate::maf_instance::arena_tree::ArenaTree;
use crate::maf_instance::arena_vertex::*;

pub type CanonId = u32;

pub struct Canonicalizer {
    leaf_table: HashMap<Label, CanonId>,
    internal_table: HashMap<(CanonId, CanonId), CanonId>,
    next_id: CanonId,
}

impl Canonicalizer {
    pub fn new() -> Self {
        Self {
            leaf_table: HashMap::new(),
            internal_table: HashMap::new(),
            next_id: 1,
        }
    }

    fn leaf_id(&mut self, label: Label) -> CanonId {
        *self.leaf_table.entry(label).or_insert_with(|| {
            let id = self.next_id;
            self.next_id += 1;
            id
        })
    }

    fn internal_id(&mut self, a: CanonId, b: CanonId) -> CanonId {
        let key = if a <= b { (a, b) } else { (b, a) };

        *self.internal_table.entry(key).or_insert_with(|| {
            let id = self.next_id;
            self.next_id += 1;
            id
        })
    }
}

impl Default for Canonicalizer {
    fn default() -> Self {
        Self::new()
    }
}

impl ArenaTree {
    pub fn canonical_ids_with(&self, canon: &mut Canonicalizer) -> Vec<CanonId> {
        let mut ids = vec![0; self.arena.len()];

        for &root in &self.roots {
            self.assign_id(root, &mut ids, canon);
        }

        ids
    }

    fn assign_id(&self, idx: Idx, ids: &mut [CanonId], canon: &mut Canonicalizer) -> CanonId {
        let i = idx as usize;

        if ids[i] != 0 {
            return ids[i];
        }

        let node = &self.arena[i];

        let id = match node.data {
            NodeData::Leaf { label } => canon.leaf_id(label),
            NodeData::Internal { left, right } => {
                let l = self.assign_id(left, ids, canon);
                let r = self.assign_id(right, ids, canon);

                canon.internal_id(l, r)
            }
        };

        ids[i] = id;
        id
    }
}
