use crate::common::generation_set::GenerationSet;
use crate::maf_instance::arena_vertex::Idx;

pub struct R1Workspace {
    /// label to root idx in "other" forest; indexed by label
    pub label_to_root: Vec<Idx>,
    /// label-set of a component; indexed by label
    pub component_labels: GenerationSet,
    /// set containing the labels of a subtree; indexed by label
    pub subtree_labels: GenerationSet,
    /// whether we've already checked a root from another forest for the current vertex; indexed by node idx
    pub checked_roots: GenerationSet,
}

impl R1Workspace {
    pub fn new(max_label: usize, max_node: usize) -> Self {
        Self {
            label_to_root: vec![Idx::MAX; max_label + 1],
            component_labels: GenerationSet::new(max_label + 1),
            subtree_labels: GenerationSet::new(max_label + 1),
            checked_roots: GenerationSet::new(max_node + 1),
        }
    }
}
