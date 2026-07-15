use super::lsi_instance::LsiInstance;
use crate::maf_instance::instance::Instance;

pub enum BranchResult {
    NotApplicable(Instance),
    Branch2([LsiInstance; 2]),
    Branch3([LsiInstance; 3]),
}
