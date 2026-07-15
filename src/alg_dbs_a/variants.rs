use clap::{Parser, ValueEnum};

#[derive(Parser)]
pub struct Args {
    #[arg(long, default_value = "all")]
    pub variant: VariantKind,
}

#[derive(Clone, ValueEnum)]
pub enum VariantKind {
    Baseline,
    BetterBranching,
    LbCalculation,
    BsOnSolution,
    Caching,
    GlobalPreprocessing,
    All,
    AllExceptCaching,
}

pub trait Variant {
    const NAME: &'static str;
    const VARIANT_ID: usize;
    const ENABLE_BETTER_BRANCHING_ON_LAST_TREE: bool;
    const ENABLE_LB_CALCULATION: bool;
    const ENABLE_BS_ON_SOLUTION: bool;
    const ENABLE_CACHING: bool;
    const ENABLE_GLOBAL_PREPROCESSING: bool;
}

pub struct Baseline;
pub struct OnlyBetterBranchingOnLastTree;
pub struct OnlyLbCalculation;
pub struct OnlyBsOnSolution;
pub struct OnlyCaching;
pub struct OnlyGlobalPreprocessing;
pub struct AllExceptCaching;
pub struct AllFeatures;

impl Variant for Baseline {
    const NAME: &'static str = "Baseline";
    const VARIANT_ID: usize = 0;
    const ENABLE_BETTER_BRANCHING_ON_LAST_TREE: bool = false;
    const ENABLE_LB_CALCULATION: bool = false;
    const ENABLE_BS_ON_SOLUTION: bool = false;
    const ENABLE_CACHING: bool = false;
    const ENABLE_GLOBAL_PREPROCESSING: bool = false;
}

impl Variant for OnlyBetterBranchingOnLastTree {
    const NAME: &'static str = "Baseline + better branching on the last tree";
    const VARIANT_ID: usize = 1;
    const ENABLE_BETTER_BRANCHING_ON_LAST_TREE: bool = true;
    const ENABLE_LB_CALCULATION: bool = false;
    const ENABLE_BS_ON_SOLUTION: bool = false;
    const ENABLE_CACHING: bool = false;
    const ENABLE_GLOBAL_PREPROCESSING: bool = false;
}

impl Variant for OnlyLbCalculation {
    const NAME: &'static str = "Baseline + lower bound calculation";
    const VARIANT_ID: usize = 2;
    const ENABLE_BETTER_BRANCHING_ON_LAST_TREE: bool = false;
    const ENABLE_LB_CALCULATION: bool = true;
    const ENABLE_BS_ON_SOLUTION: bool = false;
    const ENABLE_CACHING: bool = false;
    const ENABLE_GLOBAL_PREPROCESSING: bool = false;
}

impl Variant for OnlyBsOnSolution {
    const NAME: &'static str = "Baseline + binary search on solution size";
    const VARIANT_ID: usize = 3;
    const ENABLE_BETTER_BRANCHING_ON_LAST_TREE: bool = false;
    const ENABLE_LB_CALCULATION: bool = false;
    const ENABLE_BS_ON_SOLUTION: bool = true;
    const ENABLE_CACHING: bool = false;
    const ENABLE_GLOBAL_PREPROCESSING: bool = false;
}

impl Variant for OnlyCaching {
    const NAME: &'static str = "Baseline + caching of visited states";
    const VARIANT_ID: usize = 4;
    const ENABLE_BETTER_BRANCHING_ON_LAST_TREE: bool = false;
    const ENABLE_LB_CALCULATION: bool = false;
    const ENABLE_BS_ON_SOLUTION: bool = false;
    const ENABLE_CACHING: bool = true;
    const ENABLE_GLOBAL_PREPROCESSING: bool = false;
}

impl Variant for OnlyGlobalPreprocessing {
    const NAME: &'static str = "Baseline + global preprocessing";
    const VARIANT_ID: usize = 6;
    const ENABLE_BETTER_BRANCHING_ON_LAST_TREE: bool = false;
    const ENABLE_LB_CALCULATION: bool = false;
    const ENABLE_BS_ON_SOLUTION: bool = false;
    const ENABLE_CACHING: bool = false;
    const ENABLE_GLOBAL_PREPROCESSING: bool = true;
}

impl Variant for AllFeatures {
    const NAME: &'static str = "All features";
    const VARIANT_ID: usize = 5;
    const ENABLE_BETTER_BRANCHING_ON_LAST_TREE: bool = true;
    const ENABLE_LB_CALCULATION: bool = true;
    const ENABLE_BS_ON_SOLUTION: bool = true;
    const ENABLE_CACHING: bool = true;
    const ENABLE_GLOBAL_PREPROCESSING: bool = true;
}

impl Variant for AllExceptCaching {
    const NAME: &'static str = "All features except for caching";
    const VARIANT_ID: usize = 7;
    const ENABLE_BETTER_BRANCHING_ON_LAST_TREE: bool = true;
    const ENABLE_LB_CALCULATION: bool = true;
    const ENABLE_BS_ON_SOLUTION: bool = true;
    const ENABLE_CACHING: bool = false;
    const ENABLE_GLOBAL_PREPROCESSING: bool = true;
}
