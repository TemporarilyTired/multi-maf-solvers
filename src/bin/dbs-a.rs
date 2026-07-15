use clap::Parser;
use multi_maf_solvers::alg_dbs_a::solve;
use multi_maf_solvers::alg_dbs_a::variants::*;
use multi_maf_solvers::common::reading::PrintableSolution;
use multi_maf_solvers::common::reading::read_to_my_instance;

fn main() {
    let args = Args::parse();

    let instance = read_to_my_instance();

    println!(
        "# Succesfully read {} forests with {} leaves",
        instance.forests.len(),
        instance.num_leaves
    );

    println!("# Starting DBS-A: maximal-agreement forest depth-bounded search");
    let maf = match args.variant {
        VariantKind::Baseline => solve::<Baseline>(instance),
        VariantKind::BetterBranching => solve::<OnlyBetterBranchingOnLastTree>(instance),
        VariantKind::LbCalculation => solve::<OnlyLbCalculation>(instance),
        VariantKind::BsOnSolution => solve::<OnlyBsOnSolution>(instance),
        VariantKind::Caching => solve::<OnlyCaching>(instance),
        VariantKind::GlobalPreprocessing => solve::<OnlyGlobalPreprocessing>(instance),
        VariantKind::AllExceptCaching => solve::<AllExceptCaching>(instance),
        VariantKind::All => solve::<AllFeatures>(instance),
    };

    maf.print_newick_strings();
}
