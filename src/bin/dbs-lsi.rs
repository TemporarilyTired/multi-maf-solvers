use multi_maf_solvers::alg_dbs_lsi::solve;
use multi_maf_solvers::common::reading::PrintableSolution;
use multi_maf_solvers::common::reading::read_to_my_instance;

fn main() {
    #[cfg(feature = "logging")]
    {
        // Emit information on which features are enabled:
        #[cfg(not(feature = "LSI_DISABLE_BINARY_SEARCH_ON_SOL"))]
        println!("#s binary_search_enabled 1");
        #[cfg(not(feature = "LSI_DISABLE_CLUSTER_REDUCTION"))]
        println!("#s cluster_reduction_enabled 1");
        #[cfg(not(feature = "LSI_DISABLE_LB_CALC"))]
        println!("#s lb_calc_enabled 1");
        #[cfg(not(feature = "LSI_DISABLE_MODIFIED_BRANCHING"))]
        println!("#s modified_branching_enabled 1");
        #[cfg(not(feature = "LSI_DISABLE_GLOBAL_PREPROCESSING"))]
        println!("#s global_preprocessing_enabled 1");
    }

    let instance = read_to_my_instance();

    println!(
        "# Succesfully read {} forests with {} leaves",
        instance.forests.len(),
        instance.num_leaves
    );

    println!("# Starting LSI depth-bounded search");

    let maf = solve(instance);

    maf.print_newick_strings();
}
