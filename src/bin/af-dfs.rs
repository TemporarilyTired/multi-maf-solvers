use multi_maf_solvers::alg_af_dfs::solve;
use multi_maf_solvers::common::reading::PrintableSolution;
use multi_maf_solvers::common::reading::read_to_my_instance;

fn main() {
    let instance = read_to_my_instance();

    println!(
        "# Succesfully read {} forests with {} leaves",
        instance.forests.len(),
        instance.num_leaves
    );

    println!("# Starting Agreement Forest-DFS");

    let maf = solve(instance);

    maf.print_newick_strings();
}
