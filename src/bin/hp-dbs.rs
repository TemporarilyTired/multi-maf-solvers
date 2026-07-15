use multi_maf_solvers::alg_hp_dbs::solve_binary_search_dbs;
use multi_maf_solvers::common::reading::PrintableSolution;
use multi_maf_solvers::common::reading::read_to_my_instance;

fn main() {
    let instance = read_to_my_instance();

    println!(
        "# Succesfully read {} forests with {} leaves",
        instance.forests.len(),
        instance.num_leaves
    );

    println!("# Starting hitting pair depth bounded search with binary search on solution size");

    let maf = solve_binary_search_dbs(instance);

    maf.print_newick_strings();
}
