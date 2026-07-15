use crate::{
    maf_instance::{arena_tree::ArenaTree, arena_vertex::Label, instance::Instance},
    preprocessing::instance_ext::ReductionsInstanceExt,
};

pub fn reduce_and_solve<S>(mut instance: Instance, solve_fn: &S) -> ArenaTree
where
    S: Fn(Instance) -> ArenaTree,
{
    // INFO: apply:
    // - sibling leaf merge,
    // - single vertex tree syncing, and
    // - reduction rule 2.2.1 (from DBS-LSI)
    // - reduction rule 1 (from DBS-LSI)
    // to exhaustively to reduce the input trees to (smaller) forests
    let performed_reductions = instance.fully_reduce();

    let clusters = instance.clone().split_into_clusters();
    let mut cluster_mafs = vec![];

    for cluster in clusters {
        cluster_mafs.push(solve_cluster(cluster, solve_fn));
    }

    let mut maf = cluster_mafs
        .into_iter()
        .reduce(|maf1, maf2| maf1.join_with(maf2))
        .unwrap_or_default();

    for reduction in performed_reductions.iter().rev() {
        maf.undo_reduction(reduction.clone());
    }
    maf
}

/// Solve a fully reduced cluster of an instance
fn solve_cluster<S>(cluster: Instance, solve_fn: &S) -> ArenaTree
where
    S: Fn(Instance) -> ArenaTree,
{
    // if the instance has a useful common cluster, split it into clusters and solve separately
    // if the addition of the dummy leaf does not increase the MAF in both parts: some component
    // can span the cut edge of the subtree
    // So: test if opt(above)+opt(below) <= k, if so: return Some(_)
    // otherwise: test if opt(above)+opt(below) > k+1, if so: return None
    // At this point opt(above)+opt(below) == k+1,
    // If opt(above with dummy) == opt(above)  (we can just check for   opt(above w dummy) <= opt(above))
    // and opt(below with dummy) == opt(below)  (we can just check for   opt(below w dummy) <= opt(below))
    // Then: return Some(opt(above) + opt(below) - 1)
    // otherwise: return None

    if let Some((
        above,
        below,
        (above_w_dummy, used_dummy_above),
        (below_w_dummy, used_dummy_below),
    )) = cluster.find_clusters_w_dummy()
    {
        return solve_split_with_dummy(
            above,
            below,
            (above_w_dummy, used_dummy_above),
            (below_w_dummy, used_dummy_below),
            solve_fn,
        );
    }

    #[cfg(feature = "logging")]
    println!(
        "# minimial cluster (ord={}) leaves: \t{}",
        cluster.ord(),
        cluster.num_leaves,
    );

    solve_fn(cluster)
}

/// Solve an instance split into 4 by subtree cluster reduction
fn solve_split_with_dummy<S>(
    above: Instance,
    below: Instance,
    (above_with_dummy, used_dummy_above): (Instance, Label),
    (below_with_dummy, used_dummy_below): (Instance, Label),
    solve_fn: &S,
) -> ArenaTree
where
    S: Fn(Instance) -> ArenaTree,
{
    let mut a_dummy = reduce_and_solve(above_with_dummy, solve_fn);
    if a_dummy
        .get(a_dummy.locate_label(used_dummy_above))
        .parent
        .is_some()
    {
        let mut b_dummy = reduce_and_solve(below_with_dummy, solve_fn);
        if b_dummy
            .get(b_dummy.locate_label(used_dummy_below))
            .parent
            .is_some()
        {
            // We can construct a solution of order |a_dummy| + |b_dummy| -1
            // But this can still be one larger than |a| + |b| in the case:
            //      |a_dummy| = |a|+1 and |b_dummy| = |b|+1
            // So we need to calculate |b| and possibly |a| too
            let b = reduce_and_solve(below, solve_fn);
            if b.ord() == b_dummy.ord() {
                // |b_dummy| = |b|, thus |a_dummy| + |b_dummy| - 1 <= |a| + |b|
                return a_dummy.join_at_dummy(b_dummy, used_dummy_above, used_dummy_below);
            }
            let a = reduce_and_solve(above, solve_fn);
            if a.ord() == a_dummy.ord() {
                // |a_dummy| = |a|, thus |a_dummy| + |b_dummy| - 1 <= |a| + |b|
                return a_dummy.join_at_dummy(b_dummy, used_dummy_above, used_dummy_below);
            }
            return a.join_with(b);
        }
        // b_dummy contains an optimal solution for below
        b_dummy.remove_svt(used_dummy_below);
        let a = reduce_and_solve(above, solve_fn);
        return a.join_with(b_dummy);
    }

    // a_dummy contains an optimal solution for above
    a_dummy.remove_svt(used_dummy_above);
    let b = reduce_and_solve(below, solve_fn);
    a_dummy.join_with(b)
}
