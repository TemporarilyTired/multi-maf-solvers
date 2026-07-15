# Multi-MAF Solvers
A repository containing multiple different algorithms for solving the maximum agreement forests problem for multiple rooted binary phylogenetic trees. See the submissions to the [PACE 2026 Challenge](https://pacechallenge.org/2026/) for more (some significantly faster) algorithms for this problem.

This repository accompanies the master's thesis *"The Ingredients of Depth-Bounded Search Algorithms for the Maximum Agreement Forest Problem on Multiple Trees"* (Utrecht University, 2026), which will be published in the near future. This thesis will describe the algorithms, correctness proofs, implementation choices, and experimental results in detail.
 
## Algorithms
 
This repository contains implementations of four exact algorithms for the Maximum Agreement Forest (MAF) problem on multiple rooted binary trees:
 
| Algorithm | Binary name | Description | Reference |
|---|---|---|---|
| DBS-A | `dbs-a` | A depth-bounded search algorithm based on iterating maximal agreement forests. | Shi, Wang, Chen, Feng & Guo (2014) |
| DBS-LSI | `dbs-lsi` | A depth-bounded search algorithm based on the label-set isomorphism property. | Shi, Chen, Feng & Wang (2018) |
| HP-DBS | `hp-dbs` | A new depth-bounded search algorithm introduced in this thesis, extending DBS-LSI with "hitting pairs" to allow for better branching choices and extra reduction rules. | This thesis |
| AF-DFS | `af-dfs` | A new algorithm introduced in this thesis that iterates all agreement forests via depth-first search on their canonical representation; effective on instances with a very large optimal solution relative to the number of leaves. | This thesis |
 
## Installation
 
See [`INSTALL.md`](./INSTALL.md) for detailed build instructions, including how to select which algorithm to build and choose between variants of an algorithm.
 
## Quick start
 
```bash
# Clone the repository
git clone https://github.com/TemporarilyTired/multi-maf-solvers.git

cd multi-maf-solvers
```
 
# Build a specific algorithm in release mode
```bash
cargo build --release --bin hp-dbs
```
 
# Run it on an instance
```bash
./target/release/hp-dbs < path/to/instance.txt
```
 
## Input format
 
The algorithms accept instances in the format defined in the PACE 2026 Challenge. See the [PACE 2026 problem specification](https://pacechallenge.org/2026/format/) for details.
 
## Repository structure
 
```
.
├── src/
│   ├── bin/            # One entry point per algorithm (dbs-a.rs, dbs-lsi.rs, hp-dbs.rs, af-dfs.rs)
│   ├── alg_*/          # Module containing code for a specific algorithm (e.g., alg_af_dfs/ for af-dfs)
│   └── ...             # Shared library code (forest representation, preprocessing, etc.)

├── Cargo.toml
├── README.md
└── INSTALL.md
```

## License
 
GNU General Public License version 3
