# Installation
 
## Prerequisites
 
- [Rust](https://www.rust-lang.org/tools/install), version 1.93.1 or later
- Cargo (installed automatically alongside Rust via `rustup`)
Check your installed versions with:
 
```bash
rustc --version
cargo --version
```
 
## Cloning the repository
 
```bash
git clone https://github.com/<your-username>/multi-maf-solvers.git
cd multi-maf-solvers
```
 
## Building an algorithm
 
Each algorithm is compiled as its own binary. Use `--bin` with `cargo build` (or `cargo run`) to select which one to build:
 
```bash
cargo build --release --bin <binary-name>
```
 
Available binaries:
 
| Binary name | Algorithm |
|---|---|
| `dbs-a`   | DBS-A |
| `dbs-lsi` | DBS-LSI |
| `hp-dbs`  | HP-DBS |
| `af-dfs`  | AF-DFS |
 
Or build all binaries at once:
 
```bash
cargo build --release --bins
```
 
The resulting executables are placed in `target/release/`.
 
## Running a binary
 
```bash
./target/release/<binary-name> [OPTIONS] < instance.txt
```
 <!-- TODO:  ^ -->

## Variants

The DBS-LSI and DBS-A algorithm have different versions implemented; for DBS-LSI, the exact version can be chosen during compile time using Cargo features. For DBS-A, the binary accepts a runtime argument instead.

### DBS-LSI variants
 
For the DBS-LSI algorithm, the implemented modifications described in the thesis (e.g., cluster reduction, lower-bound pruning, etc.) can be disabled with Cargo [features](https://doc.rust-lang.org/cargo/reference/features.html) rather than always being compiled in. This keeps binary lean and lets you build the exact variant you want to run or benchmark.
 
Enable a feature with `--features`:
 
```bash
cargo build --release --bin dbs-lsi --features LSI_DISABLE_LB_CALC
```


#### Available features
 
| Feature | Description |
|---|---|
| `LSI_DISABLE_LB_CALC` | Disables/removes the lower bound computations in recursive calls of DBS-LSI. |
| `LSI_DISABLE_CLUSTER_REDUCTION` | Disables the cluster reduction rules in recursive calls of DBS-LSI. |
| `LSI_DISABLE_MODIFIED_BRANCHING` | Disables the modifications made to the original branching scheme of DBS-LSI. |
| `LSI_DISABLE_GLOBAL_PREPROCESSING` | Removes the global preprocessing step that is otherwise performed at the start of the algorithm. |
| `LSI_DISABLE_BINARY_SEARCH_ON_SOL` | Disables the binary search on the order of the MAF, instead trying all parameters k iteratively. |

### DBS-A variants
 
For the DBS-A algorithm, the implemented different variants described and tested in the thesis (e.g., caching, lower-bound pruning, etc.) can be executed using a runtime argument, called `variant`.
 
Choose a variant with `--variant`:
 
```bash
echo "\#p 2 5
((4,(3,(5,1))),2);
((4,(3,(1,2))),5);" | ./target/release/dbs-a --variant global-preprocessing
```

or, if you have an instance as a file `instance.txt`:

```bash
./target/release/dbs-a --variant all < instance.txt
```

#### Available variants
 
| Variant | Description |
|---|---|
| `baseline` | Contains none of the optional modifications made to the original algorithm. |
| `better-branching` | Slightly modifies the branching scheme. |
| `lb-calculation` | Enables the lower bound computations in recursive calls of DBS-A. |
| `bs-on-solution` | Enables the binary search on the order of the MAF, instead of trying all parameters k iteratively. |
| `caching` | Enables the caching of visited states to prune previously visited states. |
| `global-preprocessing` | Enables a global preprocessing step that is performed at the start of the algorithm. |
| `all-except-caching` | Enables all modifications mentioned above except for caching. |
| `all` | Enables all modifications mentioned. |

### Debugging features

The compilation also accepts two general Cargo features used in debugging and in the comparison between algorithms:

- `assert_validity`: this feature is only meant for debugging purposes to enable code that performs checks on the data structures to, e.g., determine the origin of an incorrect state in the algorithm. This feature severely affects the performance.
- `logging`: this feature is meant to track performance data during runtime. The resulting output is not meant to be fully human-readable and can at times lead to huge amounts of output.
