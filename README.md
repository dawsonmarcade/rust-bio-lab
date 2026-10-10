# rust-bio-lab

20 Rust exercises in 10 units. They start at "I know Python, I've seen Rust
syntax" and end at things like a deterministic parallel Ising model and
Baum–Welch training. Everything is themed around bioinformatics,
biophysics, chemistry and data science.

Every exercise is a small library crate whose tests already exist. You
make them pass. The tests were verified against reference solutions, so
if a test fails, the bug is in your code, not the test.

## Setup

```bash
# Rust 1.85+ (edition 2024). Check with:
rustc --version
# Update if needed:
rustup update stable
```

No external crates are used, so everything works offline.

## How to work

```bash
cd rust-bio-lab
cargo test -p ex01_dna_toolkit                    # run one exercise
cargo test -p ex01_dna_toolkit -- canonical       # just tests whose name contains "canonical"
cargo test -p ex01_dna_toolkit -- --nocapture     # see println! output
./progress.sh                                     # pass/fail overview of everything
```

Run commands with `-p <crate>`. A plain `cargo test` at the root builds
*every* exercise, and some exercises are **meant** not to compile until
you fix them (see below), so it will fail.

Read the doc comment at the top of each `src/lib.rs` first. It contains
the full spec, the edge-case rules and hints. Then replace the `todo!()`s.
If a function's behaviour isn't clear from its doc comment, the tests
are the spec. Read them.

Recommended tooling: VS Code / Zed / RustRover with **rust-analyzer**, and
run `cargo clippy -p <crate>` once your tests pass. Clippy teaches a lot of
idiomatic Rust.

### Exercise types

* **Fill in**: replace `todo!()` bodies (and sometimes struct fields).
  Compiles from the start; tests fail until you're done.
* **FIX IT** (ex02, ex10): the crate *doesn't compile*. Fix the errors
  without touching the tests and without spraying `.clone()`. These are
  about learning to read `rustc` errors, which are often the best
  documentation you'll get.
* A few fill-in exercises **won't compile until you've written certain
  impls or signatures**, because writing them *is* the lesson:
  ex09 (one lifetime signature), ex13 (operator overloading). ex14 uses
  *your* `Vec3` from ex13, so do ex13 first.

## Roadmap

| Unit | Exercise | What you build | Rust concepts | Difficulty |
|---|---|---|---|---|
| 1 Ownership | `ex01_dna_toolkit` | GC content, reverse complement, k-mer counting, canonical k-mers | `&str` vs `String`, borrowing slices out of inputs, `HashMap<&str,_>`, `Cow` | ★☆☆☆☆ |
| | `ex02_fix_read_qc` | FASTQ parsing & read trimming (**fix it**) | moves, `iter_mut`, borrow conflicts, `move` closures, threads | ★★☆☆☆ |
| 2 Types | `ex03_fasta_parser` | strict multi-FASTA parser with precise errors | enums with data, `Result`, `Display`, `Error`, `FromStr` | ★★☆☆☆ |
| | `ex04_translation` | codon table, 6-frame ORF finder, hydropathy | enums, `TryFrom`, exhaustive `match`, strand coordinates | ★★★☆☆ |
| 3 Traits | `ex05_online_stats` | streaming mean/variance/min-max, mergeable (Chan), quantiles | associated types, generics, tuple impls, extension traits on `[T]` | ★★★☆☆ |
| | `ex06_kmer_iterators` | 2-bit packed rolling k-mer iterator, bitwise revcomp, moving average | `impl Iterator`, borrowing structs, blanket impls, bit twiddling | ★★★☆☆ |
| 4 Data | `ex07_expression_matrix` | RNA-seq matrix: CPM, z-scores, Pearson/Spearman, top variable genes | iterator chains, closures as params, `impl Iterator + '_` | ★★★☆☆ |
| | `ex08_kmeans` | k-means++ with restarts, silhouette score | const generics `[f64; D]`, RNG state, degenerate cases | ★★★★☆ |
| 5 Lifetimes | `ex09_pdb_parser` | zero-copy PDB parser, residues, Rg, CA contacts | structs with lifetimes, two-lifetime iterators | ★★★★☆ |
| | `ex10_fix_smiles_lexer` | SMILES lexer with lifetime bugs (**fix it**) | elision rules, `'a` vs `'_`, temporaries | ★★★☆☆ |
| 6 Alignment | `ex11_global_alignment` | Needleman–Wunsch, **Hirschberg linear-space**, Levenshtein | 2-D DP, rolling rows, divide & conquer | ★★★★☆ |
| | `ex12_local_affine` | Smith–Waterman with affine gaps (Gotoh), checked vs brute force | trait generics + closures, `?Sized`, `dyn`, 3-matrix traceback | ★★★★★ |
| 7 Simulation | `ex13_vec3_geometry` | `Vec3`, bond/dihedral angles, NeRF atom placement | `std::ops` overloading, `Sum`, `Index` | ★★★☆☆ |
| | `ex14_lennard_jones_md` | LJ fluid MD: PBC, velocity Verlet, energy conservation, g(r) | path dependencies, physics-invariant testing | ★★★★☆ |
| 8 Graphs | `ex15_smiles_parser` | SMILES → molecular graph, implicit H, Hill formula, rings | hand-written parser, index graphs, BFS, union-find | ★★★★★ |
| | `ex16_phylogeny` | Newick parse/write, patristic distance, Robinson–Foulds, Neighbor-Joining | arena trees, recursive descent, `BTreeSet` | ★★★★★ |
| 9 Concurrency | `ex17_parallel_kmers` | k-mer counting 4 ways: scoped, sharded, channel pipeline, atomics | `thread::scope`, `Mutex`, `mpsc`, `Send`/`Sync`, `BinaryHeap` | ★★★★☆ |
| | `ex18_ising_parallel` | checkerboard Metropolis, **bit-identical** for any thread count | `chunks_mut`, safe shared reads, deterministic RNG streams | ★★★★★ |
| 10 Capstone | `ex19_hmm` | Viterbi, scaled forward-backward, Baum–Welch, CpG islands | numerics, validation, EM, everything | ★★★★★ |
| | `ex20_gillespie` | reaction-network DSL, Gillespie SSA, RK4 ODE limit | DSL parsing, statistics-based testing | ★★★★★ |

Difficulty ramps by unit, and also *inside* most exercises: the early
functions in a file warm you up for the last one.

## When you're stuck

1. Re-read the doc comment. The edge-case rules are spelled out.
2. Read the failing test. `assert_eq!` prints `left` (yours) vs `right` (expected).
3. Add `dbg!(...)` liberally. It prints file/line and returns its argument.
4. Write a tiny extra test of your own for the sub-step you're unsure about.
5. Only then peek at the solutions (shipped separately). Compare *after*
   your version passes too; there are many good solutions.

## Ideas for after you finish

* Re-do ex17 with `rayon` and compare timings in `--release`.
* Benchmark ex11/ex12 with `criterion`; try a SIMD-friendly
  anti-diagonal Smith–Waterman.
* Add Kabsch superposition to ex13 (needs a 3x3 SVD, or use `nalgebra`).
* Hook ex03 + ex04 + ex19 into a CLI with `clap` that reads a real genome
  FASTA and reports ORFs and CpG islands.
* Expose ex07 to Python with `pyo3` + `maturin` and call it from pandas.
