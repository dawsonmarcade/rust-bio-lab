//! # Exercise 16 — Phylogenetic trees: Newick, Neighbor-Joining, Robinson–Foulds
//!
//! **Concepts:** trees as an **arena** (`Vec<Node>` + indices) instead of
//! `Rc<RefCell<..>>`, recursive-descent parsing, recursion over indices,
//! set algebra with `BTreeSet`, an O(n³) classic algorithm.
//!
//! Python lets you build trees from objects pointing at each other. In Rust
//! that fights the borrow checker, so the idiomatic answer is an arena: all
//! nodes live in one `Vec`, and edges are `usize` indices. Parent links
//! become trivial.
//!
//! 1. `parse_newick` / `to_newick`, e.g. `((A:1,B:2)AB:3,C:4);`.
//!    Names and branch lengths are both optional on every node. Whitespace
//!    anywhere is ignored. The trailing `;` is required.
//!    `to_newick` writes lengths with Rust's default `{}` `f64` formatting
//!    (so `1.0` -> `1`, `0.25` -> `0.25`) and writes no whitespace.
//! 2. `patristic`: path length between two leaves (missing lengths = 0).
//! 3. `splits` and `robinson_foulds`: compare topologies, treating trees as
//!    **unrooted**.
//! 4. `neighbor_joining`: build a tree from a distance matrix (Saitou & Nei).
//!    For an *additive* matrix, NJ recovers the true tree exactly. The test
//!    checks that every patristic distance equals the input distance.
//!
//! ```text
//! cargo test -p ex16_phylogeny
//! ```

#![allow(unused_variables, dead_code, unused_imports, unused_mut)]

use std::collections::BTreeSet;

#[derive(Debug, Clone, PartialEq)]
pub struct Node {
    pub name: Option<String>,
    /// Length of the branch leading *to* this node from its parent.
    pub length: Option<f64>,
    pub parent: Option<usize>,
    pub children: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Tree {
    pub nodes: Vec<Node>,
    pub root: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum NewickError {
    /// 0-based byte offset into the input
    Unexpected { pos: usize, found: Option<char> },
    BadLength { pos: usize },
    MissingSemicolon,
}

pub fn parse_newick(s: &str) -> Result<Tree, NewickError> {
    todo!()
}

impl Tree {
    pub fn to_newick(&self) -> String {
        todo!()
    }

    /// Leaf names in depth-first, left-to-right order.
    /// (Unnamed leaves are skipped.)
    pub fn leaves(&self) -> Vec<&str> {
        todo!()
    }

    pub fn find(&self, name: &str) -> Option<usize> {
        self.nodes.iter().position(|n| n.name.as_deref() == Some(name))
    }

    /// Sum of branch lengths on the path between two named nodes.
    pub fn patristic(&self, a: &str, b: &str) -> Option<f64> {
        todo!()
    }

    /// Non-trivial bipartitions of the leaf set, treating the tree as unrooted.
    /// Each internal edge splits the leaves into two sides; represent the split
    /// by the side that does **not** contain the lexicographically smallest
    /// leaf name. Skip trivial splits (a side of size < 2).
    pub fn splits(&self) -> BTreeSet<BTreeSet<String>> {
        todo!()
    }

    /// Number of splits present in exactly one of the two trees. `None` if
    /// the trees have different leaf sets.
    pub fn robinson_foulds(&self, other: &Tree) -> Option<usize> {
        todo!()
    }
}

/// Neighbor-Joining.
///
/// * `Q(i,j) = (r-2) d(i,j) - S_i - S_j` over the `r` active nodes,
///   `S_i = Σ_k d(i,k)`. Join the pair with the minimal Q (ties: smallest
///   `i`, then smallest `j`, by position in the active list).
/// * Branch lengths: `δ(i,u) = d(i,j)/2 + (S_i - S_j) / (2(r-2))`,
///   `δ(j,u) = d(i,j) - δ(i,u)`.
/// * New distances: `d(u,k) = (d(i,k) + d(j,k) - d(i,j)) / 2`.
/// * When two nodes remain, connect them with branch length `d`. (Make one
///   the root and the other its child.)
///
/// Leaves get the given names; internal nodes are unnamed. Panics on fewer
/// than 2 taxa or a non-square matrix.
pub fn neighbor_joining(names: &[&str], d: &[Vec<f64>]) -> Tree {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn parse_basic() {
        let t = parse_newick("((A:1,B:2)AB:3,C:4)root;").unwrap();
        assert_eq!(t.nodes.len(), 5);
        let r = &t.nodes[t.root];
        assert_eq!(r.name.as_deref(), Some("root"));
        assert_eq!(r.children.len(), 2);
        let ab = t.find("AB").unwrap();
        assert_eq!(t.nodes[ab].length, Some(3.0));
        assert_eq!(t.nodes[t.find("A").unwrap()].parent, Some(ab));
        assert_eq!(t.leaves(), vec!["A", "B", "C"]);
    }

    #[test]
    fn roundtrip() {
        for s in [
            "((A:1,B:2):3,C:4);",
            "(A,B,(C,D)E)F;",
            "((raccoon:19.19959,bear:6.80041):0.846,((sea_lion:11.997,seal:12.003):7.52973,((monkey:100.8593,cat:47.14069):20.59201,weasel:18.87953):2.0946):3.87382,dog:25.46154);",
            "A;",
            "(,,(,));",
        ] {
            assert_eq!(parse_newick(s).unwrap().to_newick(), s);
        }
        // whitespace is ignored
        assert_eq!(parse_newick(" ( A : 0.5 ,\n B ) ; ").unwrap().to_newick(), "(A:0.5,B);");
    }

    #[test]
    fn parse_errors() {
        assert_eq!(parse_newick("(A,B)"), Err(NewickError::MissingSemicolon));
        assert_eq!(parse_newick("(A,B;"), Err(NewickError::Unexpected { pos: 4, found: Some(';') }));
        assert_eq!(parse_newick("(A:x,B);"), Err(NewickError::BadLength { pos: 3 }));
        assert_eq!(parse_newick("(A,B);C"), Err(NewickError::Unexpected { pos: 6, found: Some('C') }));
    }

    #[test]
    fn patristic_distances() {
        let t = parse_newick("((A:1,B:2):3,C:4);").unwrap();
        assert_eq!(t.patristic("A", "B"), Some(3.0));
        assert_eq!(t.patristic("A", "C"), Some(8.0));
        assert_eq!(t.patristic("C", "C"), Some(0.0));
        assert_eq!(t.patristic("A", "Z"), None);
        let t = parse_newick("((A,B:2),C:4);").unwrap();
        assert_eq!(t.patristic("A", "C"), Some(4.0));
    }

    #[test]
    fn splits_and_rf() {
        let t1 = parse_newick("((A,B),(C,D));").unwrap();
        let t2 = parse_newick("((A,C),(B,D));").unwrap();
        let t3 = parse_newick("(A,B,(C,D));").unwrap(); // same unrooted tree as t1
        let cd: BTreeSet<String> = ["C", "D"].iter().map(|s| s.to_string()).collect();
        assert_eq!(t1.splits(), BTreeSet::from([cd]));
        assert_eq!(t1.robinson_foulds(&t3), Some(0));
        assert_eq!(t1.robinson_foulds(&t2), Some(2));
        let big1 = parse_newick("(((A,B),C),(D,(E,F)));").unwrap();
        let big2 = parse_newick("(((A,C),B),(D,(E,F)));").unwrap();
        assert_eq!(big1.splits().len(), 3);
        assert_eq!(big1.robinson_foulds(&big2), Some(2));
        let other_leaves = parse_newick("((A,B),(C,X));").unwrap();
        assert_eq!(t1.robinson_foulds(&other_leaves), None);
    }

    #[test]
    fn nj_wikipedia_example() {
        // https://en.wikipedia.org/wiki/Neighbor_joining#Example
        let names = ["a", "b", "c", "d", "e"];
        let d = vec![
            vec![0.0, 5.0, 9.0, 9.0, 8.0],
            vec![5.0, 0.0, 10.0, 10.0, 9.0],
            vec![9.0, 10.0, 0.0, 8.0, 7.0],
            vec![9.0, 10.0, 8.0, 0.0, 3.0],
            vec![8.0, 9.0, 7.0, 3.0, 0.0],
        ];
        let t = neighbor_joining(&names, &d);
        for i in 0..5 {
            for j in 0..5 {
                assert!(close(t.patristic(names[i], names[j]).unwrap(), d[i][j]), "{} {}", names[i], names[j]);
            }
        }
        assert!(close(t.nodes[t.find("a").unwrap()].length.unwrap(), 2.0));
        assert!(close(t.nodes[t.find("b").unwrap()].length.unwrap(), 3.0));
        // the output is valid Newick that parses back to the same splits
        let re = parse_newick(&t.to_newick()).unwrap();
        assert_eq!(re.robinson_foulds(&t), Some(0));
    }

    #[test]
    fn nj_recovers_random_additive_trees() {
        // build a random binary tree, take its patristic matrix, run NJ, compare
        let mut seed = 12345u64;
        let mut rnd = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 33) as f64 / (1u64 << 31) as f64
        };
        for trial in 0..20 {
            let n = 4 + trial % 9;
            let mut subtrees: Vec<String> = (0..n).map(|i| format!("t{i}:{:.3}", 0.1 + rnd())).collect();
            while subtrees.len() > 1 {
                let i = (rnd() * subtrees.len() as f64) as usize % subtrees.len();
                let a = subtrees.swap_remove(i);
                let j = (rnd() * subtrees.len() as f64) as usize % subtrees.len();
                let b = subtrees.swap_remove(j);
                subtrees.push(format!("({a},{b}):{:.3}", 0.1 + rnd()));
            }
            let truth = parse_newick(&format!("{};", subtrees[0])).unwrap();
            let names: Vec<String> = (0..n).map(|i| format!("t{i}")).collect();
            let refs: Vec<&str> = names.iter().map(String::as_str).collect();
            let d: Vec<Vec<f64>> = refs
                .iter()
                .map(|a| refs.iter().map(|b| truth.patristic(a, b).unwrap()).collect())
                .collect();
            let nj = neighbor_joining(&refs, &d);
            assert_eq!(nj.robinson_foulds(&truth), Some(0), "trial {trial}: {}", truth.to_newick());
            for a in &refs {
                for b in &refs {
                    assert!((nj.patristic(a, b).unwrap() - truth.patristic(a, b).unwrap()).abs() < 1e-6);
                }
            }
        }
    }
}
