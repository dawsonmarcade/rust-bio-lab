//! # Exercise 01 — DNA toolkit
//!
//! **Concepts:** `&str` vs `String`, borrowing slices out of an input,
//! `HashMap<&str, _>` keyed by borrowed data, `Cow`, `Option`/`Result`.
//!
//! The warm-up. Every function takes borrowed input. Several of them return
//! data that *borrows from* that input instead of copying it. Pay attention
//! to which ones can do that, and why `kmer_counts` can return
//! `HashMap<&str, usize>` without a single allocation per k-mer.
//!
//! Replace every `todo!()` and run:
//!
//! ```text
//! cargo test -p ex01_dna_toolkit
//! ```
//!
//! Assume input is ASCII unless a function says otherwise.

#![allow(unused_variables, dead_code, unused_imports)]

use std::borrow::Cow;
use std::collections::HashMap;

/// Fraction of G or C among the *unambiguous* bases (A, C, G, T),
/// case-insensitive. Any other character (N, gaps, ...) is ignored entirely.
/// Returns `0.0` if there are no unambiguous bases.
pub fn gc_content(seq: &str) -> f64 {
    todo!()
}

/// Complement of a single base. A<->T, C<->G, N->N. Case is preserved
/// (`'a'` -> `'t'`). Anything else -> `None`.
pub fn complement(base: char) -> Option<char> {
    todo!()
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InvalidBase {
    /// Character index into the *original* (not reversed) sequence.
    pub position: usize,
    pub found: char,
}

/// Reverse complement. Fails on the *first* (lowest-index) invalid character.
pub fn reverse_complement(seq: &str) -> Result<String, InvalidBase> {
    todo!()
}

/// Hamming distance (number of differing characters). `None` if the strings
/// have different lengths in *characters*.
pub fn hamming(a: &str, b: &str) -> Option<usize> {
    todo!()
}

/// Start positions of every (possibly overlapping) occurrence of `motif`.
/// An empty motif matches nowhere.
pub fn find_motif(seq: &str, motif: &str) -> Vec<usize> {
    todo!()
}

/// Count every length-`k` window of `seq`. The keys must be *slices of the
/// input* (no allocation per k-mer). Only uppercase `ACGT` are valid: skip
/// any window containing another character. `k == 0` or `k > len` gives an
/// empty map.
pub fn kmer_counts(seq: &str, k: usize) -> HashMap<&str, usize> {
    todo!()
}

/// All k-mers tied for the highest count, sorted lexicographically.
/// Empty vec if there are none.
pub fn most_frequent_kmers(seq: &str, k: usize) -> Vec<&str> {
    todo!()
}

/// The *canonical* form of a k-mer: the lexicographically smaller of the
/// k-mer and its reverse complement. Input is uppercase ACGT.
///
/// Return `Cow::Borrowed` when the input already *is* canonical (no
/// allocation), `Cow::Owned` otherwise.
pub fn canonical(kmer: &str) -> Cow<'_, str> {
    todo!()
}

/// The longest run of a single repeated character, as `(char, start, len)`.
/// On ties, the earliest run wins. `None` for an empty string.
/// Characters are compared exactly (case-sensitive).
pub fn longest_homopolymer(seq: &str) -> Option<(char, usize, usize)> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gc_basic() {
        assert!((gc_content("GGCC") - 1.0).abs() < 1e-12);
        assert!((gc_content("ATGC") - 0.5).abs() < 1e-12);
        assert!((gc_content("atgcNNNN--") - 0.5).abs() < 1e-12);
        assert_eq!(gc_content(""), 0.0);
        assert_eq!(gc_content("NNN"), 0.0);
    }

    #[test]
    fn complement_cases() {
        assert_eq!(complement('A'), Some('T'));
        assert_eq!(complement('g'), Some('c'));
        assert_eq!(complement('N'), Some('N'));
        assert_eq!(complement('X'), None);
    }

    #[test]
    fn revcomp() {
        assert_eq!(reverse_complement("ATGCN").unwrap(), "NGCAT");
        assert_eq!(reverse_complement("aaCC").unwrap(), "GGtt");
        assert_eq!(reverse_complement("").unwrap(), "");
        assert_eq!(
            reverse_complement("ACXGZ"),
            Err(InvalidBase { position: 2, found: 'X' })
        );
    }

    #[test]
    fn hamming_cases() {
        assert_eq!(hamming("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT"), Some(7));
        assert_eq!(hamming("AAA", "AA"), None);
        assert_eq!(hamming("", ""), Some(0));
    }

    #[test]
    fn motif_overlapping() {
        assert_eq!(find_motif("GATATATGCATATACTT", "ATAT"), vec![1, 3, 9]);
        assert_eq!(find_motif("AAAA", "AA"), vec![0, 1, 2]);
        assert!(find_motif("ACGT", "").is_empty());
        assert!(find_motif("AC", "ACGT").is_empty());
    }

    #[test]
    fn kmers_borrow_and_skip() {
        let seq = String::from("ACGTNACGT");
        let counts = kmer_counts(&seq, 3);
        assert_eq!(counts.get("ACG"), Some(&2));
        assert_eq!(counts.get("CGT"), Some(&2));
        assert_eq!(counts.len(), 2, "windows with N must be skipped: {counts:?}");
        // keys really are slices into `seq`
        for k in counts.keys() {
            let p = k.as_ptr() as usize;
            let base = seq.as_ptr() as usize;
            assert!(p >= base && p < base + seq.len());
        }
        assert!(kmer_counts("ACGT", 0).is_empty());
        assert!(kmer_counts("ACGT", 5).is_empty());
        assert!(kmer_counts("acgt", 2).is_empty());
    }

    #[test]
    fn most_frequent() {
        let v = most_frequent_kmers("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4);
        assert_eq!(v, vec!["CATG", "GCAT"]);
        assert!(most_frequent_kmers("NNNN", 2).is_empty());
    }

    #[test]
    fn canonical_cow() {
        assert!(matches!(canonical("ACG"), Cow::Borrowed("ACG")));
        match canonical("TTT") {
            Cow::Owned(s) => assert_eq!(s, "AAA"),
            other => panic!("expected Owned, got {other:?}"),
        }
        // palindrome: borrowed
        assert!(matches!(canonical("ACGT"), Cow::Borrowed(_)));
    }

    #[test]
    fn homopolymer() {
        assert_eq!(longest_homopolymer("ACCCGTTTA"), Some(('C', 1, 3)));
        assert_eq!(longest_homopolymer("AAGGGGT"), Some(('G', 2, 4)));
        assert_eq!(longest_homopolymer("A"), Some(('A', 0, 1)));
        assert_eq!(longest_homopolymer(""), None);
    }
}
