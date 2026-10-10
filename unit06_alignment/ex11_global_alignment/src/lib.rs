//! # Exercise 11 — Global alignment: Needleman–Wunsch and Hirschberg
//!
//! **Concepts:** 2-D dynamic programming with flat `Vec`s, traceback,
//! memory-efficient DP with rolling rows, divide-and-conquer recursion
//! over slices.
//!
//! 1. `needleman_wunsch`: classic O(nm) time *and* space, with traceback and
//!    a fixed tie-breaking order so the output is deterministic.
//! 2. `nw_score_last_row`: O(min) space, returns the whole last DP row.
//! 3. `hirschberg`: an **optimal** alignment in **linear space**. This is the
//!    hard one: split `a` in half, run the score pass forwards on the top
//!    half and *backwards* (on reversed sequences) on the bottom half, find
//!    the best split point of `b`, and recurse. Base cases matter.
//! 4. `levenshtein`: unit-cost edit distance, two rows of memory.
//!
//! Scoring is linear: every gap character costs `gap` (negative).
//!
//! ```text
//! cargo test -p ex11_global_alignment
//! ```

#![allow(unused_variables, dead_code, unused_imports, unused_mut)]

#[derive(Debug, Clone, Copy)]
pub struct Scoring {
    pub match_score: i32,
    pub mismatch: i32,
    pub gap: i32,
}

impl Scoring {
    pub fn sub(&self, x: u8, y: u8) -> i32 {
        if x == y { self.match_score } else { self.mismatch }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub score: i32,
    /// `a` with `-` inserted for gaps
    pub a: String,
    /// `b` with `-` inserted for gaps; same length as `a`
    pub b: String,
}

/// Full-matrix Needleman–Wunsch.
///
/// Traceback starts at (n, m). When several moves achieve the optimum,
/// prefer **diagonal**, then **up** (consume a char of `a`, gap in `b`),
/// then **left** (gap in `a`).
pub fn needleman_wunsch(a: &[u8], b: &[u8], s: &Scoring) -> Alignment {
    todo!()
}

/// Last row of the NW matrix (`b.len() + 1` values) using O(len(b)) memory.
/// `row[j]` = best score aligning all of `a` with `b[..j]`.
pub fn nw_score_last_row(a: &[u8], b: &[u8], s: &Scoring) -> Vec<i32> {
    todo!()
}

/// An optimal global alignment in O(len(a) + len(b)) memory (excluding the
/// output). The score must equal Needleman–Wunsch's; the exact alignment may
/// differ when there are ties.
pub fn hirschberg(a: &[u8], b: &[u8], s: &Scoring) -> Alignment {
    todo!()
}

/// Levenshtein distance with two rows of memory.
pub fn levenshtein(a: &[u8], b: &[u8]) -> usize {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    const S: Scoring = Scoring { match_score: 1, mismatch: -1, gap: -2 };

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn seq(&mut self, max_len: u64) -> Vec<u8> {
            let n = self.next() % (max_len + 1);
            (0..n).map(|_| b"ACGT"[(self.next() % 4) as usize]).collect()
        }
    }

    /// Independent check: is `al` a valid alignment of a and b, with the claimed score?
    fn check(al: &Alignment, a: &[u8], b: &[u8], s: &Scoring) {
        assert_eq!(al.a.len(), al.b.len(), "{al:?}");
        let strip = |x: &str| x.bytes().filter(|&c| c != b'-').collect::<Vec<_>>();
        assert_eq!(strip(&al.a), a);
        assert_eq!(strip(&al.b), b);
        let mut sc = 0;
        for (x, y) in al.a.bytes().zip(al.b.bytes()) {
            assert!(!(x == b'-' && y == b'-'), "gap aligned to gap: {al:?}");
            sc += if x == b'-' || y == b'-' { s.gap } else { s.sub(x, y) };
        }
        assert_eq!(sc, al.score, "{al:?}");
    }

    #[test]
    fn nw_textbook() {
        let al = needleman_wunsch(b"GATTACA", b"GCATGCU", &S);
        check(&al, b"GATTACA", b"GCATGCU", &S);
        assert_eq!(al.score, -1);
    }

    #[test]
    fn nw_tie_breaking_is_deterministic() {
        let al = needleman_wunsch(b"AC", b"A", &S);
        assert_eq!((al.a.as_str(), al.b.as_str(), al.score), ("AC", "A-", -1));
        let al = needleman_wunsch(b"A", b"TA", &S);
        assert_eq!((al.a.as_str(), al.b.as_str()), ("-A", "TA"));
        // two mismatches (-2) lose to gap + match + gap (-1)
        let al = needleman_wunsch(b"AT", b"TA", &Scoring { match_score: 1, mismatch: -1, gap: -1 });
        assert_eq!((al.a.as_str(), al.b.as_str(), al.score), ("-AT", "TA-", -1));
    }

    #[test]
    fn empty_inputs() {
        let al = needleman_wunsch(b"", b"ACG", &S);
        assert_eq!((al.a.as_str(), al.b.as_str(), al.score), ("---", "ACG", -6));
        let al = hirschberg(b"ACG", b"", &S);
        assert_eq!((al.a.as_str(), al.b.as_str(), al.score), ("ACG", "---", -6));
        assert_eq!(hirschberg(b"", b"", &S).score, 0);
    }

    #[test]
    fn last_row() {
        assert_eq!(nw_score_last_row(b"", b"AC", &S), vec![0, -2, -4]);
        let full = needleman_wunsch(b"GATTACA", b"GCATGCU", &S).score;
        assert_eq!(*nw_score_last_row(b"GATTACA", b"GCATGCU", &S).last().unwrap(), full);
    }

    #[test]
    fn hirschberg_matches_nw_score_random() {
        let mut rng = Lcg(7);
        let schemes = [
            S,
            Scoring { match_score: 2, mismatch: -3, gap: -1 },
            Scoring { match_score: 0, mismatch: -1, gap: -1 },
        ];
        for _ in 0..1500 {
            let a = rng.seq(40);
            let b = rng.seq(40);
            for s in &schemes {
                let nw = needleman_wunsch(&a, &b, s);
                check(&nw, &a, &b, s);
                let h = hirschberg(&a, &b, s);
                check(&h, &a, &b, s);
                assert_eq!(h.score, nw.score, "a={:?} b={:?}", a, b);
            }
        }
    }

    #[test]
    fn hirschberg_long_sequences() {
        // 3000 x 3000 full matrix would be 36 MB of i32; Hirschberg doesn't care.
        let mut rng = Lcg(99);
        let a: Vec<u8> = (0..3000).map(|_| b"ACGT"[(rng.next() % 4) as usize]).collect();
        let mut b = a.clone();
        for i in (0..b.len()).step_by(50) {
            b[i] = b'T';
        }
        b.drain(1000..1010);
        let h = hirschberg(&a, &b, &S);
        check(&h, &a, &b, &S);
        assert_eq!(h.score, *nw_score_last_row(&a, &b, &S).last().unwrap());
    }

    #[test]
    fn edit_distance() {
        assert_eq!(levenshtein(b"kitten", b"sitting"), 3);
        assert_eq!(levenshtein(b"", b"abc"), 3);
        assert_eq!(levenshtein(b"PLEASANTLY", b"MEANLY"), 5);
        assert_eq!(levenshtein(b"same", b"same"), 0);
    }
}
