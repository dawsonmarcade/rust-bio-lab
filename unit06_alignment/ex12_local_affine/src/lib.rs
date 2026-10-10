//! # Exercise 12 — Smith–Waterman with affine gaps (Gotoh)
//!
//! **Concepts:** generics over a scoring *trait* (with closures implementing
//! it via a blanket impl), `?Sized` bounds, three coupled DP matrices,
//! state-machine traceback, linear-space score pass.
//!
//! Biological gaps are bursty: one 10-base indel is far more likely than 10
//! separate 1-base indels. Affine gap costs model this: a gap of length `L`
//! scores `gap_open + (L - 1) * gap_extend` (both negative).
//!
//! The Gotoh recurrences use three matrices:
//! * `M[i][j]`: best local alignment ending with `a[i-1]` aligned to `b[j-1]`
//! * `X[i][j]`: ... ending with `a[i-1]` aligned to a gap
//! * `Y[i][j]`: ... ending with `b[j-1]` aligned to a gap
//!
//! and `H = max(0, M, X, Y)`. Work out the recurrences yourself: when does a
//! gap *open* (from `H`) and when does it *extend* (from `X`/`Y`)?
//! The local alignment starts wherever `H` hits its maximum. Pick the
//! **first** maximal cell in row-major order (smallest `i`, then smallest
//! `j`), and trace back until `H` reaches 0.
//!
//! The tests compare you against an exhaustive brute-force search over all
//! substring pairs and all alignments on thousands of tiny random inputs,
//! so a subtle recurrence bug will not survive.
//!
//! ```text
//! cargo test -p ex12_local_affine
//! ```

#![allow(unused_variables, dead_code, unused_imports, unused_mut)]

pub trait Substitution {
    fn score(&self, x: u8, y: u8) -> i32;
}

pub struct MatchMismatch {
    pub match_score: i32,
    pub mismatch: i32,
}

/// Transitions (A<->G, C<->T) are much more common than transversions.
pub struct TransitionTransversion {
    pub match_score: i32,
    pub transition: i32,
    pub transversion: i32,
}

impl Substitution for MatchMismatch {
    fn score(&self, x: u8, y: u8) -> i32 {
        todo!()
    }
}

impl Substitution for TransitionTransversion {
    fn score(&self, x: u8, y: u8) -> i32 {
        todo!()
    }
}

/// Any closure `Fn(u8, u8) -> i32` is a substitution scheme.
/// (Question for you: why can't you *also* write
/// `impl<T: SomeOtherTrait> Substitution for T`? Look up "coherence".)
impl<F: Fn(u8, u8) -> i32> Substitution for F {
    fn score(&self, x: u8, y: u8) -> i32 {
        todo!()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Gaps {
    pub open: i32,
    pub extend: i32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocalAlignment {
    pub score: i32,
    /// half-open ranges into the inputs
    pub a_range: (usize, usize),
    pub b_range: (usize, usize),
    pub a_aligned: String,
    pub b_aligned: String,
}

/// `None` if the best score is 0 (nothing aligns positively).
pub fn smith_waterman<S: Substitution + ?Sized>(
    a: &[u8],
    b: &[u8],
    sub: &S,
    gaps: Gaps,
) -> Option<LocalAlignment> {
    todo!()
}

/// Same best score as `smith_waterman`, but O(len(b)) memory.
pub fn smith_waterman_score<S: Substitution + ?Sized>(a: &[u8], b: &[u8], sub: &S, gaps: Gaps) -> i32 {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    const MM: MatchMismatch = MatchMismatch { match_score: 2, mismatch: -1 };
    const G: Gaps = Gaps { open: -3, extend: -1 };

    /// Rescore an alignment with affine gaps.
    fn rescore<S: Substitution + ?Sized>(al: &LocalAlignment, sub: &S, g: Gaps) -> i32 {
        let (a, b) = (al.a_aligned.as_bytes(), al.b_aligned.as_bytes());
        let mut sc = 0;
        let mut prev = 0u8; // 0 none, 1 gap-in-b, 2 gap-in-a
        for (&x, &y) in a.iter().zip(b) {
            let cur = if y == b'-' { 1 } else if x == b'-' { 2 } else { 0 };
            sc += match cur {
                0 => sub.score(x, y),
                c if c == prev => g.extend,
                _ => g.open,
            };
            prev = cur;
        }
        sc
    }

    fn validate<S: Substitution + ?Sized>(al: &LocalAlignment, a: &[u8], b: &[u8], sub: &S, g: Gaps) {
        assert_eq!(al.a_aligned.len(), al.b_aligned.len());
        let strip = |x: &str| x.bytes().filter(|&c| c != b'-').collect::<Vec<_>>();
        assert_eq!(strip(&al.a_aligned), &a[al.a_range.0..al.a_range.1], "{al:?}");
        assert_eq!(strip(&al.b_aligned), &b[al.b_range.0..al.b_range.1], "{al:?}");
        assert_eq!(rescore(al, sub, g), al.score, "{al:?}");
    }

    /// Exhaustive: best global affine score over every pair of substrings.
    fn brute<S: Substitution + ?Sized>(a: &[u8], b: &[u8], sub: &S, g: Gaps) -> i32 {
        fn global<S: Substitution + ?Sized>(a: &[u8], b: &[u8], sub: &S, g: Gaps, i: usize, j: usize, prev: u8) -> i32 {
            if i == a.len() && j == b.len() {
                return 0;
            }
            let mut best = i32::MIN / 4;
            if i < a.len() && j < b.len() {
                best = best.max(sub.score(a[i], b[j]) + global(a, b, sub, g, i + 1, j + 1, 0));
            }
            if i < a.len() {
                let c = if prev == 1 { g.extend } else { g.open };
                best = best.max(c + global(a, b, sub, g, i + 1, j, 1));
            }
            if j < b.len() {
                let c = if prev == 2 { g.extend } else { g.open };
                best = best.max(c + global(a, b, sub, g, i, j + 1, 2));
            }
            best
        }
        let mut best = 0;
        for i0 in 0..a.len() {
            for i1 in i0 + 1..=a.len() {
                for j0 in 0..b.len() {
                    for j1 in j0 + 1..=b.len() {
                        best = best.max(global(&a[i0..i1], &b[j0..j1], sub, g, 0, 0, 0));
                    }
                }
            }
        }
        best
    }

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn seq(&mut self, max: u64) -> Vec<u8> {
            let n = self.next() % (max + 1);
            (0..n).map(|_| b"ACGT"[(self.next() % 4) as usize]).collect()
        }
    }

    #[test]
    fn substitution_impls() {
        assert_eq!(MM.score(b'A', b'A'), 2);
        let tt = TransitionTransversion { match_score: 1, transition: -1, transversion: -2 };
        assert_eq!(tt.score(b'A', b'G'), -1);
        assert_eq!(tt.score(b'C', b'T'), -1);
        assert_eq!(tt.score(b'A', b'C'), -2);
        let closure = |x: u8, y: u8| if x == y { 5 } else { -4 };
        assert_eq!(closure.score(b'G', b'G'), 5);
        // trait objects work too (that's what ?Sized is for)
        let dynsub: &dyn Substitution = &tt;
        assert!(smith_waterman(b"ACGT", b"ACGT", dynsub, G).is_some());
    }

    #[test]
    fn simple_core_match() {
        let al = smith_waterman(b"TTACGTAA", b"GGACGTCC", &MM, G).unwrap();
        assert_eq!(al.score, 8);
        assert_eq!((al.a_range, al.b_range), ((2, 6), (2, 6)));
        assert_eq!(al.a_aligned, "ACGT");
    }

    #[test]
    fn nothing_aligns() {
        assert_eq!(smith_waterman(b"AAAA", b"CCCC", &MM, G), None);
        assert_eq!(smith_waterman(b"", b"ACGT", &MM, G), None);
        assert_eq!(smith_waterman_score(b"AAAA", b"CCCC", &MM, G), 0);
    }

    #[test]
    fn one_long_gap_beats_many() {
        let half = b"ACGTTGCAACGT";
        let a = [&half[..], b"GGGGGG", &half[..]].concat();
        let b = [&half[..], &half[..]].concat();
        let sub = MatchMismatch { match_score: 2, mismatch: -3 };
        let g = Gaps { open: -5, extend: -1 };
        let al = smith_waterman(&a, &b, &sub, g).unwrap();
        validate(&al, &a, &b, &sub, g);
        assert_eq!(al.score, 48 - 5 - 5);
        assert!(al.b_aligned.contains("------"), "{al:?}");
    }

    #[test]
    fn matches_brute_force() {
        let mut rng = Lcg(2024);
        let schemes: Vec<(Box<dyn Substitution>, Gaps)> = vec![
            (Box::new(MatchMismatch { match_score: 2, mismatch: -1 }), Gaps { open: -3, extend: -1 }),
            (Box::new(MatchMismatch { match_score: 3, mismatch: -2 }), Gaps { open: -2, extend: -2 }),
            (
                Box::new(TransitionTransversion { match_score: 4, transition: -1, transversion: -3 }),
                Gaps { open: -4, extend: -1 },
            ),
            (Box::new(|x: u8, y: u8| if x == y { 1 } else { -1 }), Gaps { open: -1, extend: 0 }),
        ];
        for _ in 0..600 {
            let a = rng.seq(5);
            let b = rng.seq(5);
            for (sub, g) in &schemes {
                let expect = brute(&a, &b, sub.as_ref(), *g);
                let got = smith_waterman(&a, &b, sub.as_ref(), *g);
                assert_eq!(got.as_ref().map_or(0, |x| x.score), expect, "a={a:?} b={b:?}");
                if let Some(al) = &got {
                    validate(al, &a, &b, sub.as_ref(), *g);
                }
                assert_eq!(smith_waterman_score(&a, &b, sub.as_ref(), *g), expect);
            }
        }
    }

    #[test]
    fn linear_space_agrees_on_bigger_inputs() {
        let mut rng = Lcg(5);
        for _ in 0..200 {
            let a = rng.seq(120);
            let b = rng.seq(120);
            let full = smith_waterman(&a, &b, &MM, G);
            if let Some(al) = &full {
                validate(al, &a, &b, &MM, G);
            }
            assert_eq!(full.map_or(0, |x| x.score), smith_waterman_score(&a, &b, &MM, G));
        }
    }

    #[test]
    fn first_maximum_wins() {
        // "ACG" occurs twice in a; the first (row-major smallest i) end cell wins
        let al = smith_waterman(b"ACGTTTTACG", b"ACG", &MM, G).unwrap();
        assert_eq!(al.a_range, (0, 3));
    }
}
