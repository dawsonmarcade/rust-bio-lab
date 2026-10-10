//! # Exercise 06 — Bit-packed k-mers and custom iterators
//!
//! **Concepts:** implementing `Iterator` by hand, structs that borrow
//! (`KmerIter<'a>`), extension traits with a blanket impl, bit manipulation
//! on `u64`, `VecDeque`.
//!
//! Real k-mer counters (Jellyfish, KMC) pack each base into 2 bits
//! (`A=0, C=1, G=2, T=3`) so a 31-mer fits in a `u64`. A window that slides
//! one base shifts the code left by 2, ORs in the new base, and masks. That
//! makes it O(1) per base instead of O(k).
//!
//! * `encode_kmer` / `decode_kmer`: straightforward.
//! * `KmerIter`: yields `(start_position, code)` for every valid window, in
//!   **O(1) amortised per input base**. Windows containing anything other
//!   than `ACGT`/`acgt` are skipped, and the rolling state must *reset* when
//!   it hits an invalid base. The test compares you against a naive
//!   implementation on thousands of random sequences.
//! * `revcomp_code` / `canonical_code`: do it with bit operations, not by
//!   decoding to a string. Hint: complementing a 2-bit base is `3 - x`,
//!   i.e. `x ^ 0b11`.
//! * `MovingAverage`: an iterator adaptor usable on *any*
//!   `Iterator<Item = f64>` via `.moving_average(w)`.
//!
//! ```text
//! cargo test -p ex06_kmer_iterators
//! ```

#![allow(unused_variables, dead_code, unused_imports, unused_mut)]

use std::collections::VecDeque;

fn base_code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Encode up to 32 bases. First base ends up in the most significant
/// position: `"AC"` -> `0b0001`, `"CA"` -> `0b0100`.
/// `None` if `kmer.len() > 32` or there is a non-ACGT base.
pub fn encode_kmer(kmer: &[u8]) -> Option<u64> {
    todo!()
}

/// Inverse of `encode_kmer` (uppercase output).
pub fn decode_kmer(code: u64, k: usize) -> String {
    todo!()
}

/// Reverse complement of an encoded k-mer, using bit operations only.
pub fn revcomp_code(code: u64, k: usize) -> u64 {
    todo!()
}

pub fn canonical_code(code: u64, k: usize) -> u64 {
    todo!()
}

pub struct KmerIter<'a> {
    seq: &'a [u8],
    // add the rest of the state you need
}

/// Panics with the message `"k must be in 1..=32"` if `k == 0` or `k > 32`.
pub fn kmers(seq: &[u8], k: usize) -> KmerIter<'_> {
    todo!()
}

impl<'a> Iterator for KmerIter<'a> {
    /// `(start index of the window, 2-bit code)`
    type Item = (usize, u64);

    fn next(&mut self) -> Option<Self::Item> {
        todo!()
    }
}

pub struct MovingAverage<I> {
    inner: I,
    // add the rest
}

impl<I: Iterator<Item = f64>> Iterator for MovingAverage<I> {
    type Item = f64;
    /// Yields the mean of each *full* window: n inputs -> n - w + 1 outputs
    /// (none if n < w). Keep a running sum: O(1) per item.
    fn next(&mut self) -> Option<f64> {
        todo!()
    }
}

/// Extension trait: makes `.moving_average(w)` available on every
/// `Iterator<Item = f64>`. You need a *blanket impl*.
pub trait MovingAverageExt: Iterator<Item = f64> + Sized {
    /// Panics if `window == 0`.
    fn moving_average(self, window: usize) -> MovingAverage<Self>;
}

// TODO: blanket impl of MovingAverageExt for all suitable iterators
impl<I: Iterator<Item = f64>> MovingAverageExt for I {
    fn moving_average(self, window: usize) -> MovingAverage<Self> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            self.0 >> 33
        }
    }

    fn naive(seq: &[u8], k: usize) -> Vec<(usize, u64)> {
        if seq.len() < k {
            return vec![];
        }
        (0..=seq.len() - k)
            .filter_map(|i| encode_kmer(&seq[i..i + k]).map(|c| (i, c)))
            .collect()
    }

    #[test]
    fn encode_decode() {
        assert_eq!(encode_kmer(b"AC"), Some(0b0001));
        assert_eq!(encode_kmer(b"CA"), Some(0b0100));
        assert_eq!(encode_kmer(b"tgca"), Some(0b11_10_01_00));
        assert_eq!(encode_kmer(b"ANA"), None);
        assert_eq!(encode_kmer(&[b'A'; 33]), None);
        assert_eq!(encode_kmer(&[b'T'; 32]), Some(u64::MAX));
        assert_eq!(decode_kmer(0b11_10_01_00, 4), "TGCA");
        assert_eq!(decode_kmer(0, 3), "AAA");
    }

    #[test]
    fn revcomp_bits() {
        for s in ["ACG", "AAAAC", "GATTACA", "T", "ACGTACGTACGTACGTACGTACGTACGTACGT"] {
            let k = s.len();
            let rc: String = s
                .chars()
                .rev()
                .map(|c| match c { 'A' => 'T', 'T' => 'A', 'C' => 'G', _ => 'C' })
                .collect();
            let code = encode_kmer(s.as_bytes()).unwrap();
            assert_eq!(decode_kmer(revcomp_code(code, k), k), rc, "{s}");
            assert_eq!(canonical_code(code, k), code.min(encode_kmer(rc.as_bytes()).unwrap()));
        }
    }

    #[test]
    fn iter_small() {
        let got: Vec<_> = kmers(b"ACGNTTa", 2).collect();
        assert_eq!(got, vec![(0, 0b0001), (1, 0b0110), (4, 0b1111), (5, 0b1100)]);
        assert_eq!(kmers(b"AC", 3).count(), 0);
        assert_eq!(kmers(b"", 1).count(), 0);
    }

    #[test]
    fn iter_matches_naive_random() {
        let mut rng = Lcg(42);
        let alphabet = b"ACGTACGTACGTacgtN";
        for _ in 0..2000 {
            let len = (rng.next() % 80) as usize;
            let seq: Vec<u8> = (0..len).map(|_| alphabet[(rng.next() % 17) as usize]).collect();
            let k = 1 + (rng.next() % 32) as usize;
            assert_eq!(kmers(&seq, k).collect::<Vec<_>>(), naive(&seq, k), "k={k}");
        }
    }

    #[test]
    fn iter_is_lazy_and_borrows() {
        let seq = vec![b'A'; 1_000_000];
        // borrows seq, does not copy it; and only does a little work
        let first: Vec<_> = kmers(&seq, 31).take(3).collect();
        assert_eq!(first, vec![(0, 0), (1, 0), (2, 0)]);
    }

    #[test]
    #[should_panic(expected = "k must be in 1..=32")]
    fn k_zero_panics() {
        let _ = kmers(b"ACGT", 0);
    }

    #[test]
    fn moving_average() {
        let v: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0, 5.0].into_iter().moving_average(2).collect();
        assert_eq!(v, vec![1.5, 2.5, 3.5, 4.5]);
        let v: Vec<f64> = (1..=3).map(f64::from).moving_average(3).collect();
        assert_eq!(v, vec![2.0]);
        assert_eq!((1..=2).map(f64::from).moving_average(3).count(), 0);
        // chains with other adaptors (here: GC fraction in sliding windows)
        let gc = b"GGGGAAAA"
            .iter()
            .map(|&b| if b == b'G' || b == b'C' { 1.0 } else { 0.0 })
            .moving_average(4)
            .collect::<Vec<_>>();
        assert_eq!(gc, vec![1.0, 0.75, 0.5, 0.25, 0.0]);
    }
}
