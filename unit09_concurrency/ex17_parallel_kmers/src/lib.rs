//! # Exercise 17 — Parallel k-mer counting, four ways
//!
//! **Concepts:** `std::thread::scope` (borrowing data across threads
//! without `Arc`), `Mutex`, sharding to reduce contention, `mpsc` channels
//! with backpressure (`sync_channel`), `Arc<Mutex<Receiver>>` work queues,
//! atomics, `BinaryHeap` for top-N. And the `Send`/`Sync` errors you'll hit
//! along the way.
//!
//! Every function must return **exactly** the same counts as
//! `count_serial` for any thread count. The tests run them all against
//! each other.
//!
//! 1. `count_scoped`: split the reads into `n_threads` chunks, give each
//!    thread its own `HashMap`, merge at the end. The easy, fast baseline.
//! 2. `count_sharded`: one shared table split into `n_shards`
//!    `Mutex<HashMap>`s by `code % n_shards`. Lock contention kills
//!    performance if you lock per k-mer, so buffer k-mers per shard locally
//!    and flush in batches.
//! 3. `count_pipeline`: the reads arrive as an `Iterator` (think: decompressing
//!    a FASTQ.gz stream). The calling thread reads batches and sends them over a
//!    **bounded** channel; `n_workers` workers pull batches from a shared
//!    receiver and send their partial maps back on a second channel.
//! 4. `total_kmers_atomic`: count all valid k-mers with an `AtomicU64`.
//!
//! Finally `top_n` with a min-heap of size `n`: O(m log n), not a full sort.
//!
//! ```text
//! cargo test -p ex17_parallel_kmers
//! cargo test -p ex17_parallel_kmers --release -- --nocapture bench   # timings
//! ```

#![allow(unused_variables, dead_code, unused_imports, unused_mut)]

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc;
use std::sync::{Arc, Mutex};
use std::thread;

pub type Counts = HashMap<u64, u32>;

/// Canonical 2-bit k-mer codes of a read (provided; it's exercise 06's job).
pub fn canonical_kmers(seq: &[u8], k: usize) -> impl Iterator<Item = u64> + '_ {
    assert!((1..=32).contains(&k));
    let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let shift = 2 * (k as u32 - 1);
    let (mut fwd, mut rev, mut valid) = (0u64, 0u64, 0usize);
    seq.iter().filter_map(move |&b| {
        let c = match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => {
                valid = 0;
                return None;
            }
        };
        fwd = ((fwd << 2) | c) & mask;
        rev = (rev >> 2) | ((3 - c) << shift);
        valid += 1;
        (valid >= k).then(|| fwd.min(rev))
    })
}

pub fn count_serial(reads: &[Vec<u8>], k: usize) -> Counts {
    let mut m = Counts::new();
    for r in reads {
        for c in canonical_kmers(r, k) {
            *m.entry(c).or_insert(0) += 1;
        }
    }
    m
}

/// Panics if `n_threads == 0`.
pub fn count_scoped(reads: &[Vec<u8>], k: usize, n_threads: usize) -> Counts {
    todo!()
}

/// Panics if `n_threads == 0` or `n_shards == 0`.
pub fn count_sharded(reads: &[Vec<u8>], k: usize, n_threads: usize, n_shards: usize) -> Counts {
    todo!()
}

/// Pipeline. Batches of `batch_size` reads go through a
/// `sync_channel(2 * n_workers)`. Panics if `n_workers == 0` or
/// `batch_size == 0`. Note the iterator is *not* required to be `Send`:
/// it stays on the calling thread.
pub fn count_pipeline<I>(reads: I, k: usize, n_workers: usize, batch_size: usize) -> Counts
where
    I: Iterator<Item = Vec<u8>>,
{
    todo!()
}

/// Total number of valid k-mers (not distinct), using one `AtomicU64`
/// shared by `n_threads` scoped threads.
pub fn total_kmers_atomic(reads: &[Vec<u8>], k: usize, n_threads: usize) -> u64 {
    todo!()
}

/// The `n` most frequent k-mers, by count descending, ties by code ascending.
pub fn top_n(counts: &Counts, n: usize) -> Vec<(u64, u32)> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_reads(n: usize, len: usize, seed: u64) -> Vec<Vec<u8>> {
        let mut s = seed;
        let mut next = move || {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (s >> 33) as usize
        };
        // a small "genome" so k-mers repeat a lot
        let genome: Vec<u8> = (0..5000).map(|_| b"ACGT"[next() % 4]).collect();
        (0..n)
            .map(|_| {
                let start = next() % (genome.len() - len);
                let mut r = genome[start..start + len].to_vec();
                if next() % 10 == 0 {
                    r[next() % len] = b'N';
                }
                r
            })
            .collect()
    }

    #[test]
    fn canonical_helper_is_strand_symmetric() {
        let f: Vec<u64> = canonical_kmers(b"ACGTTGCA", 3).collect();
        let r: Vec<u64> = canonical_kmers(b"TGCAACGT", 3).collect();
        let mut a = f.clone();
        let mut b = r.clone();
        a.sort();
        b.sort();
        assert_eq!(a, b);
    }

    #[test]
    fn all_strategies_agree() {
        let reads = random_reads(3000, 100, 1);
        let k = 21;
        let serial = count_serial(&reads, k);
        for t in [1, 2, 3, 8, 5000] {
            assert_eq!(count_scoped(&reads, k, t), serial, "scoped t={t}");
            assert_eq!(count_sharded(&reads, k, t, 7), serial, "sharded t={t}");
            assert_eq!(count_pipeline(reads.iter().cloned(), k, t.min(16), 64), serial, "pipeline t={t}");
        }
        let expected_total: u64 = serial.values().map(|&v| v as u64).sum();
        assert_eq!(total_kmers_atomic(&reads, k, 4), expected_total);
    }

    #[test]
    fn edge_cases() {
        let empty: Vec<Vec<u8>> = vec![];
        assert!(count_scoped(&empty, 5, 4).is_empty());
        assert!(count_sharded(&empty, 5, 4, 3).is_empty());
        assert!(count_pipeline(std::iter::empty(), 5, 2, 10).is_empty());
        assert_eq!(total_kmers_atomic(&empty, 5, 2), 0);
        let short = vec![b"ACG".to_vec(), b"NNNNNNN".to_vec()];
        assert!(count_scoped(&short, 5, 2).is_empty());
        let one = vec![b"AAAAAA".to_vec()];
        assert_eq!(count_sharded(&one, 3, 4, 1), HashMap::from([(0u64, 4u32)]));
    }

    #[test]
    fn pipeline_accepts_non_send_iterator() {
        use std::rc::Rc;
        let shared = Rc::new(vec![b"ACGTACGT".to_vec(); 10]);
        let it = (0..10).map({
            let shared = Rc::clone(&shared);
            move |i| shared[i].clone()
        });
        let c = count_pipeline(it, 4, 3, 3);
        assert_eq!(c, count_serial(&shared, 4));
    }

    #[test]
    fn top_n_matches_sort() {
        let reads = random_reads(500, 80, 9);
        let counts = count_serial(&reads, 11);
        let mut all: Vec<(u64, u32)> = counts.iter().map(|(&a, &b)| (a, b)).collect();
        all.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        for n in [0, 1, 5, 50, all.len(), all.len() + 10] {
            assert_eq!(top_n(&counts, n), all[..n.min(all.len())].to_vec(), "n={n}");
        }
        let ties = HashMap::from([(5u64, 2u32), (3, 2), (9, 2), (1, 1)]);
        assert_eq!(top_n(&ties, 2), vec![(3, 2), (5, 2)]);
    }

    #[test]
    fn bench() {
        let reads = random_reads(40_000, 150, 3);
        let k = 25;
        let t0 = std::time::Instant::now();
        let s = count_serial(&reads, k);
        let ts = t0.elapsed();
        let threads = thread::available_parallelism().map_or(4, |n| n.get());
        let t0 = std::time::Instant::now();
        let p = count_scoped(&reads, k, threads);
        let tp = t0.elapsed();
        let t0 = std::time::Instant::now();
        let h = count_sharded(&reads, k, threads, 64);
        let th = t0.elapsed();
        assert_eq!(s, p);
        assert_eq!(s, h);
        println!("serial {ts:?} | scoped x{threads} {tp:?} | sharded x{threads} {th:?}");
    }
}
