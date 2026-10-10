//! # Exercise 08 — k-means++ clustering with const generics
//!
//! **Concepts:** const generics (`[f64; D]`), passing `&mut` RNG state,
//! algorithm correctness under edge cases (empty clusters, ties,
//! duplicates), struct return types.
//!
//! Clustering cells by expression profile, conformations by coordinates,
//! molecules by descriptors: k-means is everywhere. Implement:
//!
//! * `kmeans_pp_init`: k-means++ seeding (D² sampling).
//! * `kmeans`: Lloyd's algorithm from k-means++ seeds, `n_init` restarts,
//!   keep the lowest-inertia run.
//! * `silhouette`: mean silhouette score (matches `sklearn.metrics.silhouette_score`).
//!
//! Use the provided `Rng` only. The tests check behaviour, not your exact
//! random choices.
//!
//! Rules for `kmeans`:
//! * assignment: nearest centroid by squared Euclidean distance; ties go to
//!   the **lowest** centroid index.
//! * update: centroid = mean of its points. A cluster that becomes empty
//!   **keeps its previous centroid**.
//! * stop when no label changes, or after `max_iter` iterations.
//! * `inertia` = sum of squared distances from each point to its centroid.
//!
//! The hard part of this exercise is k-means++ when points are duplicated
//! and every remaining point has distance 0: you must still return `k`
//! centroids and must not divide by zero or loop forever.
//!
//! ```text
//! cargo test -p ex08_kmeans
//! ```

#![allow(unused_variables, dead_code, unused_imports, unused_mut)]

/// SplitMix64: tiny, fast, good-enough PRNG. Provided; don't change it.
#[derive(Debug, Clone)]
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed)
    }
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    /// Uniform in [0, 1).
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    /// Uniform in 0..n. Panics if n == 0.
    pub fn below(&mut self, n: usize) -> usize {
        assert!(n > 0);
        (self.next_f64() * n as f64) as usize
    }
}

pub fn sq_dist<const D: usize>(a: &[f64; D], b: &[f64; D]) -> f64 {
    todo!()
}

/// Index of the nearest centroid (ties -> lowest index) and the squared distance.
pub fn nearest<const D: usize>(p: &[f64; D], centroids: &[[f64; D]]) -> (usize, f64) {
    todo!()
}

/// k-means++ seeding. First centroid: a uniformly random point. Each next
/// centroid: a point chosen with probability proportional to its squared
/// distance to the nearest already-chosen centroid. If all those distances
/// are 0, pick uniformly at random instead.
///
/// Panics if `k == 0` or `k > data.len()`.
pub fn kmeans_pp_init<const D: usize>(data: &[[f64; D]], k: usize, rng: &mut Rng) -> Vec<[f64; D]> {
    todo!()
}

#[derive(Debug, Clone, PartialEq)]
pub struct KMeansResult<const D: usize> {
    pub labels: Vec<usize>,
    pub centroids: Vec<[f64; D]>,
    pub inertia: f64,
    /// Lloyd iterations used by the *best* run (assignment+update = 1).
    pub iterations: usize,
}

/// Run k-means `n_init` times (each with fresh k-means++ seeds from the same
/// `Rng`, seeded once with `seed`) and return the run with the lowest inertia.
pub fn kmeans<const D: usize>(
    data: &[[f64; D]],
    k: usize,
    max_iter: usize,
    n_init: usize,
    seed: u64,
) -> KMeansResult<D> {
    todo!()
}

/// Mean silhouette coefficient over all points (sklearn convention):
/// `s(i) = (b - a) / max(a, b)` using *Euclidean* (not squared) distances,
/// where `a` = mean distance to other points in the same cluster and `b` =
/// the smallest mean distance to the points of another cluster. Points in a
/// singleton cluster get `s = 0`. Returns `None` unless
/// `2 <= number of distinct labels <= n - 1`. Labels need not be contiguous.
pub fn silhouette<const D: usize>(data: &[[f64; D]], labels: &[usize]) -> Option<f64> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn blobs(seed: u64) -> (Vec<[f64; 3]>, Vec<usize>) {
        let mut rng = Rng::new(seed);
        let centers = [[0.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 10.0], [10.0, 0.0, 10.0]];
        let mut pts = vec![];
        let mut truth = vec![];
        for i in 0..200 {
            let c = i % 4;
            let p = [
                centers[c][0] + rng.next_f64() * 2.0 - 1.0,
                centers[c][1] + rng.next_f64() * 2.0 - 1.0,
                centers[c][2] + rng.next_f64() * 2.0 - 1.0,
            ];
            pts.push(p);
            truth.push(c);
        }
        (pts, truth)
    }

    /// Same partition up to relabelling?
    fn same_partition(a: &[usize], b: &[usize]) -> bool {
        use std::collections::HashMap;
        let mut ab = HashMap::new();
        let mut ba = HashMap::new();
        a.iter().zip(b).all(|(x, y)| *ab.entry(x).or_insert(y) == y && *ba.entry(y).or_insert(x) == x)
    }

    #[test]
    fn helpers() {
        assert_eq!(sq_dist(&[0.0, 0.0], &[3.0, 4.0]), 25.0);
        assert_eq!(nearest(&[1.0], &[[0.0], [2.0], [5.0]]), (0, 1.0)); // tie -> lowest
        assert_eq!(nearest(&[4.0], &[[0.0], [2.0], [5.0]]), (2, 1.0));
    }

    #[test]
    fn pp_init_picks_distinct_data_points() {
        let (pts, _) = blobs(1);
        let mut rng = Rng::new(7);
        let c = kmeans_pp_init(&pts, 4, &mut rng);
        assert_eq!(c.len(), 4);
        for x in &c {
            assert!(pts.contains(x));
        }
        // with well separated blobs, D^2 seeding almost surely hits 4 different blobs
        let blob_of = |p: &[f64; 3]| ((p[0] > 5.0) as u8) * 4 + ((p[1] > 5.0) as u8) * 2 + (p[2] > 5.0) as u8;
        let mut ids: Vec<u8> = c.iter().map(blob_of).collect();
        ids.sort();
        ids.dedup();
        assert_eq!(ids.len(), 4);
    }

    #[test]
    fn pp_init_handles_duplicates() {
        let pts = vec![[1.0, 1.0]; 5];
        let mut rng = Rng::new(3);
        let c = kmeans_pp_init(&pts, 5, &mut rng);
        assert_eq!(c, vec![[1.0, 1.0]; 5]);
        let pts = vec![[0.0], [0.0], [0.0], [9.0]];
        for s in 0..50 {
            let c = kmeans_pp_init(&pts, 2, &mut Rng::new(s));
            assert!(c.contains(&[0.0]) && c.contains(&[9.0]), "seed {s}: {c:?}");
        }
    }

    #[test]
    fn recovers_blobs() {
        let (pts, truth) = blobs(99);
        let r = kmeans(&pts, 4, 100, 5, 1234);
        assert!(same_partition(&r.labels, &truth));
        assert!(r.iterations >= 1 && r.iterations <= 100);
        // inertia equals sum of squared distances to the reported centroids
        let inertia: f64 = pts.iter().zip(&r.labels).map(|(p, &l)| sq_dist(p, &r.centroids[l])).sum();
        assert!((inertia - r.inertia).abs() < 1e-9);
        // centroids are the cluster means
        for c in 0..4 {
            let members: Vec<_> = pts.iter().zip(&r.labels).filter(|(_, l)| **l == c).map(|(p, _)| p).collect();
            for d in 0..3 {
                let m = members.iter().map(|p| p[d]).sum::<f64>() / members.len() as f64;
                assert!((m - r.centroids[c][d]).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn more_restarts_never_worse() {
        let (pts, _) = blobs(5);
        let one = kmeans(&pts, 7, 100, 1, 42);
        let many = kmeans(&pts, 7, 100, 20, 42);
        assert!(many.inertia <= one.inertia + 1e-12);
    }

    #[test]
    fn k_equals_n() {
        let pts = [[0.0], [1.0], [5.0]];
        let r = kmeans(&pts, 3, 10, 1, 0);
        assert_eq!(r.inertia, 0.0);
        let mut l = r.labels.clone();
        l.sort();
        assert_eq!(l, vec![0, 1, 2]);
    }

    #[test]
    fn silhouette_matches_sklearn() {
        let x = [[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [5.0, 5.0], [5.0, 6.0], [9.0, 9.0]];
        // sklearn.metrics.silhouette_score(X, [0,0,0,1,1,2]) == 0.6887205403884762
        let s = silhouette(&x, &[0, 0, 0, 1, 1, 2]).unwrap();
        assert!((s - 0.6887205403884762).abs() < 1e-12, "{s}");
        // labels need not be contiguous
        let s = silhouette(&x, &[7, 7, 3, 3, 40, 40]).unwrap();
        assert!((s - 0.07679668666623672).abs() < 1e-12, "{s}");
        assert_eq!(silhouette(&x, &[0; 6]), None);
        assert_eq!(silhouette(&x, &[0, 1, 2, 3, 4, 5]), None);
    }

    #[test]
    fn silhouette_prefers_true_k() {
        let (pts, _) = blobs(11);
        let s: Vec<f64> = (2..=6)
            .map(|k| silhouette(&pts, &kmeans(&pts, k, 100, 5, 9).labels).unwrap())
            .collect();
        let best_k = 2 + s.iter().enumerate().max_by(|a, b| a.1.total_cmp(b.1)).unwrap().0;
        assert_eq!(best_k, 4, "{s:?}");
    }
}
