//! # Exercise 18 — A deterministic parallel Ising model
//!
//! **Concepts:** splitting a mutable buffer between threads
//! (`chunks_mut`, `split_at_mut`), sharing read-only data with scoped
//! threads, designing an algorithm so that parallel results are
//! **bit-for-bit identical** to serial ones, per-task RNG streams.
//!
//! The 2-D Ising model is the "hello world" of statistical physics and a
//! toy model for cooperative binding, membrane lipid domains and protein
//! folding cooperativity. Spins `s = ±1` on an `L x L` periodic lattice:
//!
//! `E = -J Σ_<ij> s_i s_j - h Σ_i s_i` (each nearest-neighbour bond once).
//!
//! **Checkerboard Metropolis.** Colour site `(r, c)` black if `(r + c)` is
//! even. A sweep updates all black sites, then all white sites. Black sites
//! only have white neighbours, so all black sites can be updated
//! *simultaneously*. That is what makes it parallelisable. (Requires even L.)
//!
//! **Determinism contract.** For half-sweep `colour` (0 = black, 1 = white)
//! of sweep number `sweep`, row `r` uses its own RNG
//! `Rng::new(stream_seed(seed, sweep, colour, r))` and visits that row's
//! sites of that colour in increasing column order. For each site draw
//! `u = rng.next_f64()` **always** (even if the move is downhill), and flip
//! iff `dE <= 0 || u < exp(-beta * dE)`. Any thread split of rows then gives
//! *identical* lattices. The tests check that across 1-7 threads.
//!
//! The challenge: in `sweep_parallel`, a thread updating row `r` must read
//! rows `r-1` and `r+1`, which may belong to other threads. You can't hand
//! out `&mut` chunks *and* `&` to neighbours of the same `Vec`. Find a safe
//! way (no `unsafe`!). There are at least two good designs.
//!
//! ```text
//! cargo test -p ex18_ising_parallel
//! ```

#![allow(unused_variables, dead_code, unused_imports, unused_mut)]

use std::thread;

/// SplitMix64. Provided.
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
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Seed for one row's RNG in one half-sweep. Provided.
pub fn stream_seed(seed: u64, sweep: u64, colour: u64, row: u64) -> u64 {
    let mut r = Rng::new(seed ^ sweep.wrapping_mul(0x9E3779B97F4A7C15) ^ (colour << 63) ^ row.wrapping_mul(0xD1B54A32D192ED03));
    r.next_u64()
}

#[derive(Debug, Clone, PartialEq)]
pub struct Ising {
    pub l: usize,
    /// Row-major, each entry +1 or -1.
    pub spins: Vec<i8>,
    pub j: f64,
    pub h: f64,
}

impl Ising {
    /// All spins up. Panics unless `l` is even and >= 2 (the panic message
    /// must contain the word "even").
    pub fn new_ordered(l: usize, j: f64, h: f64) -> Self {
        todo!()
    }

    /// Random spins: site i is +1 iff `rng.next_f64() < 0.5`, row-major order.
    pub fn new_random(l: usize, j: f64, h: f64, seed: u64) -> Self {
        todo!()
    }

    pub fn get(&self, r: usize, c: usize) -> i8 {
        self.spins[r * self.l + c]
    }

    /// Total energy (each bond counted once, periodic boundaries).
    pub fn energy(&self) -> f64 {
        todo!()
    }

    /// Mean spin, in [-1, 1].
    pub fn magnetization(&self) -> f64 {
        todo!()
    }

    /// Energy change if spin (r, c) were flipped.
    pub fn delta_e(&self, r: usize, c: usize) -> f64 {
        todo!()
    }

    /// One full sweep (black then white), serially, following the contract.
    pub fn sweep_serial(&mut self, beta: f64, seed: u64, sweep: u64) {
        todo!()
    }

    /// Same result as `sweep_serial`, rows split across `n_threads` scoped
    /// threads. Panics if `n_threads == 0`.
    pub fn sweep_parallel(&mut self, beta: f64, seed: u64, sweep: u64, n_threads: usize) {
        todo!()
    }

}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn energies() {
        let s = Ising::new_ordered(8, 1.0, 0.0);
        assert_eq!(s.energy(), -2.0 * 64.0);
        assert_eq!(s.magnetization(), 1.0);
        let s = Ising::new_ordered(4, 1.0, 0.5);
        assert_eq!(s.energy(), -32.0 - 8.0);
        // checkerboard (antiferromagnetic) state
        let mut s = Ising::new_ordered(4, 1.0, 0.0);
        for r in 0..4 {
            for c in 0..4 {
                if (r + c) % 2 == 1 {
                    s.spins[r * 4 + c] = -1;
                }
            }
        }
        assert_eq!(s.energy(), 32.0);
        assert_eq!(s.magnetization(), 0.0);
    }

    #[test]
    fn delta_e_matches_energy_difference() {
        let s = Ising::new_random(6, 0.7, -0.3, 5);
        for r in 0..6 {
            for c in 0..6 {
                let mut t = s.clone();
                t.spins[r * 6 + c] *= -1;
                let diff = t.energy() - s.energy();
                assert!((s.delta_e(r, c) - diff).abs() < 1e-9, "({r},{c})");
            }
        }
    }

    #[test]
    #[should_panic(expected = "even")]
    fn odd_l_rejected() {
        Ising::new_ordered(5, 1.0, 0.0);
    }

    #[test]
    fn random_init_is_balanced() {
        let s = Ising::new_random(64, 1.0, 0.0, 1);
        assert!(s.magnetization().abs() < 0.05);
        assert!(s.spins.iter().all(|&x| x == 1 || x == -1));
    }

    #[test]
    fn parallel_is_bitwise_identical() {
        for &(l, beta) in &[(16usize, 0.3f64), (20, 0.44), (10, 1.0)] {
            let mut serial = Ising::new_random(l, 1.0, 0.1, 42);
            let mut snapshots = vec![];
            for sweep in 0..25 {
                serial.sweep_serial(beta, 7, sweep);
                snapshots.push(serial.spins.clone());
            }
            for threads in 1..=7 {
                let mut par = Ising::new_random(l, 1.0, 0.1, 42);
                for sweep in 0..25 {
                    par.sweep_parallel(beta, 7, sweep, threads);
                    assert_eq!(par.spins, snapshots[sweep as usize], "L={l} threads={threads} sweep={sweep}");
                }
            }
        }
    }

    #[test]
    fn more_threads_than_rows() {
        let mut a = Ising::new_random(4, 1.0, 0.0, 3);
        let mut b = a.clone();
        a.sweep_serial(0.5, 1, 0);
        b.sweep_parallel(0.5, 1, 0, 64);
        assert_eq!(a, b);
    }

    #[test]
    fn physics_low_and_high_temperature() {
        // ordered phase: T = 1/beta = 1.5 < T_c ≈ 2.269
        let mut cold = Ising::new_ordered(32, 1.0, 0.0);
        for sweep in 0..300 {
            cold.sweep_parallel(1.0 / 1.5, 11, sweep, 4);
        }
        assert!(cold.magnetization().abs() > 0.9, "m = {}", cold.magnetization());
        // disordered phase: T = 5
        let mut hot = Ising::new_ordered(32, 1.0, 0.0);
        let mut m_avg = 0.0;
        for sweep in 0..600 {
            hot.sweep_parallel(1.0 / 5.0, 11, sweep, 4);
            if sweep >= 100 {
                m_avg += hot.magnetization().abs() / 500.0;
            }
        }
        assert!(m_avg < 0.15, "<|m|> = {m_avg}");
    }

    #[test]
    fn energy_matches_onsager_exact_solution() {
        // Onsager's exact infinite-lattice energy per spin at T = 3 (J = 1) is -0.8173.
        let mut s = Ising::new_random(32, 1.0, 0.0, 2);
        let beta = 1.0 / 3.0;
        for sweep in 0..200 {
            s.sweep_parallel(beta, 99, sweep, 3);
        }
        let mut e = 0.0;
        for sweep in 200..1200 {
            s.sweep_parallel(beta, 99, sweep, 3);
            e += s.energy() / (32.0 * 32.0) / 1000.0;
        }
        assert!((e - (-0.8173)).abs() < 0.03, "e = {e}");
    }
}
