//! # Exercise 14 — Molecular dynamics of a Lennard-Jones fluid
//!
//! **Concepts:** using your own crate as a dependency (this crate uses
//! *your* `Vec3` from exercise 13, so finish that first), `&mut self`
//! methods that update several fields, avoiding double borrows when you
//! update `force[i]` and `force[j]` in the same loop, numerical
//! verification through physical invariants.
//!
//! Simulate argon-like atoms in **reduced units** (ε = σ = m = k_B = 1):
//!
//! * Potential: `U(r) = 4 (r⁻¹² - r⁻⁶) - U_c` for `r < r_c`, else 0, where
//!   `U_c = 4 (r_c⁻¹² - r_c⁻⁶)` shifts the potential to 0 at the cutoff.
//! * Force on `i` from `j`: `F = 24 (2 r⁻¹⁴ - r⁻⁸) (r_i - r_j)` (with the
//!   minimum-image displacement).
//! * Periodic cubic box of side `L`: positions live in `[0, L)`, and the
//!   displacement uses the **minimum image convention**.
//! * Integrator: **velocity Verlet**:
//!   `v += F dt/2; x += v dt; wrap x; recompute F; v += F dt/2`.
//! * Temperature: `T = 2 KE / (3 N)` (keep it simple: don't subtract the
//!   3 constrained COM degrees of freedom).
//!
//! The tests are physics: Newton's third law, momentum conservation and
//! **energy conservation** over 1000 steps. If your total energy drifts,
//! you have a bug, most likely in the force sign, the minimum image or the
//! order of the Verlet half-steps.
//!
//! Stretch goal (tested): the radial distribution function g(r).
//!
//! ```text
//! cargo test -p ex14_lennard_jones_md
//! ```

#![allow(unused_variables, dead_code, unused_imports, unused_mut)]

pub use ex13_vec3_geometry::Vec3;

/// SplitMix64 + Box–Muller. Provided.
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
    /// Standard normal.
    pub fn gaussian(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-300);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

/// Shortest periodic image of displacement `d` in a cubic box of side `l`.
pub fn minimum_image(d: Vec3, l: f64) -> Vec3 {
    todo!()
}

/// Wrap a position into `[0, l)` on every axis.
pub fn wrap(p: Vec3, l: f64) -> Vec3 {
    todo!()
}

/// `n³` points on a simple cubic lattice filling a box of side `l`,
/// spacing `l / n`, first point at `(a/2, a/2, a/2)`. Order: x varies
/// slowest, then y, then z fastest.
pub fn cubic_lattice(n: usize, l: f64) -> Vec<Vec3> {
    todo!()
}

/// Shifted LJ pair energy and the scalar `f` such that the force on `i` is
/// `f * (r_i - r_j)`. Takes the *squared* distance (avoid square roots in the
/// hot loop). Returns `(0, 0)` at or beyond the cutoff.
pub fn lj_pair(r2: f64, cutoff: f64) -> (f64, f64) {
    todo!()
}

#[derive(Debug, Clone)]
pub struct System {
    pub pos: Vec<Vec3>,
    pub vel: Vec<Vec3>,
    pub force: Vec<Vec3>,
    pub box_len: f64,
    pub cutoff: f64,
    /// potential energy matching the current `force`
    pub potential: f64,
}

impl System {
    /// Zero velocities. Positions are wrapped into the box. Forces and
    /// `potential` must be valid on return. Panics if `cutoff > box_len / 2`.
    pub fn new(pos: Vec<Vec3>, box_len: f64, cutoff: f64) -> Self {
        todo!()
    }

    pub fn n(&self) -> usize {
        self.pos.len()
    }

    /// Recompute `self.force` and `self.potential` (all pairs, each pair
    /// visited once, Newton's third law). Returns the potential energy.
    pub fn compute_forces(&mut self) -> f64 {
        todo!()
    }

    pub fn kinetic_energy(&self) -> f64 {
        todo!()
    }

    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy() + self.potential
    }

    pub fn temperature(&self) -> f64 {
        todo!()
    }

    pub fn momentum(&self) -> Vec3 {
        todo!()
    }

    /// Gaussian velocities, then remove the centre-of-mass velocity, then
    /// rescale so `temperature()` is exactly `t`.
    pub fn init_velocities(&mut self, t: f64, rng: &mut Rng) {
        todo!()
    }

    /// One velocity-Verlet step.
    pub fn step(&mut self, dt: f64) {
        todo!()
    }

    /// Radial distribution function with `nbins` bins on `[0, r_max)`
    /// (`r_max <= L/2`). Returns `(bin_centre, g)` for each bin, normalised
    /// so that an ideal gas gives g ≈ 1:
    /// `g = 2 * pair_count / (N * ρ * shell_volume)` with `ρ = N / L³` and
    /// `shell_volume = 4/3 π (r_hi³ - r_lo³)`.
    pub fn rdf(&self, nbins: usize, r_max: f64) -> Vec<(f64, f64)> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    fn liquid(seed: u64) -> System {
        let n_side = 6; // 216 atoms
        let rho: f64 = 0.8;
        let l = (216.0 / rho).cbrt();
        let mut s = System::new(cubic_lattice(n_side, l), l, 2.5);
        s.init_velocities(1.2, &mut Rng::new(seed));
        s
    }

    #[test]
    fn periodic_helpers() {
        let l = 10.0;
        assert_eq!(minimum_image(Vec3::new(9.0, -9.0, 4.0), l), Vec3::new(-1.0, 1.0, 4.0));
        assert_eq!(minimum_image(Vec3::new(21.0, 0.0, 0.0), l), Vec3::new(1.0, 0.0, 0.0));
        assert_eq!(wrap(Vec3::new(-1.0, 10.0, 25.5), l), Vec3::new(9.0, 0.0, 5.5));
        let w = wrap(Vec3::new(-1e-17, 0.0, 0.0), l);
        assert!(w.x >= 0.0 && w.x < l, "{w:?}");
    }

    #[test]
    fn lattice() {
        let p = cubic_lattice(2, 4.0);
        assert_eq!(p.len(), 8);
        assert_eq!(p[0], Vec3::new(1.0, 1.0, 1.0));
        assert_eq!(p[1], Vec3::new(1.0, 1.0, 3.0));
        assert_eq!(p[7], Vec3::new(3.0, 3.0, 3.0));
    }

    #[test]
    fn pair_potential() {
        let rc: f64 = 2.5;
        let shift = 4.0 * (rc.powi(-12) - rc.powi(-6));
        let (u, f) = lj_pair(1.0, rc);
        assert!(close(u, -shift, 1e-12));
        assert!(close(f, 24.0, 1e-12));
        let rmin = 2f64.powf(1.0 / 6.0);
        let (u, f) = lj_pair(rmin * rmin, rc);
        assert!(close(u, -1.0 - shift, 1e-12));
        assert!(close(f, 0.0, 1e-12));
        assert_eq!(lj_pair(6.25, rc), (0.0, 0.0));
        // f * r should equal -dU/dr: finite-difference check
        let r = 1.3;
        let h = 1e-6;
        let du = (lj_pair((r + h) * (r + h), rc).0 - lj_pair((r - h) * (r - h), rc).0) / (2.0 * h);
        assert!(close(lj_pair(r * r, rc).1 * r, -du, 1e-6));
    }

    #[test]
    fn two_atoms_across_boundary() {
        let l = 10.0;
        let s = System::new(vec![Vec3::new(0.25, 5.0, 5.0), Vec3::new(9.75, 5.0, 5.0)], l, 2.5);
        // they are 0.5 apart through the boundary -> strong repulsion, pushes atom 0 in +x
        assert!(s.force[0].x > 1000.0);
        assert!(close(s.force[0].x, -s.force[1].x, 1e-9));
        assert!(close(s.force[0].y, 0.0, 1e-12));
    }

    #[test]
    fn newtons_third_law_and_momentum() {
        let mut s = liquid(1);
        let fsum: Vec3 = s.force.iter().sum();
        assert!(fsum.norm() < 1e-9, "{fsum:?}");
        assert!(s.momentum().norm() < 1e-9);
        assert!(close(s.temperature(), 1.2, 1e-12));
        for _ in 0..200 {
            s.step(0.004);
        }
        assert!(s.momentum().norm() < 1e-8, "{:?}", s.momentum());
        for p in &s.pos {
            for d in 0..3 {
                assert!(p[d] >= 0.0 && p[d] < s.box_len);
            }
        }
    }

    #[test]
    fn energy_conservation() {
        let mut s = liquid(7);
        // let the lattice melt first
        for _ in 0..200 {
            s.step(0.002);
        }
        let e0 = s.total_energy();
        let mut emax: f64 = 0.0;
        for _ in 0..1000 {
            s.step(0.002);
            emax = emax.max((s.total_energy() - e0).abs());
        }
        let per_atom = emax / s.n() as f64;
        assert!(per_atom < 2e-3, "energy drift {per_atom} per atom");
        // the lattice melted and the system equilibrated to a sensible temperature
        assert!(s.temperature() > 0.3 && s.temperature() < 1.5, "T = {}", s.temperature());
    }

    #[test]
    fn smaller_timestep_conserves_better() {
        let drift = |dt: f64| {
            let mut s = liquid(3);
            for _ in 0..(0.4 / dt) as usize {
                s.step(dt);
            }
            let e0 = s.total_energy();
            let mut m: f64 = 0.0;
            for _ in 0..(1.0 / dt) as usize {
                s.step(dt);
                m = m.max((s.total_energy() - e0).abs());
            }
            m
        };
        let (d1, d2) = (drift(0.005), drift(0.0025));
        // velocity Verlet error scales ~ dt^2
        assert!(d2 < d1 * 0.5, "{d1} vs {d2}");
    }

    #[test]
    fn rdf_of_liquid() {
        let mut s = liquid(11);
        for _ in 0..400 {
            s.step(0.004);
        }
        let g = s.rdf(30, s.box_len / 2.0);
        // excluded volume: nothing very close
        assert!(g.iter().filter(|(r, _)| *r < 0.8).all(|(_, g)| *g == 0.0));
        // first-neighbour peak around r ≈ 1.1, clearly above 1
        let (rpeak, gpeak) = g.iter().cloned().max_by(|a, b| a.1.total_cmp(&b.1)).unwrap();
        assert!(rpeak > 0.95 && rpeak < 1.3 && gpeak > 1.8, "peak {gpeak} at {rpeak}");
        // tends to 1 at long range
        let tail: Vec<f64> = g.iter().filter(|(r, _)| *r > 2.5).map(|x| x.1).collect();
        let avg = tail.iter().sum::<f64>() / tail.len() as f64;
        assert!((avg - 1.0).abs() < 0.15, "tail average {avg}");
    }
}
