//! # Exercise 13 — `Vec3`, operator overloading and molecular geometry
//!
//! **Concepts:** `std::ops` traits (`Add`, `Sub`, `Mul`, `Neg`,
//! `AddAssign`, `Index`, `IndexMut`), implementing a trait for a foreign
//! type on the left-hand side (`f64 * Vec3`), `Sum`, `Copy` types,
//! `#[must_use]`.
//!
//! You'll build the vector type used by the molecular-dynamics exercise
//! next, then use it for the three internal coordinates every structural
//! biologist uses: bond length, bond angle and **dihedral (torsion)
//! angle**, the φ/ψ angles of a Ramachandran plot.
//!
//! The dihedral sign convention matters: use the IUPAC convention (the same
//! as MDTraj / Biopython's `calc_dihedral`), result in `(-π, π]`. Use
//! `atan2` (not `acos`) so you get the sign and good precision.
//!
//! **Note:** the tests use `+`, `-`, `*`, `[]`, `.sum()` on `Vec3`, so this
//! crate's *tests* won't compile until you've written the operator impls.
//! Treat the compiler errors as your to-do list.
//!
//! ```text
//! cargo test -p ex13_vec3_geometry
//! ```

#![allow(unused_variables, dead_code, unused_imports)]

use std::iter::Sum;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Neg, Sub, SubAssign};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const ZERO: Vec3 = Vec3 { x: 0.0, y: 0.0, z: 0.0 };

    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    #[must_use]
    pub fn dot(self, o: Vec3) -> f64 {
        todo!()
    }

    #[must_use]
    pub fn cross(self, o: Vec3) -> Vec3 {
        todo!()
    }

    #[must_use]
    pub fn norm_sq(self) -> f64 {
        todo!()
    }

    #[must_use]
    pub fn norm(self) -> f64 {
        todo!()
    }

    /// Unit vector; `None` for the zero vector.
    #[must_use]
    pub fn normalized(self) -> Option<Vec3> {
        todo!()
    }

    /// Apply `f` to each component.
    #[must_use]
    pub fn map(self, f: impl Fn(f64) -> f64) -> Vec3 {
        todo!()
    }
}

// TODO: implement
//   Vec3 + Vec3, Vec3 - Vec3, Vec3 * f64, f64 * Vec3, -Vec3,
//   +=, -=, v[0..3] (read and write; panic with a message containing
//   "out of range" for any other index),
//   Sum<Vec3> and Sum<&Vec3> (so `.iter().sum()` works)

/// |b - a|
pub fn distance(a: Vec3, b: Vec3) -> f64 {
    todo!()
}

/// Angle a-b-c at vertex `b`, radians in [0, π].
pub fn bond_angle(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    todo!()
}

/// Dihedral angle p0-p1-p2-p3 in radians, range (-π, π].
pub fn dihedral(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3) -> f64 {
    todo!()
}

/// Geometric centre (unweighted). `None` if empty.
pub fn centroid(points: &[Vec3]) -> Option<Vec3> {
    todo!()
}

/// Root-mean-square deviation between two equally long point sets
/// **without** superposition. `None` if lengths differ or empty.
pub fn rmsd(a: &[Vec3], b: &[Vec3]) -> Option<f64> {
    todo!()
}

/// Place atom D given A, B, C and internal coordinates: |CD| = `bond`,
/// angle B-C-D = `angle`, dihedral A-B-C-D = `torsion` (radians).
/// This is the NeRF algorithm used to build proteins from φ/ψ angles.
/// Must be the exact inverse of `distance`/`bond_angle`/`dihedral`.
pub fn place_atom(a: Vec3, b: Vec3, c: Vec3, bond: f64, angle: f64, torsion: f64) -> Vec3 {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn operators() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(-1.0, 0.5, 2.0);
        assert_eq!(a + b, Vec3::new(0.0, 2.5, 5.0));
        assert_eq!(a - b, Vec3::new(2.0, 1.5, 1.0));
        assert_eq!(a * 2.0, Vec3::new(2.0, 4.0, 6.0));
        assert_eq!(2.0 * a, a * 2.0);
        assert_eq!(-a, Vec3::new(-1.0, -2.0, -3.0));
        let mut c = a;
        c += b;
        c -= a;
        assert_eq!(c, b);
        c[2] = 9.0;
        assert_eq!((c[0], c[1], c[2]), (-1.0, 0.5, 9.0));
        let v = vec![a, b, a];
        assert_eq!(v.iter().sum::<Vec3>(), Vec3::new(1.0, 4.5, 8.0));
        assert_eq!(v.into_iter().sum::<Vec3>(), Vec3::new(1.0, 4.5, 8.0));
    }

    #[test]
    #[should_panic(expected = "out of range")]
    fn index_out_of_range() {
        let _ = Vec3::ZERO[3];
    }

    #[test]
    fn products_and_norms() {
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 1.0, 0.0);
        assert_eq!(x.cross(y), Vec3::new(0.0, 0.0, 1.0));
        assert_eq!(y.cross(x), Vec3::new(0.0, 0.0, -1.0));
        assert_eq!(Vec3::new(3.0, 4.0, 0.0).norm(), 5.0);
        assert_eq!(Vec3::ZERO.normalized(), None);
        assert!(close(Vec3::new(2.0, -3.0, 6.0).normalized().unwrap().norm(), 1.0));
        assert_eq!(x.map(|v| v + 1.0), Vec3::new(2.0, 1.0, 1.0));
    }

    #[test]
    fn angles() {
        let o = Vec3::ZERO;
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 2.0, 0.0);
        assert!(close(bond_angle(x, o, y), PI / 2.0));
        assert!(close(bond_angle(x, o, x * 3.0), 0.0));
        assert!(close(bond_angle(x, o, -x), PI));
        // water: 104.5 degrees
        let h1 = Vec3::new(0.9572, 0.0, 0.0);
        let t = 104.5f64.to_radians();
        let h2 = Vec3::new(0.9572 * t.cos(), 0.9572 * t.sin(), 0.0);
        assert!(close(bond_angle(h1, o, h2).to_degrees(), 104.5));
    }

    #[test]
    fn dihedral_sign_convention() {
        let p0 = Vec3::new(1.0, 0.0, 0.0);
        let p1 = Vec3::new(0.0, 0.0, 0.0);
        let p2 = Vec3::new(0.0, 0.0, 1.0);
        // cis
        assert!(close(dihedral(p0, p1, p2, Vec3::new(1.0, 0.0, 1.0)), 0.0));
        // trans -> exactly +π, never -π
        assert_eq!(dihedral(p0, p1, p2, Vec3::new(-1.0, 0.0, 1.0)), PI);
        // ±90°, with the standard (IUPAC / Biopython / MDTraj) sign
        assert!(close(dihedral(p0, p1, p2, Vec3::new(0.0, 1.0, 1.0)), PI / 2.0));
        assert!(close(dihedral(p0, p1, p2, Vec3::new(0.0, -1.0, 1.0)), -PI / 2.0));
        // an obtuse negative torsion
        let d = dihedral(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(-1.0, 1.0, 1.0),
        );
        assert!(close(d, -3.0 * PI / 4.0), "{d}");
    }

    #[test]
    fn centroid_and_rmsd() {
        let pts = [Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0), Vec3::new(1.0, 3.0, 0.0)];
        assert_eq!(centroid(&pts), Some(Vec3::new(1.0, 1.0, 0.0)));
        assert_eq!(centroid(&[]), None);
        let shifted: Vec<Vec3> = pts.iter().map(|&p| p + Vec3::new(0.0, 0.0, 2.0)).collect();
        assert!(close(rmsd(&pts, &shifted).unwrap(), 2.0));
        assert_eq!(rmsd(&pts, &shifted[..2]), None);
    }

    #[test]
    fn nerf_roundtrip() {
        let a = Vec3::new(1.2, -0.3, 0.5);
        let b = Vec3::new(0.1, 0.4, -0.2);
        let c = Vec3::new(-0.5, 1.6, 0.3);
        let cases: [(f64, f64, f64); 4] = [(1.53, 111.0, -60.0), (1.33, 120.0, 180.0), (1.0, 90.0, 37.0), (2.0, 150.0, -179.0)];
        for (bond, ang_deg, tor_deg) in cases {
            let d = place_atom(a, b, c, bond, ang_deg.to_radians(), tor_deg.to_radians());
            assert!(close(distance(c, d), bond));
            assert!(close(bond_angle(b, c, d).to_degrees(), ang_deg));
            let t = dihedral(a, b, c, d).to_degrees();
            let diff = ((t - tor_deg + 540.0) % 360.0) - 180.0;
            assert!(diff.abs() < 1e-7, "{t} vs {tor_deg}");
        }
    }
}
