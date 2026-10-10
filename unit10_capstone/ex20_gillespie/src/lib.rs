//! # Exercise 20 — Stochastic chemical kinetics: a mini reaction-network simulator
//!
//! **Concepts:** a small DSL parser with good errors, data modelling with
//! index-based references (species ids), the Gillespie direct method, an
//! RK4 ODE solver for the deterministic limit, binary search, and
//! validating a stochastic simulator *statistically*.
//!
//! Inside a cell, many species exist at copy numbers of 1-100, so
//! deterministic ODEs are wrong. Noise is the physics. Gillespie's SSA
//! simulates every single reaction event exactly.
//!
//! **Network language** (one reaction per line):
//! ```text
//! # comments and blank lines are ignored; text after '#' too
//! 0 -> mRNA @ 10          # constitutive transcription (0 = nothing)
//! mRNA -> mRNA + P @ 5    # translation
//! mRNA -> 0 @ 1
//! 2 P -> P2 @ 0.001       # dimerisation: coefficient, space optional ("2P")
//! P2 -> 2 P @ 0.1
//! ```
//! * Species names: a letter followed by letters, digits or `_`.
//! * Terms are separated by `+`. `A + A` is the same as `2 A`.
//! * Species are numbered in order of first appearance in the file
//!   (left to right, top to bottom).
//! * The rate after `@` must parse as a finite, non-negative `f64`.
//!
//! **Propensity** (stochastic mass action): `rate * Π_s C(x_s, n_s)` where
//! `C` is the binomial coefficient. That's why `2 A` with 5 molecules has
//! 10 possible pairs, not 25.
//!
//! **ODE limit**: `d[x]/dt` from reaction rate `rate * Π_s x_s^n_s / n_s!`
//! (consistent with the propensity when counts are large).
//!
//! ```text
//! cargo test -p ex20_gillespie
//! ```

#![allow(unused_variables, dead_code, unused_imports, unused_mut)]

use std::collections::HashMap;

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
    /// Uniform in [0, 1).
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Reaction {
    /// `(species id, stoichiometric coefficient)`, merged per species, in
    /// order of first appearance within the side.
    pub reactants: Vec<(usize, u32)>,
    pub products: Vec<(usize, u32)>,
    pub rate: f64,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Network {
    pub species: Vec<String>,
    pub reactions: Vec<Reaction>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ParseErrorKind {
    MissingArrow,
    MissingRate,
    BadRate(String),
    BadTerm(String),
}

#[derive(Debug, Clone, PartialEq)]
pub struct ParseError {
    /// 1-based line number
    pub line: usize,
    pub kind: ParseErrorKind,
}

pub fn parse_network(src: &str) -> Result<Network, ParseError> {
    todo!()
}

/// Binomial coefficient as f64 (0 if k > n).
pub fn choose(n: u64, k: u32) -> f64 {
    todo!()
}

pub fn propensity(r: &Reaction, x: &[u64]) -> f64 {
    todo!()
}

/// Apply one firing of `r` to `x`. Panics (debug) on underflow; that would
/// mean you fired a reaction with zero propensity.
pub fn fire(r: &Reaction, x: &mut [u64]) {
    todo!()
}

#[derive(Debug, Clone, PartialEq)]
pub struct Trajectory {
    /// strictly increasing; `times[0] == 0.0`
    pub times: Vec<f64>,
    /// `states[i]` holds from `times[i]` until `times[i + 1]`
    pub states: Vec<Vec<u64>>,
}

impl Trajectory {
    /// State in effect at time `t` (the last event at or before `t`).
    /// Use binary search. Panics if `t < 0`.
    pub fn state_at(&self, t: f64) -> &[u64] {
        todo!()
    }
}

/// Gillespie direct method from t = 0 to `t_end`. Each step:
/// 1. `a0 = Σ propensities`; if `a0 == 0` stop (nothing can happen).
/// 2. `tau = -ln(1 - u1) / a0` with `u1 = rng.next_f64()`.
/// 3. if `t + tau > t_end` stop.
/// 4. choose reaction `j` with probability `a_j / a0` using `u2 = rng.next_f64()`
///    (the first `j` whose cumulative sum exceeds `u2 * a0`), fire it, record.
pub fn gillespie(net: &Network, init: &[u64], t_end: f64, rng: &mut Rng) -> Trajectory {
    todo!()
}

/// Deterministic mass-action ODE, classic RK4 with step `dt` (the last step
/// is shortened to land exactly on `t_end`). Returns the state at `t_end`.
pub fn ode_rk4(net: &Network, x0: &[f64], t_end: f64, dt: f64) -> Vec<f64> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    const GENE: &str = "
# toy gene expression
0 -> mRNA @ 10          # transcription
mRNA -> mRNA + P @ 5
mRNA -> 0 @ 1

2 P -> P2 @ 0.001
P2 -> 2P @ 0.1
P + P -> P2 @ 0
";

    #[test]
    fn parses_network() {
        let n = parse_network(GENE).unwrap();
        assert_eq!(n.species, vec!["mRNA", "P", "P2"]);
        assert_eq!(n.reactions.len(), 6);
        assert_eq!(n.reactions[0], Reaction { reactants: vec![], products: vec![(0, 1)], rate: 10.0 });
        assert_eq!(n.reactions[1].products, vec![(0, 1), (1, 1)]);
        assert_eq!(n.reactions[3].reactants, vec![(1, 2)]);
        assert_eq!(n.reactions[4].products, vec![(1, 2)]);
        assert_eq!(n.reactions[5].reactants, vec![(1, 2)]); // P + P merged
        assert_eq!(n.reactions[5].rate, 0.0);
        assert_eq!(parse_network("").unwrap(), Network::default());
    }

    #[test]
    fn parse_errors() {
        use ParseErrorKind::*;
        let e = |s: &str| parse_network(s).unwrap_err();
        assert_eq!(e("A -> B"), ParseError { line: 1, kind: MissingRate });
        assert_eq!(e("\nA => B @ 1"), ParseError { line: 2, kind: MissingArrow });
        assert_eq!(e("A -> B @ fast"), ParseError { line: 1, kind: BadRate("fast".into()) });
        assert_eq!(e("A -> B @ -1"), ParseError { line: 1, kind: BadRate("-1".into()) });
        assert_eq!(e("A -> B @ inf"), ParseError { line: 1, kind: BadRate("inf".into()) });
        assert_eq!(e("A + -> B @ 1"), ParseError { line: 1, kind: BadTerm("".into()) });
        assert_eq!(e("3x$ -> B @ 1"), ParseError { line: 1, kind: BadTerm("3x$".into()) });
        assert_eq!(e("1A -> 0A @ 1"), ParseError { line: 1, kind: BadTerm("0A".into()) });
    }

    #[test]
    fn propensities() {
        assert_eq!(choose(5, 2), 10.0);
        assert_eq!(choose(1, 2), 0.0);
        assert_eq!(choose(7, 0), 1.0);
        let n = parse_network("2 A -> B @ 1\nA + C -> 0 @ 0.5\n0 -> A @ 3\n3A -> 0 @ 2").unwrap();
        let x = [5, 0, 4]; // A, B, C
        assert_eq!(propensity(&n.reactions[0], &x), 10.0);
        assert_eq!(propensity(&n.reactions[1], &x), 10.0);
        assert_eq!(propensity(&n.reactions[2], &x), 3.0);
        assert_eq!(propensity(&n.reactions[3], &x), 20.0);
        assert_eq!(propensity(&n.reactions[0], &[1, 0, 0]), 0.0);
        let mut y = x;
        fire(&n.reactions[0], &mut y);
        assert_eq!(y, [3, 1, 4]);
    }

    #[test]
    fn trajectory_structure_and_absorption() {
        let n = parse_network("A -> 0 @ 1").unwrap();
        let tr = gillespie(&n, &[5], 1e9, &mut Rng::new(1));
        assert_eq!(tr.times.len(), 6);
        assert_eq!(tr.states.last().unwrap(), &vec![0]);
        assert!(tr.times.windows(2).all(|w| w[0] < w[1]));
        assert_eq!(tr.state_at(0.0), &[5]);
        assert_eq!(tr.state_at(tr.times[2]), &[3]);
        assert_eq!(tr.state_at((tr.times[2] + tr.times[3]) / 2.0), &[3]);
        assert_eq!(tr.state_at(1e10), &[0]);
        let short = gillespie(&n, &[1000], 0.01, &mut Rng::new(2));
        assert!(*short.times.last().unwrap() <= 0.01);
    }

    #[test]
    fn deterministic_given_seed() {
        let n = parse_network(GENE).unwrap();
        let a = gillespie(&n, &[0, 0, 0], 5.0, &mut Rng::new(42));
        let b = gillespie(&n, &[0, 0, 0], 5.0, &mut Rng::new(42));
        assert_eq!(a, b);
    }

    #[test]
    fn birth_death_is_poisson() {
        // stationary distribution of 0 -> X @ k, X -> 0 @ g is Poisson(k/g)
        let n = parse_network("0 -> X @ 10\nX -> 0 @ 1").unwrap();
        let tr = gillespie(&n, &[0], 3000.0, &mut Rng::new(7));
        let samples: Vec<f64> = (50..3000).map(|t| tr.state_at(t as f64)[0] as f64).collect();
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        let var = samples.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (samples.len() - 1) as f64;
        assert!((mean - 10.0).abs() < 0.5, "mean {mean}");
        assert!((var - 10.0).abs() < 1.5, "var {var}");
    }

    #[test]
    fn conservation_law() {
        // A + 2*B is conserved by dimerisation/dissociation
        let n = parse_network("A + A -> B @ 0.01\nB -> 2 A @ 0.5").unwrap();
        let tr = gillespie(&n, &[100, 0], 50.0, &mut Rng::new(3));
        assert!(tr.times.len() > 100);
        assert!(tr.states.iter().all(|s| s[0] + 2 * s[1] == 100));
    }

    #[test]
    fn decay_mean_matches_exponential() {
        let n = parse_network("A -> 0 @ 0.5").unwrap();
        let mut rng = Rng::new(11);
        let runs = 200;
        let mean: f64 = (0..runs)
            .map(|_| gillespie(&n, &[1000], 2.0, &mut rng).states.last().unwrap()[0] as f64)
            .sum::<f64>()
            / runs as f64;
        let expect = 1000.0 * (-1.0f64).exp();
        assert!((mean - expect).abs() < 5.0, "{mean} vs {expect}");
    }

    #[test]
    fn ode_limit() {
        let n = parse_network("0 -> X @ 10\nX -> 0 @ 1").unwrap();
        let x = ode_rk4(&n, &[0.0], 3.0, 0.01);
        assert!((x[0] - 10.0 * (1.0 - (-3.0f64).exp())).abs() < 1e-8);
        // dt that doesn't divide t_end
        let x = ode_rk4(&n, &[0.0], 1.0, 0.3);
        assert!((x[0] - 10.0 * (1.0 - (-1.0f64).exp())).abs() < 1e-3);
        // 2A -> B with rate k: d[A]/dt = -2 * k * A^2 / 2 = -k A^2 -> A(t) = A0 / (1 + k A0 t)
        let n = parse_network("2 A -> B @ 0.01").unwrap();
        let x = ode_rk4(&n, &[100.0, 0.0], 2.0, 0.001);
        assert!((x[0] - 100.0 / 3.0).abs() < 1e-6, "{x:?}");
        assert!((x[0] + 2.0 * x[1] - 100.0).abs() < 1e-9);
    }

    #[test]
    fn ssa_mean_approaches_ode_for_large_numbers() {
        let n = parse_network("2 A -> B @ 0.0001").unwrap();
        let mut rng = Rng::new(5);
        let runs = 30;
        let mean_a = (0..runs)
            .map(|_| gillespie(&n, &[2000, 0], 5.0, &mut rng).states.last().unwrap()[0] as f64)
            .sum::<f64>()
            / runs as f64;
        let ode = ode_rk4(&n, &[2000.0, 0.0], 5.0, 0.001)[0];
        assert!((mean_a - ode).abs() / ode < 0.01, "{mean_a} vs {ode}");
    }
}
