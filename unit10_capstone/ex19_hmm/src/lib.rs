//! # Exercise 19 — Hidden Markov Models: Viterbi, Forward–Backward, Baum–Welch
//!
//! **Concepts:** everything so far. Validated constructors returning
//! `Result`, numerically stable algorithms (log space and scaling),
//! nested `Vec` vs flat storage, careful indexing, an EM algorithm whose
//! correctness is checked by an invariant (likelihood never decreases).
//!
//! HMMs power gene finders, profile HMMs (Pfam/HMMER), CpG-island
//! detection, ion-channel state inference and much more. Implement:
//!
//! * `Hmm::new`: validate shapes, non-negative entries, rows summing to 1
//!   (±1e-9).
//! * `viterbi`: most probable state path and its **natural-log**
//!   probability. Work in log space. Ties: prefer the lower state index.
//! * `forward`: *scaled* forward variables `α̂_t` (each row sums to 1) and
//!   scale factors `c_t`, so `ln P(obs) = Σ ln c_t`. Plain probabilities
//!   underflow to 0 after a few hundred symbols, and the tests use
//!   sequences of 5000.
//! * `backward` (scaled with the same `c_t`) and `posterior` (γ).
//! * `baum_welch_step`: one EM re-estimation over **multiple** sequences.
//!   If a state has zero expected occupancy, keep its old row(s).
//! * `train`: iterate until the log-likelihood gain is `< tol`.
//! * `cpg_islands`: use the provided two-state model and Viterbi to find
//!   islands in a DNA sequence.
//!
//! The tests check you against brute-force enumeration over every possible
//! state path on short sequences.
//!
//! ```text
//! cargo test -p ex19_hmm
//! ```

#![allow(unused_variables, dead_code, unused_imports, unused_mut)]

#[derive(Debug, Clone, PartialEq)]
pub struct Hmm {
    /// `start[i]` = P(state_0 = i)
    pub start: Vec<f64>,
    /// `trans[i][j]` = P(state_{t+1} = j | state_t = i)
    pub trans: Vec<Vec<f64>>,
    /// `emit[i][o]` = P(obs_t = o | state_t = i)
    pub emit: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum HmmError {
    /// No states, no symbols, or inconsistent dimensions.
    Shape,
    Negative,
    /// `which` is "start", "trans" or "emit"; `row` is 0 for start.
    NotStochastic { which: &'static str, row: usize },
}

impl Hmm {
    pub fn new(start: Vec<f64>, trans: Vec<Vec<f64>>, emit: Vec<Vec<f64>>) -> Result<Self, HmmError> {
        todo!()
    }

    pub fn n_states(&self) -> usize {
        self.start.len()
    }

    pub fn n_symbols(&self) -> usize {
        self.emit[0].len()
    }

    /// `(path, ln P(path, obs))`. Empty `obs` -> `(vec![], 0.0)`.
    /// If every path has probability 0, the log-probability is `-inf` and
    /// the path may be anything of the right length.
    pub fn viterbi(&self, obs: &[usize]) -> (Vec<usize>, f64) {
        todo!()
    }

    /// Scaled forward pass: `(alpha_hat, c)` with `alpha_hat[t]` summing to 1
    /// and `c[t]` the normaliser at step t. Panics if the sequence has
    /// probability 0.
    pub fn forward(&self, obs: &[usize]) -> (Vec<Vec<f64>>, Vec<f64>) {
        todo!()
    }

    /// Scaled backward pass using the forward scale factors `c`:
    /// `beta_hat[T-1] = 1`, `beta_hat[t] = (Σ_j a_ij b_j(o_{t+1}) beta_hat[t+1][j]) / c[t+1]`.
    pub fn backward(&self, obs: &[usize], c: &[f64]) -> Vec<Vec<f64>> {
        todo!()
    }

    /// Natural log of P(obs). Empty sequence -> 0.0.
    pub fn log_likelihood(&self, obs: &[usize]) -> f64 {
        todo!()
    }

    /// `gamma[t][i]` = P(state_t = i | obs).
    pub fn posterior(&self, obs: &[usize]) -> Vec<Vec<f64>> {
        todo!()
    }

    /// One Baum–Welch iteration over all sequences (empty ones are ignored).
    pub fn baum_welch_step(&self, seqs: &[Vec<usize>]) -> Hmm {
        todo!()
    }

    /// Run Baum–Welch until the total log-likelihood improves by less than
    /// `tol` or `max_iter` steps were taken. Returns the final model and the
    /// log-likelihood history (first entry = initial model).
    pub fn train(&self, seqs: &[Vec<usize>], max_iter: usize, tol: f64) -> (Hmm, Vec<f64>) {
        todo!()
    }
}

/// Two-state CpG model. State 0 = background, 1 = island. Symbols A,C,G,T = 0..4.
pub fn cpg_model() -> Hmm {
    Hmm::new(
        vec![0.99, 0.01],
        vec![vec![0.999, 0.001], vec![0.004, 0.996]],
        vec![vec![0.3, 0.2, 0.2, 0.3], vec![0.15, 0.35, 0.35, 0.15]],
    )
    .unwrap()
}

/// Half-open intervals of the input where the Viterbi path of `cpg_model()`
/// is in state 1 for at least `min_len` consecutive positions. Bases other
/// than ACGT (case-insensitive) are treated as... you decide, but document
/// it, and make sure an `N` doesn't panic.
pub fn cpg_islands(seq: &[u8], min_len: usize) -> Vec<(usize, usize)> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn casino() -> Hmm {
        // Durbin et al.: fair die (0) vs loaded die (1); symbols 0..6 = faces 1..6
        Hmm::new(
            vec![0.5, 0.5],
            vec![vec![0.95, 0.05], vec![0.1, 0.9]],
            vec![vec![1.0 / 6.0; 6], vec![0.1, 0.1, 0.1, 0.1, 0.1, 0.5]],
        )
        .unwrap()
    }

    fn three_state() -> Hmm {
        Hmm::new(
            vec![0.6, 0.3, 0.1],
            vec![vec![0.7, 0.2, 0.1], vec![0.3, 0.5, 0.2], vec![0.0, 0.4, 0.6]],
            vec![vec![0.5, 0.4, 0.1], vec![0.1, 0.3, 0.6], vec![0.3, 0.3, 0.4]],
        )
        .unwrap()
    }

    /// Every path: (path, joint probability).
    fn all_paths(h: &Hmm, obs: &[usize]) -> Vec<(Vec<usize>, f64)> {
        let n = h.n_states();
        let mut out = vec![];
        let total = n.pow(obs.len() as u32);
        for code in 0..total {
            let mut path = vec![0; obs.len()];
            let mut c = code;
            for t in 0..obs.len() {
                path[t] = c % n;
                c /= n;
            }
            let mut p = h.start[path[0]] * h.emit[path[0]][obs[0]];
            for t in 1..obs.len() {
                p *= h.trans[path[t - 1]][path[t]] * h.emit[path[t]][obs[t]];
            }
            out.push((path, p));
        }
        out
    }

    struct Lcg(u64);
    impl Lcg {
        fn f(&mut self) -> f64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
        fn pick(&mut self, p: &[f64]) -> usize {
            let mut u = self.f();
            for (i, &x) in p.iter().enumerate() {
                if u < x {
                    return i;
                }
                u -= x;
            }
            p.len() - 1
        }
    }

    fn sample(h: &Hmm, len: usize, rng: &mut Lcg) -> (Vec<usize>, Vec<usize>) {
        let mut s = rng.pick(&h.start);
        let (mut states, mut obs) = (vec![], vec![]);
        for _ in 0..len {
            states.push(s);
            obs.push(rng.pick(&h.emit[s]));
            s = rng.pick(&h.trans[s]);
        }
        (states, obs)
    }

    #[test]
    fn validation() {
        assert_eq!(Hmm::new(vec![], vec![], vec![]), Err(HmmError::Shape));
        assert_eq!(
            Hmm::new(vec![1.0], vec![vec![1.0, 0.0]], vec![vec![1.0]]),
            Err(HmmError::Shape)
        );
        assert_eq!(
            Hmm::new(vec![0.5, 0.5], vec![vec![0.5, 0.5], vec![0.2, 0.7]], vec![vec![1.0], vec![1.0]]),
            Err(HmmError::NotStochastic { which: "trans", row: 1 })
        );
        assert_eq!(
            Hmm::new(vec![1.5, -0.5], vec![vec![0.5, 0.5]; 2], vec![vec![1.0]; 2]),
            Err(HmmError::Negative)
        );
        assert_eq!(
            Hmm::new(vec![0.9], vec![vec![1.0]], vec![vec![1.0]]),
            Err(HmmError::NotStochastic { which: "start", row: 0 })
        );
    }

    #[test]
    fn likelihood_matches_brute_force() {
        let mut rng = Lcg(1);
        for h in [casino(), three_state()] {
            for len in 1..=7 {
                let (_, obs) = sample(&h, len, &mut rng);
                let p: f64 = all_paths(&h, &obs).iter().map(|x| x.1).sum();
                assert!((h.log_likelihood(&obs) - p.ln()).abs() < 1e-10, "len {len}");
            }
        }
        assert_eq!(casino().log_likelihood(&[]), 0.0);
    }

    #[test]
    fn viterbi_matches_brute_force() {
        let mut rng = Lcg(2);
        for h in [casino(), three_state()] {
            for len in 1..=7 {
                let (_, obs) = sample(&h, len, &mut rng);
                let paths = all_paths(&h, &obs);
                let best = paths.iter().map(|x| x.1).fold(0.0, f64::max);
                let (vp, lp) = h.viterbi(&obs);
                assert!((lp - best.ln()).abs() < 1e-10);
                let p_of_vp = paths.iter().find(|x| x.0 == vp).unwrap().1;
                assert!((p_of_vp - best).abs() < 1e-15 * best.max(1e-300) + 1e-300);
            }
        }
        assert_eq!(casino().viterbi(&[]), (vec![], 0.0));
    }

    #[test]
    fn posterior_matches_brute_force() {
        let h = three_state();
        let obs = [0, 2, 2, 1, 0, 2];
        let paths = all_paths(&h, &obs);
        let total: f64 = paths.iter().map(|x| x.1).sum();
        let g = h.posterior(&obs);
        for t in 0..obs.len() {
            assert!((g[t].iter().sum::<f64>() - 1.0).abs() < 1e-12);
            for s in 0..3 {
                let p: f64 = paths.iter().filter(|x| x.0[t] == s).map(|x| x.1).sum::<f64>() / total;
                assert!((g[t][s] - p).abs() < 1e-12, "t={t} s={s}");
            }
        }
    }

    #[test]
    fn long_sequences_do_not_underflow() {
        let h = casino();
        let (states, obs) = sample(&h, 5000, &mut Lcg(3));
        let ll = h.log_likelihood(&obs);
        assert!(ll.is_finite() && ll < -5000.0, "{ll}");
        let (path, lp) = h.viterbi(&obs);
        assert!(lp.is_finite() && lp < ll);
        // decoding recovers most hidden states (the casino is a hard problem: ~80%)
        let acc = path.iter().zip(&states).filter(|(a, b)| a == b).count() as f64 / 5000.0;
        assert!(acc > 0.7, "viterbi accuracy {acc}");
        let g = h.posterior(&obs);
        let acc = g.iter().zip(&states).filter(|(p, s)| (p[1] > 0.5) as usize == **s).count() as f64 / 5000.0;
        assert!(acc > 0.7, "posterior accuracy {acc}");
    }

    #[test]
    fn baum_welch_increases_likelihood_and_learns() {
        let truth = casino();
        let mut rng = Lcg(4);
        let seqs: Vec<Vec<usize>> = (0..20).map(|_| sample(&truth, 500, &mut rng).1).collect();
        let init = Hmm::new(
            vec![0.5, 0.5],
            vec![vec![0.8, 0.2], vec![0.2, 0.8]],
            vec![vec![0.2, 0.15, 0.15, 0.15, 0.15, 0.2], vec![0.1, 0.15, 0.15, 0.15, 0.15, 0.3]],
        )
        .unwrap();
        let (model, hist) = init.train(&seqs, 200, 1e-6);
        assert!(hist.len() >= 3);
        for w in hist.windows(2) {
            assert!(w[1] >= w[0] - 1e-6, "likelihood decreased: {w:?}");
        }
        // the trained model is a valid HMM
        assert!(Hmm::new(model.start.clone(), model.trans.clone(), model.emit.clone()).is_ok());
        // the loaded state's six is far more likely than the fair one's
        let six_loaded = model.emit[1][5].max(model.emit[0][5]);
        assert!((six_loaded - 0.5).abs() < 0.1, "{model:?}");
    }

    #[test]
    fn baum_welch_ignores_empty_and_keeps_unused_rows() {
        let h = Hmm::new(
            vec![1.0, 0.0],
            vec![vec![1.0, 0.0], vec![0.5, 0.5]],
            vec![vec![0.5, 0.5], vec![0.9, 0.1]],
        )
        .unwrap();
        let next = h.baum_welch_step(&[vec![0, 1, 1, 0], vec![]]);
        assert_eq!(next.trans[1], vec![0.5, 0.5]); // state 1 never visited
        assert_eq!(next.emit[1], vec![0.9, 0.1]);
        assert_eq!(next.emit[0], vec![0.5, 0.5]);
        assert_eq!(h.baum_welch_step(&[]), h);
    }

    #[test]
    fn finds_a_planted_cpg_island() {
        let m = cpg_model();
        let mut rng = Lcg(5);
        let mut seq = vec![];
        for (state, len) in [(0, 3000), (1, 800), (0, 3000)] {
            for _ in 0..len {
                seq.push(b"ACGT"[rng.pick(&m.emit[state])]);
            }
        }
        seq[100] = b'N';
        let islands = cpg_islands(&seq, 100);
        assert_eq!(islands.len(), 1, "{islands:?}");
        let (s, e) = islands[0];
        assert!(s.abs_diff(3000) < 150 && e.abs_diff(3800) < 150, "{islands:?}");
    }
}
