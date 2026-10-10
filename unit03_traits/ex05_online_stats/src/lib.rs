//! # Exercise 05 — Streaming statistics with traits
//!
//! **Concepts:** traits with associated types, generic functions with trait
//! bounds, `Default`, implementing a trait for tuples (composition),
//! extension traits on `[T]`, numerically stable algorithms.
//!
//! When you stream millions of measurements (read qualities, expression
//! values, MD energies), you can't keep them all in memory. You compute
//! statistics *online*. And when you split the work across threads (unit
//! 9), you need to *merge* partial results.
//!
//! * `Mean`, `Variance` (Welford's algorithm), `MinMax` implement
//!   `OnlineStat`.
//! * `(A, B)` implements `OnlineStat` whenever `A` and `B` do, so you can
//!   compute several statistics in a single pass.
//! * `Mergeable` combines two partial results. For variance, look up
//!   *Chan et al.'s parallel algorithm*. The test checks against a dataset
//!   where the naive `E[x²] - E[x]²` formula is badly wrong.
//! * `SliceStats` adds `.median()` / `.quantile()` to any slice of numbers.
//!
//! ```text
//! cargo test -p ex05_online_stats
//! ```

#![allow(unused_variables, dead_code, unused_imports)]

pub trait OnlineStat {
    type Output;
    fn push(&mut self, x: f64);
    fn finish(&self) -> Self::Output;
}

pub trait Mergeable {
    /// After `a.merge(b)`, `a` must equal what you'd get by pushing all of
    /// `a`'s and `b`'s values into one accumulator.
    fn merge(&mut self, other: Self);
}

/// Arithmetic mean. `finish` -> `None` when empty.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Mean {
    // your fields here
}

/// **Sample** variance (divide by n - 1) via Welford. `None` if n < 2.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Variance {
    // your fields here
}

/// `(min, max)`, `None` when empty. NaN inputs are ignored.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MinMax {
    // your fields here
}

impl OnlineStat for Mean {
    type Output = Option<f64>;
    fn push(&mut self, x: f64) {
        todo!()
    }
    fn finish(&self) -> Option<f64> {
        todo!()
    }
}

impl Mergeable for Mean {
    fn merge(&mut self, other: Self) {
        todo!()
    }
}

impl OnlineStat for Variance {
    type Output = Option<f64>;
    fn push(&mut self, x: f64) {
        todo!()
    }
    fn finish(&self) -> Option<f64> {
        todo!()
    }
}

impl Mergeable for Variance {
    fn merge(&mut self, other: Self) {
        todo!()
    }
}

impl OnlineStat for MinMax {
    type Output = Option<(f64, f64)>;
    fn push(&mut self, x: f64) {
        todo!()
    }
    fn finish(&self) -> Option<(f64, f64)> {
        todo!()
    }
}

impl Mergeable for MinMax {
    fn merge(&mut self, other: Self) {
        todo!()
    }
}

impl<A: OnlineStat, B: OnlineStat> OnlineStat for (A, B) {
    type Output = (A::Output, B::Output);
    fn push(&mut self, x: f64) {
        todo!()
    }
    fn finish(&self) -> Self::Output {
        todo!()
    }
}

impl<A: Mergeable, B: Mergeable> Mergeable for (A, B) {
    fn merge(&mut self, other: Self) {
        todo!()
    }
}

/// Feed every item of `iter` into a fresh `S` and return its output.
/// Must work for `Vec<f64>`, `&[f32]` iterators of copied values, `u8`s, ...
pub fn collect_stat<S, I>(iter: I) -> S::Output
where
    S: OnlineStat + Default,
    I: IntoIterator,
    I::Item: Into<f64>,
{
    todo!()
}

/// Index of the first maximum. Items that are not comparable with the
/// current max (e.g. NaN) are skipped. `None` for an empty iterator or if
/// nothing was comparable.
pub fn argmax<T, I>(iter: I) -> Option<usize>
where
    T: PartialOrd,
    I: IntoIterator<Item = T>,
{
    todo!()
}

pub trait SliceStats {
    fn mean(&self) -> Option<f64>;
    /// `None` if empty or any value is NaN.
    fn median(&self) -> Option<f64>;
    /// Linear interpolation, same as `numpy.quantile(..., method="linear")`:
    /// sort, `h = (n - 1) * q`, interpolate between `x[floor(h)]` and
    /// `x[floor(h) + 1]`. `None` if empty, any NaN, or `q` outside `[0, 1]`.
    fn quantile(&self, q: f64) -> Option<f64>;
}

impl<T: Copy + Into<f64>> SliceStats for [T] {
    fn mean(&self) -> Option<f64> {
        todo!()
    }
    fn median(&self) -> Option<f64> {
        todo!()
    }
    fn quantile(&self, q: f64) -> Option<f64> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * (1.0 + b.abs())
    }

    #[test]
    fn mean_and_variance() {
        let xs = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        assert_eq!(collect_stat::<Mean, _>(xs), Some(5.0));
        let v = collect_stat::<Variance, _>(xs).unwrap();
        assert!(close(v, 32.0 / 7.0, 1e-12));
        assert_eq!(collect_stat::<Variance, _>([1.0]), None);
        assert_eq!(collect_stat::<Mean, _>(Vec::<f64>::new()), None);
    }

    #[test]
    fn generic_inputs() {
        let quals: Vec<u8> = vec![30, 40, 20];
        assert_eq!(collect_stat::<Mean, _>(quals.iter().copied()), Some(30.0));
        let f: [f32; 3] = [1.0, 2.0, 3.0];
        assert_eq!(collect_stat::<MinMax, _>(f), Some((1.0, 3.0)));
    }

    #[test]
    fn minmax_ignores_nan() {
        assert_eq!(collect_stat::<MinMax, _>([f64::NAN, 3.0, -1.0, f64::NAN]), Some((-1.0, 3.0)));
        assert_eq!(collect_stat::<MinMax, _>([f64::NAN]), None);
    }

    #[test]
    fn tuples_compose() {
        let ((mean, var), mm) =
            collect_stat::<((Mean, Variance), MinMax), _>([1.0, 2.0, 3.0, 4.0]);
        assert_eq!(mean, Some(2.5));
        assert!(close(var.unwrap(), 5.0 / 3.0, 1e-12));
        assert_eq!(mm, Some((1.0, 4.0)));
    }

    #[test]
    fn welford_is_stable() {
        // Huge offset: the naive E[x^2]-E[x]^2 formula returns garbage here.
        let xs: Vec<f64> = (0..1000).map(|i| 1e9 + (i % 10) as f64).collect();
        let v = collect_stat::<Variance, _>(xs.iter().copied()).unwrap();
        let exact = {
            let m = 4.5;
            xs.iter().map(|x| (x - 1e9 - m).powi(2)).sum::<f64>() / 999.0
        };
        assert!(close(v, exact, 1e-6), "{v} vs {exact}");
    }

    #[test]
    fn merging_matches_single_pass() {
        let xs: Vec<f64> = (0..997).map(|i| ((i * 7919) % 1013) as f64 * 0.37 + 1e6).collect();
        let whole = {
            let mut s = (Variance::default(), (Mean::default(), MinMax::default()));
            xs.iter().for_each(|&x| s.push(x));
            s.finish()
        };
        for split in [0, 1, 13, 500, 996, 997] {
            let mut a = (Variance::default(), (Mean::default(), MinMax::default()));
            let mut b = a.clone();
            xs[..split].iter().for_each(|&x| a.push(x));
            xs[split..].iter().for_each(|&x| b.push(x));
            a.merge(b);
            let got = a.finish();
            assert!(close(got.0.unwrap(), whole.0.unwrap(), 1e-9), "split {split}");
            assert!(close(got.1.0.unwrap(), whole.1.0.unwrap(), 1e-12), "split {split}");
            assert_eq!(got.1.1, whole.1.1);
        }
    }

    #[test]
    fn argmax_cases() {
        assert_eq!(argmax([3, 9, 2, 9]), Some(1));
        assert_eq!(argmax([f64::NAN, 1.0, 5.0, f64::NAN]), Some(2));
        assert_eq!(argmax(Vec::<i32>::new()), None);
        assert_eq!(argmax([f64::NAN]), None);
        assert_eq!(argmax(["b", "c", "a"]), Some(1));
    }

    #[test]
    fn slice_quantiles() {
        let xs = [7.0, 1.0, 3.0, 5.0];
        assert_eq!(xs.median(), Some(4.0));
        assert_eq!(xs.quantile(0.0), Some(1.0));
        assert_eq!(xs.quantile(1.0), Some(7.0));
        assert!(close(xs.quantile(0.25).unwrap(), 2.5, 1e-12));
        assert!(close(xs.quantile(0.9).unwrap(), 6.4, 1e-12));
        let ints: &[u16] = &[5, 1, 3];
        assert_eq!(ints.median(), Some(3.0));
        assert_eq!(ints.mean(), Some(3.0));
        assert_eq!([1.0, f64::NAN].median(), None);
        assert_eq!(xs.quantile(1.5), None);
        assert_eq!(Vec::<f64>::new().median(), None);
    }
}
