//! # Exercise 07 — A tiny RNA-seq expression toolkit
//!
//! **Concepts:** iterator chains (`zip`, `map`, `fold`, `enumerate`,
//! `filter`), closures as arguments (`impl Fn`), returning
//! `impl Iterator + '_`, row-major matrix layout, careful handling of ties
//! and degenerate cases.
//!
//! You'll implement the core of a typical bulk RNA-seq exploratory
//! analysis, roughly what you'd do with pandas/numpy/scipy in Python:
//!
//! 1. parse a TSV count matrix (genes x samples)
//! 2. normalise to log2 counts-per-million
//! 3. filter genes with a closure
//! 4. z-score genes
//! 5. Pearson & Spearman correlation (Spearman needs *average ranks for ties*,
//!    like `scipy.stats.rankdata`)
//! 6. a gene-gene correlation matrix and the top-N most variable genes
//!
//! Try to write the numeric parts as iterator chains rather than index
//! loops. It's good practice, and the optimiser makes them just as fast.
//!
//! ```text
//! cargo test -p ex07_expression_matrix
//! ```

#![allow(unused_variables, dead_code, unused_imports, unused_mut)]

use std::collections::HashSet;

#[derive(Debug, Clone, PartialEq)]
pub enum MatrixError {
    /// No header line, or a header with no samples.
    Empty,
    /// 1-based line number; `expected`/`found` count *values* (not the gene name).
    Ragged { line: usize, expected: usize, found: usize },
    /// 1-based line and 1-based column (the gene name is column 1).
    BadNumber { line: usize, col: usize },
    DuplicateGene(String),
}

#[derive(Debug, Clone, PartialEq)]
pub struct ExprMatrix {
    pub genes: Vec<String>,
    pub samples: Vec<String>,
    /// Row-major: value for (gene g, sample s) is `data[g * n_samples + s]`.
    pub data: Vec<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    Pearson,
    Spearman,
}

impl ExprMatrix {
    /// Header: `<anything>\tS1\tS2...`; rows: `GENE\tv1\tv2...`.
    /// Fields are tab-separated and trimmed. Blank lines are skipped.
    /// A row with zero genes is allowed (just a header).
    pub fn from_tsv(text: &str) -> Result<Self, MatrixError> {
        todo!()
    }

    pub fn n_genes(&self) -> usize {
        self.genes.len()
    }

    pub fn n_samples(&self) -> usize {
        self.samples.len()
    }

    pub fn row(&self, g: usize) -> &[f64] {
        todo!()
    }

    /// Lazily iterate over one sample's values (a strided column).
    pub fn column(&self, s: usize) -> impl Iterator<Item = f64> + '_ {
        if true {
            todo!()
        }
        std::iter::empty() // placeholder so the signature type-checks
    }

    /// `log2(count / library_size * 1e6 + pseudocount)` where `library_size`
    /// is that sample's column sum. A sample with library size 0 gets CPM 0.
    pub fn log2_cpm(&self, pseudocount: f64) -> ExprMatrix {
        todo!()
    }

    /// Keep genes for which `keep(name, row)` is true (original order).
    pub fn filter_genes<F>(&self, keep: F) -> ExprMatrix
    where
        F: Fn(&str, &[f64]) -> bool,
    {
        todo!()
    }

    /// In place: each gene row -> (x - mean) / sample_std (n - 1).
    /// Rows with zero variance (or a single sample) become all zeros.
    pub fn zscore_rows(&mut self) {
        todo!()
    }

    /// Gene-gene correlation, `n_genes x n_genes`, row-major. Undefined
    /// correlations (constant gene) are `f64::NAN`, *including* on the
    /// diagonal. Symmetric. Compute each pair once.
    pub fn correlation_matrix(&self, method: Method) -> Vec<f64> {
        todo!()
    }

    /// Names of the `n` genes with the largest sample variance, in
    /// descending variance; ties broken by gene name ascending.
    pub fn top_variable(&self, n: usize) -> Vec<&str> {
        todo!()
    }
}

/// Pearson correlation. `None` if lengths differ, fewer than 2 points,
/// or either input has zero variance.
pub fn pearson(a: &[f64], b: &[f64]) -> Option<f64> {
    todo!()
}

/// 1-based ranks; tied values get the average of the ranks they span.
/// `[10, 20, 20, 5]` -> `[2.0, 3.5, 3.5, 1.0]`. Input has no NaNs.
pub fn rank_average(xs: &[f64]) -> Vec<f64> {
    todo!()
}

pub fn spearman(a: &[f64], b: &[f64]) -> Option<f64> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    const TSV: &str = "gene\tctrl1\tctrl2\ttreat1\ttreat2
BRCA1\t10\t20\t30\t40
TP53\t100\t100\t100\t100
GAPDH\t500\t480\t520\t510
MYC\t0\t5\t50\t500

ACTB\t1\t1\t2\t2
";

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn parse_ok() {
        let m = ExprMatrix::from_tsv(TSV).unwrap();
        assert_eq!(m.n_genes(), 5);
        assert_eq!(m.samples, vec!["ctrl1", "ctrl2", "treat1", "treat2"]);
        assert_eq!(m.row(3), &[0.0, 5.0, 50.0, 500.0]);
        assert_eq!(m.column(1).collect::<Vec<_>>(), vec![20.0, 100.0, 480.0, 5.0, 1.0]);
    }

    #[test]
    fn parse_errors() {
        assert_eq!(ExprMatrix::from_tsv(""), Err(MatrixError::Empty));
        assert_eq!(ExprMatrix::from_tsv("gene\n"), Err(MatrixError::Empty));
        assert_eq!(
            ExprMatrix::from_tsv("g\ta\tb\nX\t1\n"),
            Err(MatrixError::Ragged { line: 2, expected: 2, found: 1 })
        );
        assert_eq!(
            ExprMatrix::from_tsv("g\ta\tb\n\nX\t1\tfoo\n"),
            Err(MatrixError::BadNumber { line: 3, col: 3 })
        );
        assert_eq!(
            ExprMatrix::from_tsv("g\ta\nX\t1\nX\t2\n"),
            Err(MatrixError::DuplicateGene("X".into()))
        );
        assert_eq!(ExprMatrix::from_tsv("g\ta\n").unwrap().n_genes(), 0);
    }

    #[test]
    fn cpm() {
        let m = ExprMatrix::from_tsv("g\ta\tb\nX\t1\t0\nY\t3\t0\n").unwrap();
        let l = m.log2_cpm(1.0);
        assert!(close(l.row(0)[0], (250_000.0f64 + 1.0).log2()));
        assert!(close(l.row(1)[0], (750_000.0f64 + 1.0).log2()));
        assert_eq!(l.row(0)[1], 0.0); // empty library -> log2(0 + 1)
    }

    #[test]
    fn filter_with_closure() {
        let m = ExprMatrix::from_tsv(TSV).unwrap();
        let min_total = 50.0;
        let f = m.filter_genes(|name, row| !name.starts_with("GAP") && row.iter().sum::<f64>() >= min_total);
        assert_eq!(f.genes, vec!["BRCA1", "TP53", "MYC"]);
        assert_eq!(f.row(2), m.row(3));
    }

    #[test]
    fn zscores() {
        let mut m = ExprMatrix::from_tsv(TSV).unwrap();
        m.zscore_rows();
        let r = m.row(0); // 10,20,30,40: mean 25, sd sqrt(500/3)
        let sd = (500.0f64 / 3.0).sqrt();
        assert!(close(r[0], -15.0 / sd) && close(r[3], 15.0 / sd));
        assert_eq!(m.row(1), &[0.0; 4]);
    }

    #[test]
    fn pearson_and_spearman() {
        assert!(close(pearson(&[1.0, 2.0, 3.0], &[2.0, 4.0, 6.0]).unwrap(), 1.0));
        assert!(close(pearson(&[1.0, 2.0, 3.0], &[3.0, 2.0, 1.0]).unwrap(), -1.0));
        assert_eq!(pearson(&[1.0, 1.0], &[1.0, 2.0]), None);
        assert_eq!(pearson(&[1.0], &[1.0]), None);
        assert_eq!(pearson(&[1.0, 2.0], &[1.0]), None);
        // scipy.stats.pearsonr([1,2,3,4,5],[2,1,4,3,7]) -> 0.8241633836921339
        assert!(close(pearson(&[1., 2., 3., 4., 5.], &[2., 1., 4., 3., 7.]).unwrap(), 0.8241633836921339));
        // monotone but non-linear -> spearman 1
        assert!(close(spearman(&[1., 2., 3., 4.], &[1., 8., 27., 1000.]).unwrap(), 1.0));
        // scipy.stats.spearmanr([1,2,2,3,4],[5,6,7,8,7]) -> 0.7631578947368421
        assert!(close(spearman(&[1., 2., 2., 3., 4.], &[5., 6., 7., 8., 7.]).unwrap(), 0.7631578947368421));
    }

    #[test]
    fn ranks_with_ties() {
        assert_eq!(rank_average(&[10.0, 20.0, 20.0, 5.0]), vec![2.0, 3.5, 3.5, 1.0]);
        assert_eq!(rank_average(&[7.0, 7.0, 7.0]), vec![2.0, 2.0, 2.0]);
        assert!(rank_average(&[]).is_empty());
    }

    #[test]
    fn corr_matrix() {
        let m = ExprMatrix::from_tsv(TSV).unwrap();
        let c = m.correlation_matrix(Method::Spearman);
        let g = m.n_genes();
        assert_eq!(c.len(), g * g);
        assert!(close(c[3], 1.0)); // BRCA1 vs MYC monotone
        assert!(c[g + 1].is_nan()); // TP53 constant
        assert!(c[g].is_nan());
        for i in 0..g {
            for j in 0..g {
                let (a, b) = (c[i * g + j], c[j * g + i]);
                assert!(a == b || (a.is_nan() && b.is_nan()));
            }
        }
        let p = m.correlation_matrix(Method::Pearson);
        assert!(close(p[0], 1.0));
        assert!(p[3] < 0.9); // MYC is not linear in BRCA1
    }

    #[test]
    fn most_variable() {
        let m = ExprMatrix::from_tsv(TSV).unwrap();
        assert_eq!(m.top_variable(3), vec!["MYC", "GAPDH", "BRCA1"]);
        let t = ExprMatrix::from_tsv("g\ta\tb\nZ\t1\t2\nA\t1\t2\nM\t0\t0\n").unwrap();
        assert_eq!(t.top_variable(10), vec!["A", "Z", "M"]);
    }
}
