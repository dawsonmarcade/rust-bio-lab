//! # Exercise 02 — FIX IT: sequencing-read QC
//!
//! **Concepts:** moves vs borrows, `iter()` vs `iter_mut()`, borrowing while
//! mutating, returning references to locals, `move` closures and threads.
//!
//! This crate **does not compile**, and once it does, one test still fails
//! because of a logic bug. Your job:
//!
//! 1. Make it compile. Read every compiler error carefully: `rustc` usually
//!    tells you exactly what is wrong. Fix the *root cause*, not the symptom.
//!    Sprinkling `.clone()` everywhere is not allowed, and neither is
//!    changing the tests. Signatures (return types, `self` vs `&self`,
//!    `&` vs `&mut`) may change where that is the real fix.
//! 2. Make all tests pass. There is one silent logic bug that compiles
//!    fine but gives wrong answers. Read the FASTQ spec in the docs.
//!
//! ```text
//! cargo test -p ex02_fix_read_qc
//! ```
//!
//! There are 6 compile problems and 1 logic bug.

#![allow(dead_code)]

#[derive(Debug, Clone, PartialEq)]
pub struct Read {
    pub id: String,
    pub seq: String,
    /// Phred quality scores (already decoded, i.e. 0..=41ish, NOT ASCII).
    pub qual: Vec<u8>,
}

impl Read {
    pub fn mean_quality(self) -> f64 {
        if self.qual.is_empty() {
            return 0.0;
        }
        self.qual.iter().map(|&q| q as f64).sum::<f64>() / self.qual.len() as f64
    }
}

/// Parse 4-line FASTQ records. Quality characters use the Sanger / Illumina
/// 1.8+ encoding: `phred = ascii_code - 33`.
pub fn parse_fastq(text: &str) -> Vec<Read> {
    let lines: Vec<&str> = text.lines().filter(|l| !l.trim().is_empty()).collect();
    lines
        .chunks(4)
        .filter(|c| c.len() == 4)
        .map(|c| Read {
            id: c[0].trim_start_matches('@').to_string(),
            seq: c[1].to_string(),
            qual: c[3].bytes().map(|b| b.saturating_sub(64)).collect(),
        })
        .collect()
}

/// First `k` bases of `seq`, uppercased.
pub fn first_kmer_upper(seq: &str, k: usize) -> &str {
    let up = seq.to_uppercase();
    &up[..k]
}

pub fn read_ids(reads: &[Read]) -> Vec<&str> {
    reads.iter().map(|r| r.id).collect()
}

/// Truncate each read at the first base whose quality is `< min_q`, then
/// drop reads that became empty.
pub fn trim_reads(reads: &mut Vec<Read>, min_q: u8) {
    for r in reads.iter() {
        if let Some(cut) = r.qual.iter().position(|&q| q < min_q) {
            r.seq.truncate(cut);
            r.qual.truncate(cut);
        }
    }
    reads.retain(|r| !r.seq.is_empty());
}

/// Returns the id of the read with the highest mean quality *among the reads
/// present before* `extra` is appended, then appends `extra`.
pub fn best_read_then_push(reads: &mut Vec<Read>, extra: Read) -> String {
    let best = reads
        .iter()
        .max_by(|a, b| a.mean_quality().total_cmp(&b.mean_quality()))
        .expect("non-empty");
    reads.push(extra);
    best.id.clone()
}

/// Count total bases on a background thread.
pub fn count_bases_in_background(reads: Vec<Read>) -> usize {
    let handle = std::thread::spawn(|| reads.iter().map(|r| r.seq.len()).sum());
    handle.join().unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    const FASTQ: &str = "@r1
ACGTACGT
+
IIIIII#I
@r2
TTTTGGGG
+
55555555
@r3
CCCC
+
!!!!
";

    #[test]
    fn parses_sanger_qualities() {
        let reads = parse_fastq(FASTQ);
        assert_eq!(reads.len(), 3);
        assert_eq!(reads[0].id, "r1");
        assert_eq!(reads[0].qual, vec![40, 40, 40, 40, 40, 40, 2, 40]);
        assert_eq!(reads[1].qual, vec![20; 8]);
        assert_eq!(reads[2].qual, vec![0; 4]);
    }

    #[test]
    fn mean_quality_does_not_consume() {
        let reads = parse_fastq(FASTQ);
        let r = &reads[1];
        assert!((r.mean_quality() - 20.0).abs() < 1e-12);
        assert_eq!(r.id, "r2"); // still usable afterwards
    }

    #[test]
    fn first_kmer() {
        assert_eq!(first_kmer_upper("acgtNN", 3), "ACG");
    }

    #[test]
    fn ids_borrowed() {
        let reads = parse_fastq(FASTQ);
        assert_eq!(read_ids(&reads), vec!["r1", "r2", "r3"]);
    }

    #[test]
    fn trimming() {
        let mut reads = parse_fastq(FASTQ);
        trim_reads(&mut reads, 10);
        assert_eq!(reads.len(), 2);
        assert_eq!(reads[0].seq, "ACGTAC");
        assert_eq!(reads[1].seq, "TTTTGGGG");
    }

    #[test]
    fn best_then_push() {
        let mut reads = parse_fastq(FASTQ);
        let extra = Read { id: "perfect".into(), seq: "A".into(), qual: vec![60] };
        assert_eq!(best_read_then_push(&mut reads, extra), "r1");
        assert_eq!(reads.len(), 4);
        assert_eq!(reads[3].id, "perfect");
    }

    #[test]
    fn background_count() {
        let reads = parse_fastq(FASTQ);
        assert_eq!(count_bases_in_background(reads), 20);
    }
}
