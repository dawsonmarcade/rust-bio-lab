//! # Exercise 04 — Genetic code, ORFs and hydropathy
//!
//! **Concepts:** fieldless enums with derives, `TryFrom`, `Display`,
//! exhaustive `match`, enums carrying data (`Strand`), sorting with
//! `Ord`, working with `&[u8]`.
//!
//! Part A: an `AminoAcid` enum with conversions.
//! Part B: translation in a reading frame.
//! Part C: open reading frames on **both strands**, which is the tricky bit.
//! Part D: a Kyte–Doolittle hydropathy sliding window.
//!
//! The standard codon table is given as `CODON_TABLE`: a 64-byte string of
//! one-letter codes indexed by `16*b1 + 4*b2 + b3` with base order
//! **T=0, C=1, A=2, G=3** (so `TTT` -> index 0 -> `F`, `GGG` -> 63 -> `G`).
//! `*` means stop.
//!
//! ```text
//! cargo test -p ex04_translation
//! ```

#![allow(unused_variables, dead_code, unused_imports)]

use std::fmt;

pub const CODON_TABLE: &[u8; 64] =
    b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum AminoAcid {
    Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile,
    Leu, Lys, Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val,
    Stop,
}

use AminoAcid::*;

impl AminoAcid {
    /// One-letter code -> amino acid (`'*'` is `Stop`). Case-insensitive.
    pub fn from_one_letter(c: char) -> Option<Self> {
        todo!()
    }

    /// Uppercase one-letter code, `'*'` for stop.
    pub fn one_letter(self) -> char {
        todo!()
    }

    /// Kyte–Doolittle hydropathy index; `None` for `Stop`.
    ///
    /// A 1.8, R -4.5, N -3.5, D -3.5, C 2.5, Q -3.5, E -3.5, G -0.4,
    /// H -3.2, I 4.5, L 3.8, K -3.9, M 1.9, F 2.8, P -1.6, S -0.8,
    /// T -0.7, W -0.9, Y -1.3, V 4.2
    pub fn kyte_doolittle(self) -> Option<f64> {
        todo!()
    }
}

/// `Display` prints the three-letter code (`Ala`, `Arg`, ...), and `Ter` for stop.
impl fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        todo!()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CodonError {
    WrongLength(usize),
    /// The offending byte. `U`/`u` is accepted as `T`; lowercase is fine.
    InvalidBase(u8),
}

impl TryFrom<&[u8]> for AminoAcid {
    type Error = CodonError;
    fn try_from(codon: &[u8]) -> Result<Self, Self::Error> {
        todo!()
    }
}

/// Translate `seq` starting at offset `frame` (0, 1 or 2). A trailing partial
/// codon is ignored. Stops are included as `Stop` (translation does *not*
/// halt at them). Panics if `frame > 2`.
pub fn translate(seq: &[u8], frame: usize) -> Result<Vec<AminoAcid>, CodonError> {
    todo!()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Orf {
    /// 0-based, half-open interval **in forward-strand coordinates**, and
    /// covering the stop codon. For reverse-strand ORFs, `start` is where the
    /// stop codon's reverse complement begins on the forward strand.
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
    /// One-letter protein, starting with `M`, *excluding* the stop.
    pub protein: String,
}

/// Find ORFs on all 6 frames (3 forward + 3 on the reverse complement).
///
/// Input is uppercase `ACGT` only.
///
/// Within a frame, scan codon by codon. An ORF opens at the **first** `ATG`
/// after the previous stop codon (or after the frame start) and closes at
/// the next in-frame stop codon. Nested `ATG`s inside an open ORF do *not*
/// start new ORFs. An ORF that never reaches a stop is discarded. Keep ORFs
/// whose protein length (amino acids, without stop) is `>= min_aa`.
///
/// Result is sorted by `(start, end, strand)`.
pub fn find_orfs(seq: &[u8], min_aa: usize) -> Vec<Orf> {
    todo!()
}

/// Mean Kyte–Doolittle score of every full window of size `window`.
/// Returns an empty vec if `window == 0` or `window > protein.len()`.
/// Returns `None` if the protein contains `Stop`.
/// Should be O(n), not O(n * window).
pub fn hydropathy_profile(protein: &[AminoAcid], window: usize) -> Option<Vec<f64>> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn prot(s: &str) -> Vec<AminoAcid> {
        s.chars().map(|c| AminoAcid::from_one_letter(c).unwrap()).collect()
    }

    #[test]
    fn letters_roundtrip() {
        for &b in CODON_TABLE.iter() {
            let aa = AminoAcid::from_one_letter(b as char).unwrap();
            assert_eq!(aa.one_letter(), b as char);
        }
        assert_eq!(AminoAcid::from_one_letter('w'), Some(Trp));
        assert_eq!(AminoAcid::from_one_letter('B'), None);
    }

    #[test]
    fn display_three_letter() {
        assert_eq!(Trp.to_string(), "Trp");
        assert_eq!(Stop.to_string(), "Ter");
        assert_eq!(format!("{}-{}", Met, Lys), "Met-Lys");
    }

    #[test]
    fn codons() {
        assert_eq!(AminoAcid::try_from(&b"ATG"[..]), Ok(Met));
        assert_eq!(AminoAcid::try_from(&b"uga"[..]), Ok(Stop));
        assert_eq!(AminoAcid::try_from(&b"TTT"[..]), Ok(Phe));
        assert_eq!(AminoAcid::try_from(&b"GGG"[..]), Ok(Gly));
        assert_eq!(AminoAcid::try_from(&b"TGG"[..]), Ok(Trp));
        assert_eq!(AminoAcid::try_from(&b"AT"[..]), Err(CodonError::WrongLength(2)));
        assert_eq!(AminoAcid::try_from(&b"ANG"[..]), Err(CodonError::InvalidBase(b'N')));
    }

    #[test]
    fn translate_frames() {
        let s = b"ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG";
        let p: String = translate(s, 0).unwrap().iter().map(|a| a.one_letter()).collect();
        assert_eq!(p, "MAIVMGR*KGAR*");
        let p1 = translate(s, 1).unwrap();
        assert_eq!(p1.len(), 12);
        assert!(translate(b"AC", 2).unwrap().is_empty());
        assert_eq!(translate(b"ATGXXX", 0), Err(CodonError::InvalidBase(b'X')));
    }

    #[test]
    fn simple_forward_orf() {
        let orfs = find_orfs(b"ATGAAATAG", 2);
        assert_eq!(
            orfs,
            vec![Orf { start: 0, end: 9, strand: Strand::Forward, protein: "MK".into() }]
        );
        assert!(find_orfs(b"ATGAAATAG", 3).is_empty());
    }

    #[test]
    fn nested_atg_is_one_orf() {
        let orfs = find_orfs(b"ATGATGAAATAA", 1);
        assert_eq!(orfs.len(), 1);
        assert_eq!(orfs[0].protein, "MMK");
    }

    #[test]
    fn unterminated_orf_discarded() {
        assert!(find_orfs(b"ATGAAAAAAAAA", 1).is_empty());
    }

    #[test]
    fn reverse_strand_coordinates() {
        // reverse complement of ATGCCCTGA, padded on both sides
        let orfs = find_orfs(b"GGTCAGGGCATGG", 1);
        assert_eq!(
            orfs,
            vec![Orf { start: 2, end: 11, strand: Strand::Reverse, protein: "MP".into() }]
        );
    }

    #[test]
    fn rosalind_style() {
        let seq = b"AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG";
        let mut proteins: Vec<String> = find_orfs(seq, 1).into_iter().map(|o| o.protein).collect();
        proteins.sort();
        // two distinct one-residue ORFs ("M") exist at different places.
        // MTPRLGLESLLE is nested inside MGMTPRLGLESLLE, so it is NOT separate
        assert_eq!(proteins, vec!["M", "M", "MGMTPRLGLESLLE", "MLLGSFRLIPKETLIQVAGSSPCNLS"]);
        // orfs come back sorted by start
        let orfs = find_orfs(seq, 1);
        assert!(orfs.windows(2).all(|w| (w[0].start, w[0].end) <= (w[1].start, w[1].end)));
    }

    #[test]
    fn hydropathy() {
        let p = prot("IVLKKD");
        let h = hydropathy_profile(&p, 3).unwrap();
        let expect = [
            (4.5 + 4.2 + 3.8) / 3.0,
            (4.2 + 3.8 - 3.9) / 3.0,
            (3.8 - 3.9 - 3.9) / 3.0,
            (-3.9 - 3.9 - 3.5) / 3.0,
        ];
        assert_eq!(h.len(), 4);
        for (a, b) in h.iter().zip(expect) {
            assert!((a - b).abs() < 1e-9);
        }
        assert_eq!(hydropathy_profile(&p, 0), Some(vec![]));
        assert_eq!(hydropathy_profile(&p, 7), Some(vec![]));
        assert_eq!(hydropathy_profile(&prot("MK*"), 1), None);
    }
}
