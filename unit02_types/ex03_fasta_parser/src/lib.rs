//! # Exercise 03 — A strict FASTA parser
//!
//! **Concepts:** structs, enums with data, `Result` and custom error types,
//! `impl Display` + `impl std::error::Error`, `FromStr`, line-oriented
//! parsing with good error positions.
//!
//! Real-world FASTA files are messy, and a good parser reports *where* things
//! went wrong. Implement `parse_fasta` according to these rules:
//!
//! * Lines are numbered from **1**. Columns are numbered from **1** and count
//!   characters in the line (after removing a trailing `\r`).
//! * Trailing whitespace on any line is ignored. Lines that are empty or
//!   whitespace-only are skipped everywhere.
//! * A header line starts with `>`. The **id** is the text immediately after
//!   `>` up to the first whitespace (so `>` followed directly by whitespace or
//!   by nothing means an empty id). The **description** is the rest of the
//!   line, trimmed, or `None` if that is empty.
//! * A sequence may span many lines; the lines are concatenated.
//!   Allowed sequence characters are ASCII letters, `*` and `-`. They are
//!   stored **uppercased**. Anything else (including leading spaces or
//!   internal spaces) is `InvalidChar` at the first offending column.
//! * Errors are reported in the order a single top-to-bottom pass would hit
//!   them. When a header line is reached, check in this order:
//!   1. the *previous* record has an empty sequence -> `EmptySequence`
//!   2. the new id is empty -> `EmptyId`
//!   3. the new id was already used -> `DuplicateId`
//! * At end of input, the last record is checked for `EmptySequence`.
//! * Input with no non-blank lines is `EmptyInput`.
//! * Sequence data before any header is `MissingHeader`.
//!
//! Then implement `write_fasta` (line wrapping) and `FromStr` for `Record`
//! (parse exactly one record).
//!
//! ```text
//! cargo test -p ex03_fasta_parser
//! ```

#![allow(unused_variables, dead_code, unused_imports, unused_mut)]

use std::collections::HashSet;
use std::fmt;
use std::str::FromStr;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub id: String,
    pub description: Option<String>,
    pub seq: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FastaError {
    EmptyInput,
    MissingHeader { line: usize },
    EmptyId { line: usize },
    InvalidChar { line: usize, col: usize, ch: char },
    EmptySequence { id: String },
    DuplicateId { id: String },
    /// Only used by `FromStr for Record`: the input held a number of
    /// records other than exactly one.
    NotExactlyOne { found: usize },
}

impl fmt::Display for FastaError {
    /// Human-readable messages. Tests only check that each message contains
    /// the relevant numbers / ids, so phrase them however you like.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        todo!()
    }
}

impl std::error::Error for FastaError {}

pub fn parse_fasta(input: &str) -> Result<Vec<Record>, FastaError> {
    todo!()
}

/// Serialize records. Header is `>id` or `>id description`. Sequence lines
/// are wrapped to at most `width` characters (`width == 0` means no
/// wrapping). Every line, including the last, ends with `\n`.
pub fn write_fasta(records: &[Record], width: usize) -> String {
    todo!()
}

impl FromStr for Record {
    type Err = FastaError;

    /// Parse text containing **exactly one** record.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_multi_record_with_wrapping_and_crlf() {
        let input = ">sp|P1 Insulin  human \r\nMALW\r\nmrll\r\n\r\n>seq2\nAC-GT*\n";
        let recs = parse_fasta(input).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, "sp|P1");
        assert_eq!(recs[0].description.as_deref(), Some("Insulin  human"));
        assert_eq!(recs[0].seq, "MALWMRLL");
        assert_eq!(recs[1].id, "seq2");
        assert_eq!(recs[1].description, None);
        assert_eq!(recs[1].seq, "AC-GT*");
    }

    #[test]
    fn empty_inputs() {
        assert_eq!(parse_fasta(""), Err(FastaError::EmptyInput));
        assert_eq!(parse_fasta("\n   \n\t\n"), Err(FastaError::EmptyInput));
    }

    #[test]
    fn missing_header() {
        assert_eq!(parse_fasta("\n\nACGT\n>x\nA"), Err(FastaError::MissingHeader { line: 3 }));
    }

    #[test]
    fn empty_id() {
        assert_eq!(parse_fasta(">a\nA\n>   desc only\nC"), Err(FastaError::EmptyId { line: 3 }));
    }

    #[test]
    fn invalid_char_position() {
        assert_eq!(
            parse_fasta(">a\nACGT\nAC GT\n"),
            Err(FastaError::InvalidChar { line: 3, col: 3, ch: ' ' })
        );
        assert_eq!(
            parse_fasta(">a\n  ACGT\n"),
            Err(FastaError::InvalidChar { line: 2, col: 1, ch: ' ' })
        );
        assert_eq!(
            parse_fasta(">a\nAC1\n"),
            Err(FastaError::InvalidChar { line: 2, col: 3, ch: '1' })
        );
    }

    #[test]
    fn error_ordering() {
        // previous record empty is reported before the duplicate id
        assert_eq!(
            parse_fasta(">a\n>a\nACGT"),
            Err(FastaError::EmptySequence { id: "a".into() })
        );
        assert_eq!(
            parse_fasta(">a\nA\n>b\nC\n>a\nG"),
            Err(FastaError::DuplicateId { id: "a".into() })
        );
        assert_eq!(
            parse_fasta(">a\nA\n>last\n\n  \n"),
            Err(FastaError::EmptySequence { id: "last".into() })
        );
    }

    #[test]
    fn display_messages_carry_context() {
        let e = FastaError::InvalidChar { line: 12, col: 7, ch: '#' };
        let s = e.to_string();
        assert!(s.contains("12") && s.contains('7'), "{s}");
        let e: Box<dyn std::error::Error> = Box::new(FastaError::DuplicateId { id: "chrX".into() });
        assert!(e.to_string().contains("chrX"));
    }

    #[test]
    fn roundtrip_and_wrap() {
        let recs = vec![
            Record { id: "a".into(), description: Some("first one".into()), seq: "ACGTACGTAC".into() },
            Record { id: "b".into(), description: None, seq: "MKV".into() },
        ];
        let out = write_fasta(&recs, 4);
        assert_eq!(out, ">a first one\nACGT\nACGT\nAC\n>b\nMKV\n");
        assert_eq!(parse_fasta(&out).unwrap(), recs);
        assert_eq!(write_fasta(&recs[1..], 0), ">b\nMKV\n");
    }

    #[test]
    fn from_str_exactly_one() {
        let r: Record = ">x hello\nacgt\n".parse().unwrap();
        assert_eq!(r.seq, "ACGT");
        let err = ">x\nA\n>y\nC".parse::<Record>().unwrap_err();
        assert_eq!(err, FastaError::NotExactlyOne { found: 2 });
        assert_eq!("".parse::<Record>().unwrap_err(), FastaError::EmptyInput);
    }
}
