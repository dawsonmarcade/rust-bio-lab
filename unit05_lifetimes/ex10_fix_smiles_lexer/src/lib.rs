//! # Exercise 10 — FIX IT: a SMILES lexer with lifetime bugs
//!
//! **Concepts:** lifetime elision rules and when they pick the *wrong*
//! lifetime, structs holding references, `'a` vs `'_`, closures capturing
//! references, returning references to temporaries.
//!
//! This crate does **not compile**. Every error is a lifetime or borrow
//! problem. Fix them without changing the tests, without `'static` hacks,
//! and without `.to_string()`/`.clone()` *inside the lexer* (the lexer
//! must stay zero-copy: tokens borrow from the source string).
//!
//! Before each fix, write down in a comment *why* the compiler rejected it.
//! That is the real exercise.
//!
//! Lexing rules (already implemented correctly):
//! * `Cl`, `Br`: two-letter atoms. Any other uppercase letter: one-letter
//!   atom. Lowercase `b c n o p s`: aromatic atom.
//! * `[ ... ]`: a bracket atom, the whole thing including brackets is one
//!   token.
//! * `- = # :` bonds, `(` `)` branches, digits ring closures, `.` dot.
//! * Anything else is skipped.
//!
//! ```text
//! cargo test -p ex10_fix_smiles_lexer
//! ```
//!
//! There are 6 problems.

#![allow(dead_code)]

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Atom,
    Bond,
    BranchOpen,
    BranchClose,
    Ring,
    Dot,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Token {
    pub text: &str,
    pub kind: Kind,
    pub pos: usize,
}

pub struct Lexer<'a> {
    src: &'a str,
    pos: usize,
}

impl<'a> Lexer<'a> {
    pub fn new(src: &'a str) -> Self {
        Lexer { src, pos: 0 }
    }

    pub fn next_token(&mut self) -> Option<Token<'_>> {
        let bytes = self.src.as_bytes();
        while self.pos < bytes.len() {
            let start = self.pos;
            let c = bytes[start];
            let (len, kind) = match c {
                b'C' if bytes.get(start + 1) == Some(&b'l') => (2, Kind::Atom),
                b'B' if bytes.get(start + 1) == Some(&b'r') => (2, Kind::Atom),
                b'A'..=b'Z' | b'b' | b'c' | b'n' | b'o' | b'p' | b's' => (1, Kind::Atom),
                b'[' => match self.src[start..].find(']') {
                    Some(end) => (end + 1, Kind::Atom),
                    None => (bytes.len() - start, Kind::Atom),
                },
                b'-' | b'=' | b'#' | b':' => (1, Kind::Bond),
                b'(' => (1, Kind::BranchOpen),
                b')' => (1, Kind::BranchClose),
                b'0'..=b'9' => (1, Kind::Ring),
                b'.' => (1, Kind::Dot),
                _ => {
                    self.pos += 1;
                    continue;
                }
            };
            self.pos += len;
            return Some(Token { text: &self.src[start..start + len], kind, pos: start });
        }
        None
    }
}

pub fn tokenize(src: &str) -> Vec<Token<'_>> {
    let mut lexer = Lexer::new(src);
    let mut out = Vec::new();
    while let Some(t) = lexer.next_token() {
        out.push(t);
    }
    out
}

/// Whichever token text is longer (first one on ties).
pub fn longer_text(a: &Token, b: &Token) -> &str {
    if a.text.len() >= b.text.len() { a.text } else { b.text }
}

/// Uppercased atom symbols (brackets kept as-is).
pub fn upper_atoms<'a>(tokens: &[Token<'a>]) -> Vec<&'a str> {
    tokens
        .iter()
        .filter(|t| t.kind == Kind::Atom)
        .map(|t| t.text.to_uppercase().as_str())
        .collect()
}

/// Returns a predicate that keeps tokens whose text starts with `prefix`.
pub fn starts_with_filter(prefix: &str) -> impl Fn(&Token) -> bool {
    |t| t.text.starts_with(prefix)
}

pub struct SymbolTable {
    names: Vec<String>,
}

impl SymbolTable {
    pub fn new() -> Self {
        SymbolTable { names: vec!["C".into()] }
    }
    /// Push `name` and return the length of the *first* stored name.
    pub fn push_and_first_len(&mut self, name: &str) -> usize {
        let first = &self.names[0];
        self.names.push(name.to_string());
        first.len()
    }
    pub fn len(&self) -> usize {
        self.names.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn texts(src: &str) -> Vec<&str> {
        tokenize(src).into_iter().map(|t| t.text).collect()
    }

    #[test]
    fn lexes_aspirin() {
        assert_eq!(
            texts("CC(=O)Oc1ccccc1C(=O)O"),
            vec!["C", "C", "(", "=", "O", ")", "O", "c", "1", "c", "c", "c", "c", "c", "1", "C", "(", "=", "O", ")", "O"]
        );
    }

    #[test]
    fn two_letter_and_brackets() {
        assert_eq!(texts("ClCBr.[NH4+]"), vec!["Cl", "C", "Br", ".", "[NH4+]"]);
        let t = tokenize("C[O-]");
        assert_eq!(t[1], Token { text: "[O-]", kind: Kind::Atom, pos: 1 });
    }

    #[test]
    fn tokens_outlive_lexer_and_are_zero_copy() {
        let src = String::from("c1ccccc1");
        let toks = {
            let mut lx = Lexer::new(&src);
            let a = lx.next_token().unwrap();
            let b = lx.next_token().unwrap(); // `a` still alive here
            vec![a, b]
        }; // lexer dropped
        assert_eq!(toks[1].kind, Kind::Ring);
        assert_eq!(toks[0].text.as_ptr(), src.as_ptr());
    }

    #[test]
    fn longer() {
        let src = "Cl[Fe+2]";
        let t = tokenize(src);
        let s;
        {
            let (a, b) = (t[0], t[1]);
            s = longer_text(&a, &b);
        }
        assert_eq!(s, "[Fe+2]");
    }

    #[test]
    fn uppercase_symbols() {
        let t = tokenize("c1ncc[nH]1");
        assert_eq!(upper_atoms(&t), vec!["C", "N", "C", "C", "[NH]"]);
    }

    #[test]
    fn closure_filter() {
        let t = tokenize("CCl.CBr.Cl");
        let prefix = String::from("Cl");
        let f = starts_with_filter(&prefix);
        assert_eq!(t.iter().filter(|x| f(x)).count(), 2);
    }

    #[test]
    fn symbol_table() {
        let mut s = SymbolTable::new();
        assert_eq!(s.push_and_first_len("Br"), 1);
        assert_eq!(s.len(), 2);
    }
}
