//! # Exercise 15 — Parse SMILES into a molecular graph
//!
//! **Concepts:** writing a real parser by hand (a state machine over
//! bytes), graphs as index-based edge lists (no `Rc<RefCell>` needed),
//! BFS, union-find or DFS for connected components, `HashMap`/`BTreeMap`,
//! rich error enums.
//!
//! SMILES is how chemists write molecules as text: aspirin is
//! `CC(=O)Oc1ccccc1C(=O)O`. Support this subset:
//!
//! * **Organic-subset atoms** (outside brackets): `B C N O P S F Cl Br I`
//!   and aromatic `b c n o p s`. Any other letter outside brackets is
//!   `UnexpectedChar`.
//! * **Bracket atoms**: `[` optional isotope digits (ignored), element symbol
//!   (an uppercase letter plus an optional lowercase letter, e.g. `Fe`, `Na`,
//!   or a single aromatic lowercase letter), optional `H` with optional count,
//!   optional charge (`+`, `-`, `+2`, `-3`, `++`, `--`), then `]`. A bracket
//!   atom has exactly the hydrogens written: `[CH4]` has 4 and `[C]` has 0.
//!   A missing `]` is `UnclosedBracket { pos }` where `pos` is the `[`.
//! * **Bonds** `-` `=` `#` `:`. With no symbol, the bond between two
//!   aromatic atoms is `Aromatic`, otherwise `Single`.
//! * **Branches** `(...)`, nestable.
//! * **Ring closures**: a digit `1`-`9`, or `%` + two digits. A bond symbol
//!   may appear at the opening, at the closing, or both (must agree, else
//!   `RingBondMismatch`).
//! * `.` separates disconnected components.
//!
//! **Implicit hydrogens** (organic subset only): let `s` = sum of bond orders
//! (aromatic bond counts 1). Default valences: B 3; C 4; N 3, 5; O 2;
//! P 3, 5; S 2, 4, 6; F/Cl/Br/I 1. Let `v` = the smallest default valence
//! `>= s` (if none, H = 0). Then H = `v - s`, minus 1 more if the atom is
//! aromatic `b`, `c`, `n` or `p` (its π bond), never going below 0.
//! (Check: benzene `c` has s = 2, v = 4, H = 4 - 2 - 1 = 1. Pyridine `n` has
//! s = 2, v = 3, H = 0.)
//!
//! **Formula** in Hill order: C first, then H, then everything else
//! alphabetically. If there is no carbon, *everything* (H included) is
//! alphabetical. Counts of 1 are omitted.
//!
//! Positions in errors are 0-based byte offsets.
//!
//! ```text
//! cargo test -p ex15_smiles_parser
//! ```

#![allow(unused_variables, dead_code, unused_imports, unused_mut)]

use std::collections::{BTreeMap, HashMap, VecDeque};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
}

impl BondOrder {
    /// Contribution to valence (aromatic counts as 1).
    pub fn valence(self) -> u32 {
        match self {
            BondOrder::Single | BondOrder::Aromatic => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Atom {
    /// Capitalised symbol: `"C"`, `"Cl"`, `"N"` (also for aromatic `n`).
    pub element: String,
    pub aromatic: bool,
    pub charge: i8,
    /// `Some(h)` for bracket atoms, `None` for organic-subset atoms.
    pub bracket_h: Option<u8>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Bond {
    pub a: usize,
    pub b: usize,
    pub order: BondOrder,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    /// In the order they were created while reading left to right.
    pub bonds: Vec<Bond>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SmilesError {
    Empty,
    UnexpectedChar { pos: usize, ch: char },
    UnclosedBracket { pos: usize },
    UnmatchedParen { pos: usize },
    UnclosedBranch,
    UnclosedRing { label: u32 },
    RingBondMismatch { label: u32 },
    /// A bond symbol with no atom after it (end of input, `)` or `.`).
    DanglingBond { pos: usize },
}

pub fn parse_smiles(s: &str) -> Result<Molecule, SmilesError> {
    todo!()
}

fn atomic_mass(e: &str) -> Option<f64> {
    Some(match e {
        "H" => 1.008,
        "B" => 10.81,
        "C" => 12.011,
        "N" => 14.007,
        "O" => 15.999,
        "F" => 18.998,
        "Na" => 22.990,
        "P" => 30.974,
        "S" => 32.06,
        "Cl" => 35.45,
        "K" => 39.098,
        "Fe" => 55.845,
        "Br" => 79.904,
        "I" => 126.904,
        _ => return None,
    })
}

impl Molecule {
    /// `(neighbour index, bond order)` for every bond touching atom `i`,
    /// in bond-creation order.
    pub fn neighbors(&self, i: usize) -> impl Iterator<Item = (usize, BondOrder)> + '_ {
        if true {
            todo!()
        }
        std::iter::empty() // placeholder so the signature type-checks
    }

    pub fn hydrogens(&self, i: usize) -> u32 {
        todo!()
    }

    pub fn formula(&self) -> String {
        todo!()
    }

    /// Average molecular weight including implicit H. `None` if any element
    /// is not in the mass table.
    pub fn molecular_weight(&self) -> Option<f64> {
        todo!()
    }

    pub fn total_charge(&self) -> i32 {
        todo!()
    }

    pub fn components(&self) -> usize {
        todo!()
    }

    /// Number of independent rings (cyclomatic number): bonds − atoms + components.
    pub fn ring_count(&self) -> usize {
        todo!()
    }

    /// Size of the smallest ring containing atom `i`, or `None`.
    pub fn smallest_ring(&self, i: usize) -> Option<usize> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn formula(s: &str) -> String {
        parse_smiles(s).unwrap_or_else(|e| panic!("{s}: {e:?}")).formula()
    }

    #[test]
    fn simple_aliphatics() {
        let m = parse_smiles("CCO").unwrap();
        assert_eq!(m.atoms.len(), 3);
        assert_eq!(m.bonds, vec![
            Bond { a: 0, b: 1, order: BondOrder::Single },
            Bond { a: 1, b: 2, order: BondOrder::Single }
        ]);
        assert_eq!(m.formula(), "C2H6O");
        assert!((m.molecular_weight().unwrap() - 46.069).abs() < 1e-9);
        assert_eq!(formula("C#C"), "C2H2");
        assert_eq!(formula("C#N"), "CHN");
        assert_eq!(formula("ClCCl"), "CH2Cl2");
        assert_eq!(formula("CC(C)(C)C"), "C5H12");
        assert_eq!(formula("O=C=O"), "CO2");
        assert_eq!(formula("CS(=O)(=O)C"), "C2H6O2S");
        assert_eq!(formula("CP(=O)(O)O"), "CH5O3P");
    }

    #[test]
    fn branches_and_neighbors() {
        let m = parse_smiles("CC(=O)O").unwrap();
        let n: Vec<_> = m.neighbors(1).collect();
        assert_eq!(n, vec![(0, BondOrder::Single), (2, BondOrder::Double), (3, BondOrder::Single)]);
        assert_eq!(m.formula(), "C2H4O2");
        assert_eq!(formula("CC(=O)Oc1ccccc1C(=O)O"), "C9H8O4"); // aspirin
    }

    #[test]
    fn aromatics() {
        assert_eq!(formula("c1ccccc1"), "C6H6");
        assert_eq!(formula("Cc1ccccc1"), "C7H8");
        assert_eq!(formula("c1ccncc1"), "C5H5N");
        assert_eq!(formula("c1ccoc1"), "C4H4O");
        assert_eq!(formula("c1ccsc1"), "C4H4S");
        assert_eq!(formula("c1cc[nH]c1"), "C4H5N");
        assert_eq!(formula("c1ccc2ccccc2c1"), "C10H8");
        let m = parse_smiles("c1ccccc1").unwrap();
        assert!(m.bonds.iter().all(|b| b.order == BondOrder::Aromatic));
        // explicit single bond between aromatic atoms (biphenyl)
        let m = parse_smiles("c1ccccc1-c1ccccc1").unwrap();
        assert_eq!(m.formula(), "C12H10");
        assert_eq!(m.bonds.iter().filter(|b| b.order == BondOrder::Single).count(), 1);
    }

    #[test]
    fn caffeine() {
        let m = parse_smiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C").unwrap();
        assert_eq!(m.formula(), "C8H10N4O2");
        assert!((m.molecular_weight().unwrap() - 194.194).abs() < 0.01);
        assert_eq!(m.ring_count(), 2);
        assert_eq!(m.components(), 1);
    }

    #[test]
    fn rings() {
        let m = parse_smiles("C1CCCCC1").unwrap();
        assert_eq!((m.formula().as_str(), m.ring_count()), ("C6H12", 1));
        assert_eq!(m.smallest_ring(3), Some(6));
        assert_eq!(formula("C=1CCCCC1"), "C6H10");
        assert_eq!(formula("C1CCCCC=1"), "C6H10");
        assert_eq!(formula("C%10CCC%10"), "C4H8");
        let cubane = parse_smiles("C12C3C4C1C5C2C3C45").unwrap();
        assert_eq!(cubane.formula(), "C8H8");
        assert_eq!(cubane.ring_count(), 5);
        assert!((0..8).all(|i| cubane.smallest_ring(i) == Some(4)));
        let naph = parse_smiles("c1ccc2ccccc2c1").unwrap();
        assert_eq!(naph.ring_count(), 2);
        assert_eq!(naph.smallest_ring(4), Some(6));
        let m = parse_smiles("C1CC1CCC").unwrap();
        assert_eq!(m.smallest_ring(4), None);
        assert_eq!(m.smallest_ring(1), Some(3));
        // ring label reuse after closing
        assert_eq!(parse_smiles("C1CC1C1CC1").unwrap().ring_count(), 2);
    }

    #[test]
    fn brackets_and_charges() {
        let m = parse_smiles("[NH4+].[Cl-]").unwrap();
        assert_eq!(m.formula(), "ClH4N");
        assert_eq!(m.components(), 2);
        assert_eq!(m.total_charge(), 0);
        let m = parse_smiles("CC(=O)[O-]").unwrap();
        assert_eq!(m.formula(), "C2H3O2");
        assert_eq!(m.total_charge(), -1);
        assert_eq!(m.atoms[3], Atom { element: "O".into(), aromatic: false, charge: -1, bracket_h: Some(0) });
        let fe = parse_smiles("[Fe+2]").unwrap();
        assert_eq!((fe.atoms[0].element.as_str(), fe.atoms[0].charge), ("Fe", 2));
        assert_eq!(parse_smiles("[Fe++]").unwrap().atoms[0].charge, 2);
        assert_eq!(parse_smiles("[13CH4]").unwrap().formula(), "CH4");
        assert_eq!(parse_smiles("[C]").unwrap().formula(), "C");
        assert_eq!(parse_smiles("[Na+].[Cl-]").unwrap().formula(), "ClNa");
        assert_eq!(parse_smiles("[Xx]").unwrap().molecular_weight(), None);
    }

    #[test]
    fn errors() {
        use SmilesError::*;
        assert_eq!(parse_smiles(""), Err(Empty));
        assert_eq!(parse_smiles("C(C"), Err(UnclosedBranch));
        assert_eq!(parse_smiles("CC)C"), Err(UnmatchedParen { pos: 2 }));
        assert_eq!(parse_smiles("C1CC"), Err(UnclosedRing { label: 1 }));
        assert_eq!(parse_smiles("C1CC2"), Err(UnclosedRing { label: 1 }));
        assert_eq!(parse_smiles("CX"), Err(UnexpectedChar { pos: 1, ch: 'X' }));
        assert_eq!(parse_smiles("C="), Err(DanglingBond { pos: 1 }));
        assert_eq!(parse_smiles("CC(=)C"), Err(DanglingBond { pos: 3 }));
        assert_eq!(parse_smiles("C=1CCCCC#1"), Err(RingBondMismatch { label: 1 }));
        assert_eq!(parse_smiles("C[NH4"), Err(UnclosedBracket { pos: 1 }));
        assert_eq!(parse_smiles("=C"), Err(UnexpectedChar { pos: 0, ch: '=' }));
        assert_eq!(parse_smiles("(C)"), Err(UnexpectedChar { pos: 0, ch: '(' }));
    }
}
