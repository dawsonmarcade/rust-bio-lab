//! # Exercise 09 — Zero-copy PDB parsing
//!
//! **Concepts:** structs with lifetime parameters, *two* lifetimes in one
//! type (`Residues<'s, 'a>`), iterators that yield borrowed slices,
//! returning `impl Iterator + '_`, fixed-width byte parsing.
//!
//! A PDB file can be hundreds of MB. Instead of allocating a `String` for
//! every atom name, `Atom<'a>` *borrows* its text fields from the input
//! buffer. The `'a` lifetime says "this atom cannot outlive the text it was
//! parsed from".
//!
//! PDB ATOM/HETATM columns (**1-based, inclusive**):
//!
//! | field      | cols  |   | field       | cols  |
//! |------------|-------|---|-------------|-------|
//! | record     | 1-6   |   | x           | 31-38 |
//! | serial     | 7-11  |   | y           | 39-46 |
//! | atom name  | 13-16 |   | z           | 47-54 |
//! | residue    | 18-20 |   | occupancy   | 55-60 |
//! | chain      | 22    |   | B-factor    | 61-66 |
//! | res. seq   | 23-26 |   | element     | 77-78 |
//!
//! Parsing rules:
//! * Only lines starting with `ATOM  ` or `HETATM` are atoms; ignore others.
//! * A line shorter than 54 characters is `Truncated`.
//! * Text fields are trimmed. Occupancy/B-factor may be missing (short line)
//!   or blank -> `0.0`.
//! * Element: cols 77-78 trimmed; if missing or blank, use the first ASCII
//!   letter of the atom name.
//! * A field that fails to parse -> `BadField { line, field }` where `field`
//!   is one of `"serial" "res_seq" "x" "y" "z" "b_factor"`. Lines are 1-based.
//!
//! ```text
//! cargo test -p ex09_pdb_parser
//! ```

#![allow(unused_variables, dead_code, unused_imports, unused_mut)]

#[derive(Debug, Clone, PartialEq)]
pub struct Atom<'a> {
    pub serial: u32,
    pub name: &'a str,
    pub res_name: &'a str,
    pub chain: char,
    pub res_seq: i32,
    pub pos: [f64; 3],
    pub b_factor: f64,
    pub element: &'a str,
    pub hetero: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum PdbError {
    Truncated { line: usize },
    BadField { line: usize, field: &'static str },
}

#[derive(Debug, Clone, PartialEq)]
pub struct Structure<'a> {
    pub atoms: Vec<Atom<'a>>,
}

pub fn parse_pdb(text: &str) -> Result<Structure<'_>, PdbError> {
    todo!()
}

/// Iterator over runs of consecutive atoms sharing `(chain, res_seq, res_name)`.
/// `'s` is how long we borrow the `Structure`; `'a` is the lifetime of the
/// text inside the atoms. Yields slices *into* `Structure::atoms`.
pub struct Residues<'s, 'a> {
    rest: &'s [Atom<'a>],
    // anything else you need
}

impl<'s, 'a> Iterator for Residues<'s, 'a> {
    type Item = &'s [Atom<'a>];
    fn next(&mut self) -> Option<Self::Item> {
        todo!()
    }
}

fn atomic_mass(element: &str) -> f64 {
    match element.to_ascii_uppercase().as_str() {
        "H" => 1.008,
        "C" => 12.011,
        "N" => 14.007,
        "O" => 15.999,
        "P" => 30.974,
        "S" => 32.06,
        _ => 0.0,
    }
}

fn three_to_one(res: &str) -> char {
    match res {
        "ALA" => 'A', "ARG" => 'R', "ASN" => 'N', "ASP" => 'D', "CYS" => 'C',
        "GLN" => 'Q', "GLU" => 'E', "GLY" => 'G', "HIS" => 'H', "ILE" => 'I',
        "LEU" => 'L', "LYS" => 'K', "MET" => 'M', "PHE" => 'F', "PRO" => 'P',
        "SER" => 'S', "THR" => 'T', "TRP" => 'W', "TYR" => 'Y', "VAL" => 'V',
        _ => 'X',
    }
}

impl<'a> Structure<'a> {
    pub fn residues(&self) -> Residues<'_, 'a> {
        todo!()
    }

    /// Alpha carbons: non-hetero atoms named `CA`.
    pub fn ca_atoms(&self) -> impl Iterator<Item = &Atom<'a>> + '_ {
        if true {
            todo!()
        }
        std::iter::empty() // placeholder so the signature type-checks
    }

    /// Mass-weighted centre using the `atomic_mass` table (unknown
    /// elements weigh 0). `None` if total mass is 0.
    pub fn center_of_mass(&self) -> Option<[f64; 3]> {
        todo!()
    }

    /// Mass-weighted radius of gyration: sqrt( Σ m_i |r_i - com|² / Σ m_i ).
    pub fn radius_of_gyration(&self) -> Option<f64> {
        todo!()
    }

    /// One-letter sequence per chain, in order of first appearance, built from
    /// residues that contain a `CA` atom. Unknown residues -> `'X'`.
    /// Returns `(chain, sequence)` pairs.
    pub fn sequences(&self) -> Vec<(char, String)> {
        todo!()
    }

    /// Pairs of CA atoms `(i, j)` with `i` before `j` in the file, distance
    /// `<= cutoff`, and **not** sequence-local: skip pairs on the same chain
    /// with `|res_seq_i - res_seq_j| < 3`. Returned as references, so the
    /// caller can read any field. O(n²) is fine.
    pub fn ca_contacts(&self, cutoff: f64) -> Vec<(&Atom<'a>, &Atom<'a>)> {
        todo!()
    }
}

/// Lifetime puzzle: return the name of the atom with the highest B-factor
/// across *two different structures* parsed from two different texts.
/// Figure out what lifetime annotations this signature needs so that the
/// result can outlive both `Structure` values but not the texts. (It won't
/// compile until you fix the signature!)
pub fn hottest_atom_name(a: &Structure, b: &Structure) -> Option<&str> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[allow(clippy::too_many_arguments)]
    fn atom_line(rec: &str, serial: u32, name: &str, res: &str, chain: char, seq: i32, xyz: [f64; 3], b: f64, el: &str) -> String {
        let name_field = if name.len() < 4 { format!(" {name:<3}") } else { name.to_string() };
        format!(
            "{rec:<6}{serial:>5} {name_field:<4} {res:>3} {chain}{seq:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{b:>6.2}          {el:>2}",
            xyz[0], xyz[1], xyz[2], 1.0
        )
    }

    fn sample() -> String {
        let mut lines = vec!["HEADER    TEST".to_string(), "REMARK  nothing".to_string()];
        let mut serial = 1;
        // chain A: 5 residues along x, 3.8 Å apart; chain B: 2 residues near residue A1
        let res = ["MET", "ALA", "GLY", "LYS", "TRP"];
        for (i, r) in res.iter().enumerate() {
            let x = i as f64 * 3.8;
            lines.push(atom_line("ATOM", serial, "N", r, 'A', i as i32 + 1, [x - 1.0, 0.0, 0.0], 10.0, "N"));
            serial += 1;
            lines.push(atom_line("ATOM", serial, "CA", r, 'A', i as i32 + 1, [x, 0.0, 0.0], 20.0 + i as f64, "C"));
            serial += 1;
        }
        lines.push(atom_line("ATOM", serial, "CA", "SER", 'B', 1, [0.0, 5.0, 0.0], 15.0, "C"));
        serial += 1;
        lines.push(atom_line("ATOM", serial, "CA", "XYZ", 'B', 2, [15.2, 4.0, 0.0], 15.0, "C"));
        serial += 1;
        lines.push(atom_line("HETATM", serial, "O", "HOH", 'W', 1, [100.0, 100.0, 100.0], 99.0, ""));
        lines.push("TER".to_string());
        lines.push("END".to_string());
        lines.join("\n")
    }

    #[test]
    fn parses_fields() {
        let text = sample();
        let s = parse_pdb(&text).unwrap();
        assert_eq!(s.atoms.len(), 13);
        let a = &s.atoms[1];
        assert_eq!((a.serial, a.name, a.res_name, a.chain, a.res_seq), (2, "CA", "MET", 'A', 1));
        assert_eq!(a.pos, [0.0, 0.0, 0.0]);
        assert_eq!(a.b_factor, 20.0);
        assert_eq!(a.element, "C");
        let w = s.atoms.last().unwrap();
        assert!(w.hetero);
        assert_eq!(w.element, "O"); // inferred from name
        // zero-copy: names point into the original text
        let base = text.as_ptr() as usize;
        let p = a.name.as_ptr() as usize;
        assert!(p >= base && p < base + text.len());
    }

    #[test]
    fn short_lines_and_errors() {
        let ok = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000";
        let s = parse_pdb(ok).unwrap();
        assert_eq!(s.atoms[0].b_factor, 0.0);
        assert_eq!(s.atoms[0].element, "C");
        assert_eq!(parse_pdb("HEADER\nATOM      1  CA  ALA A   1       1.000   2.000"), Err(PdbError::Truncated { line: 2 }));
        let bad = "ATOM      1  CA  ALA A   1       1.000   2.0x0   3.000";
        assert_eq!(parse_pdb(bad), Err(PdbError::BadField { line: 1, field: "y" }));
        let bad = "ATOM      1  CA  ALA A   ?       1.000   2.000   3.000";
        assert_eq!(parse_pdb(bad), Err(PdbError::BadField { line: 1, field: "res_seq" }));
    }

    #[test]
    fn residues_iterator() {
        let text = sample();
        let s = parse_pdb(&text).unwrap();
        let res: Vec<&[Atom]> = s.residues().collect();
        assert_eq!(res.len(), 8);
        assert_eq!(res[0].len(), 2);
        assert_eq!(res[4][1].res_name, "TRP");
        assert_eq!(res[7][0].res_name, "HOH");
        // slices point into s.atoms
        assert!(std::ptr::eq(res[1].as_ptr(), &s.atoms[2]));
    }

    #[test]
    fn sequences_by_chain() {
        let text = sample();
        let s = parse_pdb(&text).unwrap();
        assert_eq!(s.sequences(), vec![('A', "MAGKW".to_string()), ('B', "SX".to_string())]);
    }

    #[test]
    fn geometry() {
        let text = [
            "ATOM      1  C   XXX A   1       0.000   0.000   0.000",
            "ATOM      2  C   XXX A   1       2.000   0.000   0.000",
            "ATOM      3 ZN   XXX A   1      50.000   0.000   0.000  1.00  0.00          ZN",
        ]
        .join("\n");
        let s = parse_pdb(&text).unwrap();
        assert_eq!(s.center_of_mass(), Some([1.0, 0.0, 0.0]));
        assert!((s.radius_of_gyration().unwrap() - 1.0).abs() < 1e-12);
        let empty = parse_pdb("HEADER only").unwrap();
        assert_eq!(empty.center_of_mass(), None);
        assert_eq!(empty.radius_of_gyration(), None);
    }

    #[test]
    fn contacts() {
        let text = sample();
        let s = parse_pdb(&text).unwrap();
        let c: Vec<(char, i32, char, i32)> = s
            .ca_contacts(6.0)
            .iter()
            .map(|(a, b)| (a.chain, a.res_seq, b.chain, b.res_seq))
            .collect();
        // A1-A2, A1-A3 etc. are sequence-local; B1 is 5 Å from A1 but 6.28 Å
        // from A2; B2 is 4 Å from A5 and 5.52 Å from A4.
        assert_eq!(c, vec![('A', 1, 'B', 1), ('A', 4, 'B', 2), ('A', 5, 'B', 2)]);
        // 7 CAs -> 21 pairs, minus 7 local pairs in chain A and 1 in chain B
        assert_eq!(s.ca_contacts(100.0).len(), 21 - 8);
    }

    #[test]
    fn lifetimes_across_structures() {
        let t1 = sample();
        let t2 = "ATOM      1  OXT GLY A   1       0.000   0.000   0.000  1.00 99.50           O".to_string();
        let name;
        {
            let s1 = parse_pdb(&t1).unwrap();
            let s2 = parse_pdb(&t2).unwrap();
            name = hottest_atom_name(&s1, &s2);
        } // s1 and s2 dropped here, but `name` borrows the *texts*
        assert_eq!(name, Some("OXT")); // B = 99.50 beats the water's 99.00
    }
}
