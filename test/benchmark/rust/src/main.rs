//! Brute-force off-target ground-truth generator (Rust port).
//!
//! Exhaustive, budget-bounded Needleman-Wunsch aligner. For a given guide
//! (spacer+PAM, IUPAC allowed) it scans a per-chromosome FASTA on both strands
//! and reports every alignment within the mismatch / DNA-gap / RNA-gap budget,
//! applying PAM + edge-gap filters.
//!
//! This is a drop-in, faster replacement for `generate_brute_force.py`,
//! producing byte-for-byte identical output (as a set of rows). std-only,
//! no external crate dependencies.

use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Read, Write};

// ---------------------------------------------------------------------------
// IUPAC matching
// ---------------------------------------------------------------------------
// Each IUPAC code is represented as a bitmask over the 5 concrete symbols
// {A, C, G, T, U}. iupac_match == (mask(a) & mask(b) != 0).
//
// Mirrors the Python iupac_dict exactly (T and U both include {T, U}).
const A: u8 = 1 << 0;
const C: u8 = 1 << 1;
const G: u8 = 1 << 2;
const T: u8 = 1 << 3;
const U: u8 = 1 << 4;

#[inline]
fn iupac_mask(b: u8) -> u8 {
    // Case-insensitive: uppercase the ASCII byte.
    let up = b.to_ascii_uppercase();
    match up {
        b'A' => A,
        b'C' => C,
        b'G' => G,
        b'T' => T | U,
        b'U' => T | U,
        b'R' => A | G,
        b'Y' => C | T,
        b'S' => G | C,
        b'W' => A | T,
        b'K' => G | T,
        b'M' => A | C,
        b'B' => C | G | T,
        b'D' => A | G | T,
        b'H' => A | C | T,
        b'V' => A | C | G,
        b'N' => A | C | G | T,
        _ => 0, // unknown code -> empty set, matches iupac_dict.get(..., set())
    }
}

#[inline]
fn iupac_match(a: u8, b: u8) -> bool {
    (iupac_mask(a) & iupac_mask(b)) != 0
}

// ---------------------------------------------------------------------------
// Alignment result
// ---------------------------------------------------------------------------
#[derive(Clone)]
struct Alignment {
    aln1: Vec<u8>, // RNA aligned string
    aln2: Vec<u8>, // DNA aligned string
    mismatches: usize,
    rna_gaps: usize,
    dna_gaps: usize,
}

// A count triple (mismatch, dna_gap, rna_gap).
type Cnt = (i32, i32, i32);
// A move: (kind, pi, pj). kind: 0=diag, 1=up, 2=left.
type Move = (u8, usize, usize);

// ---------------------------------------------------------------------------
// needleman_wunsch_all_alignments
// ---------------------------------------------------------------------------
fn needleman_wunsch_all_alignments(
    rna: &[u8],
    dna: &[u8],
    max_mismatches: i32,
    max_dna_gaps: i32,
    max_rna_gaps: i32,
) -> Vec<Alignment> {
    let m = rna.len();
    let n = dna.len();

    // dp_score[i][j]: minimal cost; f64::INFINITY sentinel.
    let mut dp_score = vec![vec![f64::INFINITY; n + 1]; m + 1];
    // dp_cells[i][j]: map cnt -> list of moves.
    let mut dp_cells: Vec<Vec<HashMap<Cnt, Vec<Move>>>> =
        vec![vec![HashMap::new(); n + 1]; m + 1];

    dp_score[0][0] = 0.0;
    dp_cells[0][0].insert((0, 0, 0), Vec::new());

    for i in 1..=m {
        if (i as i32) <= max_dna_gaps {
            dp_score[i][0] = i as f64;
            dp_cells[i][0].insert((0, i as i32, 0), vec![(1u8, i - 1, 0)]);
        }
    }
    for j in 1..=n {
        if (j as i32) <= max_rna_gaps {
            dp_score[0][j] = j as f64;
            dp_cells[0][j].insert((0, 0, j as i32), vec![(2u8, 0, j - 1)]);
        }
    }

    for i in 1..=m {
        for j in 1..=n {
            // candidates: (cost, cnt, move_kind, pi, pj)
            let mut candidates: Vec<(f64, Cnt, u8, usize, usize)> = Vec::new();

            // diag
            let mis: i32 = if iupac_match(rna[i - 1], dna[j - 1]) { 0 } else { 1 };
            let base_diag = dp_score[i - 1][j - 1];
            for &(mm, dg, rg) in dp_cells[i - 1][j - 1].keys() {
                let new = (mm + mis, dg, rg);
                if new.0 <= max_mismatches {
                    candidates.push((base_diag + mis as f64, new, 0u8, i - 1, j - 1));
                }
            }

            // up (DNA gap): consume RNA base, gap in DNA
            let base_up = dp_score[i - 1][j];
            for &(mm, dg, rg) in dp_cells[i - 1][j].keys() {
                let new = (mm, dg + 1, rg);
                if new.1 <= max_dna_gaps {
                    candidates.push((base_up + 1.0, new, 1u8, i - 1, j));
                }
            }

            // left (RNA gap): consume DNA base, gap in RNA
            let base_left = dp_score[i][j - 1];
            for &(mm, dg, rg) in dp_cells[i][j - 1].keys() {
                let (new, cost) = if i == m {
                    ((mm, dg, rg), base_left)
                } else {
                    ((mm, dg, rg + 1), base_left + 1.0)
                };
                if new.2 <= max_rna_gaps {
                    candidates.push((cost, new, 2u8, i, j - 1));
                }
            }

            if candidates.is_empty() {
                continue;
            }

            let mut best = f64::INFINITY;
            for c in &candidates {
                if c.0 < best {
                    best = c.0;
                }
            }
            dp_score[i][j] = best;

            let cell = &mut dp_cells[i][j];
            for (_cost, cnt, mv, pi, pj) in candidates {
                cell.entry(cnt).or_insert_with(Vec::new).push((mv, pi, pj));
            }
        }
    }

    let mut alignments: Vec<Alignment> = Vec::new();
    let mut seen: HashSet<(Vec<u8>, Vec<u8>, usize, usize, usize)> = HashSet::new();

    // Collect the final cnts before borrowing dp_cells in recursion.
    let final_cnts: Vec<Cnt> = dp_cells[m][n].keys().copied().collect();

    // Iterative backtrack via explicit stack to avoid recursion overhead / stack
    // depth issues. State: (i, j, cnt, a1, a2) where a1/a2 are built in reverse.
    struct Frame {
        i: usize,
        j: usize,
        cnt: Cnt,
        a1: Vec<u8>,
        a2: Vec<u8>,
    }

    for final_cnt in final_cnts {
        let mut stack: Vec<Frame> = vec![Frame {
            i: m,
            j: n,
            cnt: final_cnt,
            a1: Vec::new(),
            a2: Vec::new(),
        }];

        while let Some(fr) = stack.pop() {
            if fr.i == 0 && fr.j == 0 {
                // Reverse accumulated strings.
                let mut aln1 = fr.a1.clone();
                aln1.reverse();
                let mut aln2 = fr.a2.clone();
                aln2.reverse();

                // mismatches: over aligned non-gap columns.
                let mut mismatches = 0usize;
                for k in 0..aln1.len() {
                    let x = aln1[k];
                    let y = aln2[k];
                    if x != b'-' && y != b'-' && !iupac_match(x, y) {
                        mismatches += 1;
                    }
                }
                // dna_gaps = aln2.count('-')
                let dna_gaps = aln2.iter().filter(|&&c| c == b'-').count();
                // rna_gaps = aln1.rstrip('-').count('-')
                let rstrip_len = {
                    let mut e = aln1.len();
                    while e > 0 && aln1[e - 1] == b'-' {
                        e -= 1;
                    }
                    e
                };
                let rna_gaps = aln1[..rstrip_len].iter().filter(|&&c| c == b'-').count();

                let key = (
                    aln1.clone(),
                    aln2.clone(),
                    mismatches,
                    dna_gaps,
                    rna_gaps,
                );
                if seen.insert(key) {
                    alignments.push(Alignment {
                        aln1,
                        aln2,
                        mismatches,
                        rna_gaps,
                        dna_gaps,
                    });
                }
                continue;
            }

            let moves = match dp_cells[fr.i][fr.j].get(&fr.cnt) {
                Some(v) => v,
                None => continue,
            };
            for &(kind, pi, pj) in moves {
                match kind {
                    0 => {
                        // diag
                        let mis: i32 =
                            if iupac_match(rna[fr.i - 1], dna[fr.j - 1]) { 0 } else { 1 };
                        let mut a1 = fr.a1.clone();
                        a1.push(rna[fr.i - 1]);
                        let mut a2 = fr.a2.clone();
                        a2.push(dna[fr.j - 1]);
                        stack.push(Frame {
                            i: pi,
                            j: pj,
                            cnt: (fr.cnt.0 - mis, fr.cnt.1, fr.cnt.2),
                            a1,
                            a2,
                        });
                    }
                    1 => {
                        // up
                        let mut a1 = fr.a1.clone();
                        a1.push(rna[fr.i - 1]);
                        let mut a2 = fr.a2.clone();
                        a2.push(b'-');
                        stack.push(Frame {
                            i: pi,
                            j: pj,
                            cnt: (fr.cnt.0, fr.cnt.1 - 1, fr.cnt.2),
                            a1,
                            a2,
                        });
                    }
                    _ => {
                        // left
                        let prev = if fr.i == m {
                            fr.cnt
                        } else {
                            (fr.cnt.0, fr.cnt.1, fr.cnt.2 - 1)
                        };
                        let mut a1 = fr.a1.clone();
                        a1.push(b'-');
                        let mut a2 = fr.a2.clone();
                        a2.push(dna[fr.j - 1]);
                        stack.push(Frame {
                            i: pi,
                            j: pj,
                            cnt: prev,
                            a1,
                            a2,
                        });
                    }
                }
            }
        }
    }

    alignments
}

// ---------------------------------------------------------------------------
// reverse_complement
// ---------------------------------------------------------------------------
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        let up = b.to_ascii_uppercase();
        let c = match up {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            b'N' => b'N',
            other => other, // complement.get(b, b): unknown maps to itself (uppercased)
        };
        out.push(c);
    }
    out
}

// ---------------------------------------------------------------------------
// keep_alignment
// ---------------------------------------------------------------------------
fn keep_alignment(rna_out: &[u8], dna_out: &[u8], pam: &[u8], pam_5prime: bool) -> bool {
    if rna_out.is_empty() || rna_out[0] == b'-' || rna_out[rna_out.len() - 1] == b'-' {
        return false;
    }
    if !pam.is_empty() {
        // d = dna_out.upper(); seg = first/last len(pam) chars.
        let plen = pam.len();
        let dlen = dna_out.len();
        let seg: &[u8] = if pam_5prime {
            if dlen < plen {
                return false;
            }
            &dna_out[..plen]
        } else {
            if dlen < plen {
                return false;
            }
            &dna_out[dlen - plen..]
        };
        if seg.len() < plen {
            return false;
        }
        for &s in seg {
            if s == b'-' {
                return false;
            }
        }
        for k in 0..plen {
            // pam is already uppercase; seg is uppercased via iupac_match.
            if !iupac_match(pam[k], seg[k]) {
                return false;
            }
        }
    }
    true
}

// ---------------------------------------------------------------------------
// scan
// ---------------------------------------------------------------------------
#[allow(clippy::too_many_arguments)]
fn scan<W: Write>(
    rna: &[u8],
    dna: &[u8],
    strand: &str,
    max_dna_gaps: i32,
    max_rna_gaps: i32,
    max_mismatches: i32,
    chrom: &str,
    writer: &mut W,
    pam: &[u8],
    pam_5prime: bool,
) {
    let rna_len = rna.len();
    let dna_len = dna.len();
    // range(0, len(DNA) - len(RNA) + max_dna_gaps + 1)
    // Python: if the upper bound is <= 0, the loop simply doesn't execute.
    let upper = dna_len as i64 - rna_len as i64 + max_dna_gaps as i64 + 1;
    if upper <= 0 {
        return;
    }
    let upper = upper as usize;

    for starting_pos in 0..upper {
        let slice_end = std::cmp::min(starting_pos + rna_len + max_rna_gaps as usize, dna_len);
        let dna_slice = &dna[starting_pos..slice_end];

        // Skip windows containing 'N' or 'n'.
        let mut has_n = false;
        for &b in dna_slice {
            if b == b'N' || b == b'n' {
                has_n = true;
                break;
            }
        }
        if has_n {
            continue;
        }

        let aligns = needleman_wunsch_all_alignments(
            rna,
            dna_slice,
            max_mismatches,
            max_dna_gaps,
            max_rna_gaps,
        );
        if aligns.is_empty() {
            continue;
        }

        for al in &aligns {
            // actual_len = len(r_aln.rstrip('-'))
            let actual_len = {
                let mut e = al.aln1.len();
                while e > 0 && al.aln1[e - 1] == b'-' {
                    e -= 1;
                }
                e
            };
            let r_out = &al.aln1[..actual_len];
            let d_out = &al.aln2[..actual_len];

            if !keep_alignment(r_out, d_out, pam, pam_5prime) {
                continue;
            }

            let (start, end): (i64, i64) = if strand == "+" {
                (starting_pos as i64, (starting_pos + actual_len) as i64)
            } else {
                let dg = al.dna_gaps as i64;
                (
                    dna_len as i64 - starting_pos as i64 - actual_len as i64 + dg,
                    dna_len as i64 - starting_pos as i64 + dg,
                )
            };

            // Write TSV row.
            writer
                .write_all(chrom.as_bytes())
                .and_then(|_| writer.write_all(b"\t"))
                .and_then(|_| writer.write_all(r_out))
                .and_then(|_| writer.write_all(b"\t"))
                .and_then(|_| writer.write_all(d_out))
                .and_then(|_| write!(
                    writer,
                    "\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    strand, start, end, al.mismatches, al.rna_gaps, al.dna_gaps
                ))
                .expect("write failed");
        }
    }
}

// ---------------------------------------------------------------------------
// FASTA reader (single record): skip '>' header, concatenate seq lines, keep case.
// ---------------------------------------------------------------------------
fn read_single_fasta(path: &str) -> std::io::Result<Vec<u8>> {
    let mut file = File::open(path)?;
    let mut raw = Vec::new();
    file.read_to_end(&mut raw)?;

    let mut seq = Vec::with_capacity(raw.len());
    let mut i = 0;
    let n = raw.len();
    let mut started_record = false;
    while i < n {
        // Read one line.
        let line_start = i;
        while i < n && raw[i] != b'\n' {
            i += 1;
        }
        let mut line_end = i;
        // strip trailing '\r'
        if line_end > line_start && raw[line_end - 1] == b'\r' {
            line_end -= 1;
        }
        if i < n {
            i += 1; // skip '\n'
        }
        let line = &raw[line_start..line_end];
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if started_record {
                // SeqIO.read reads a single record; stop at the next header.
                break;
            }
            started_record = true;
            continue;
        }
        seq.extend_from_slice(line);
    }
    Ok(seq)
}

// ---------------------------------------------------------------------------
// CLI parsing
// ---------------------------------------------------------------------------
struct Args {
    fasta: String,
    rna: String,
    max_mismatches: i32,
    max_dna_gaps: i32,
    max_rna_gaps: i32,
    chrom: String,
    output: String,
    pam: String,
    pam_5prime: bool,
}

fn parse_args() -> Args {
    let mut fasta: Option<String> = None;
    let mut rna: Option<String> = None;
    let mut max_mismatches: i32 = 4;
    let mut max_dna_gaps: i32 = 1;
    let mut max_rna_gaps: i32 = 1;
    let mut chrom: Option<String> = None;
    let mut output: Option<String> = None;
    let mut pam = String::new();
    let mut pam_5prime = false;

    let argv: Vec<String> = std::env::args().collect();
    let mut i = 1;
    while i < argv.len() {
        let a = argv[i].as_str();
        macro_rules! next_val {
            () => {{
                i += 1;
                if i >= argv.len() {
                    eprintln!("missing value for {}", a);
                    std::process::exit(2);
                }
                argv[i].clone()
            }};
        }
        match a {
            "--fasta" => fasta = Some(next_val!()),
            "--rna" => rna = Some(next_val!()),
            "--max-mismatches" => {
                max_mismatches = next_val!().parse().expect("invalid --max-mismatches")
            }
            "--max-dna-gaps" => {
                max_dna_gaps = next_val!().parse().expect("invalid --max-dna-gaps")
            }
            "--max-rna-gaps" => {
                max_rna_gaps = next_val!().parse().expect("invalid --max-rna-gaps")
            }
            "--chrom" => chrom = Some(next_val!()),
            "--output" => output = Some(next_val!()),
            "--pam" => pam = next_val!(),
            "--pam-5prime" => pam_5prime = true,
            other => {
                eprintln!("unknown argument: {}", other);
                std::process::exit(2);
            }
        }
        i += 1;
    }

    Args {
        fasta: fasta.expect("--fasta required"),
        rna: rna.expect("--rna required"),
        max_mismatches,
        max_dna_gaps,
        max_rna_gaps,
        chrom: chrom.expect("--chrom required"),
        output: output.expect("--output required"),
        pam,
        pam_5prime,
    }
}

fn main() {
    let args = parse_args();

    let dna = read_single_fasta(&args.fasta).expect("failed to read FASTA");
    let rna = args.rna.as_bytes().to_vec();
    // pam=args.pam.upper()
    let pam: Vec<u8> = args.pam.to_ascii_uppercase().into_bytes();

    let fout = File::create(&args.output).expect("failed to create output");
    let mut writer = BufWriter::new(fout);

    // Header
    writeln!(
        writer,
        "CHR\tRNA\tDNA\tStrand\tStart\tEND\tmismatches\tgaps_in_RNA\tgaps_in_DNA"
    )
    .expect("write header failed");

    // + strand
    scan(
        &rna,
        &dna,
        "+",
        args.max_dna_gaps,
        args.max_rna_gaps,
        args.max_mismatches,
        &args.chrom,
        &mut writer,
        &pam,
        args.pam_5prime,
    );

    // - strand
    let rc = reverse_complement(&dna);
    scan(
        &rna,
        &rc,
        "-",
        args.max_dna_gaps,
        args.max_rna_gaps,
        args.max_mismatches,
        &args.chrom,
        &mut writer,
        &pam,
        args.pam_5prime,
    );

    writer.flush().expect("flush failed");
}
