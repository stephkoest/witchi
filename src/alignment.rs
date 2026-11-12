//! Alignment input/output with Biopython semantics.
//!
//! `read_alignment` mirrors `AlignmentReader.run` (`Bio.AlignIO.read` followed by
//! `np.char.upper`), `detect_sequence_type` mirrors `SequenceTypeDetector.detect`
//! and `write_fasta` mirrors `AlignIO.write(..., "fasta")` of the pruned records.

use anyhow::{Context, Result, bail};
use std::borrow::Cow;
use std::fs;
use std::io::Write;
use std::path::Path;

use rayon::prelude::*;

use crate::pyfmt::{py_rstrip, py_split};
use crate::{fasta, nexus};

/// Identification of a sequence record, as Biopython's SeqRecord holds it.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub id: String,
    pub description: String,
}

/// A T x L alignment of upper-cased ASCII residues.
///
/// The residues are stored once, column-major (`columns[j * nrows + t]`), which
/// is the layout the chi-square kernels use: [`crate::chi2::Encoded::from_alignment`]
/// encodes them in place and takes them over, so a large alignment never
/// exists twice in memory.
#[derive(Debug, Clone)]
pub struct Alignment {
    pub records: Vec<Record>,
    pub ncols: usize,
    /// Column-major residues; empty once encoded.
    pub columns: Vec<u8>,
    /// Which byte values occur among the residues.
    pub present: [bool; 256],
}

impl Alignment {
    pub fn nrows(&self) -> usize {
        self.records.len()
    }

    /// Row `t` (gathered from the column-major storage; for small uses).
    pub fn row(&self, t: usize) -> Vec<u8> {
        let n = self.nrows();
        (0..self.ncols).map(|j| self.columns[j * n + t]).collect()
    }

    /// Build from (record, raw sequence) pairs: checks equal lengths
    /// (`MultipleSeqAlignment`), ASCII-only residues, and upper-cases.
    pub fn from_records(entries: Vec<(Record, Vec<u8>)>) -> Result<Self> {
        let mut builder = AlignmentBuilder::new(entries.len());
        for (record, seq) in entries {
            builder.push(record, &seq)?;
        }
        builder.finish()
    }
}

/// Rows buffered before they are transposed into the column-major matrix.
const BUILD_TILE: usize = 256;
/// Columns per parallel transpose task.
const BUILD_COLS_PER_TASK: usize = 64;

/// Builds an [`Alignment`] record by record, with the number of records
/// known up front, writing residues straight into column-major storage
/// through a small row tile.
///
/// The checks and their order match `Alignment::from_records` of a list:
/// no records, zero columns, then per record unequal length and non-ASCII.
pub struct AlignmentBuilder {
    nrows: usize,
    ncols: Option<usize>,
    records: Vec<Record>,
    columns: Vec<u8>,
    tile: Vec<u8>,
    tile_rows: usize,
}

impl AlignmentBuilder {
    pub fn new(nrows: usize) -> Self {
        AlignmentBuilder {
            nrows,
            ncols: None,
            records: Vec::with_capacity(nrows),
            columns: Vec::new(),
            tile: Vec::new(),
            tile_rows: 0,
        }
    }

    /// Add the next record and its raw (not yet upper-cased) sequence.
    pub fn push(&mut self, record: Record, seq: &[u8]) -> Result<()> {
        let ncols = match self.ncols {
            Some(l) => l,
            None => {
                if seq.is_empty() {
                    // numpy cannot upper-case the resulting float array (TypeError in Python).
                    bail!("alignment has zero columns");
                }
                let l = seq.len();
                self.ncols = Some(l);
                // Zero-filled pages are only touched as they are written.
                self.columns = vec![0u8; self.nrows * l];
                self.tile = vec![0u8; BUILD_TILE.min(self.nrows) * l];
                l
            }
        };
        if seq.len() != ncols {
            bail!("Sequences must all be the same length");
        }
        if !seq.is_ascii() {
            let b = *seq.iter().find(|b| !b.is_ascii()).expect("non-ASCII byte");
            bail!(
                "sequence {:?} contains a non-ASCII byte 0x{b:02x}; only ASCII alignments are supported",
                record.id
            );
        }
        if self.records.len() >= self.nrows {
            bail!("more records than counted ({})", self.nrows);
        }
        let row = &mut self.tile[self.tile_rows * ncols..(self.tile_rows + 1) * ncols];
        row.copy_from_slice(seq);
        row.make_ascii_uppercase();
        self.records.push(record);
        self.tile_rows += 1;
        if self.tile_rows == BUILD_TILE {
            self.flush();
        }
        Ok(())
    }

    /// Transpose the buffered rows into the column-major matrix, in parallel
    /// over blocks of columns (each block is a contiguous part of the matrix).
    fn flush(&mut self) {
        let (Some(ncols), true) = (self.ncols, self.tile_rows > 0) else {
            return;
        };
        let (t_n, rows) = (self.nrows, self.tile_rows);
        let t0 = self.records.len() - rows;
        let tile = &self.tile;
        self.columns.par_chunks_mut(t_n * BUILD_COLS_PER_TASK).enumerate().for_each(|(blk, dst)| {
            let j0 = blk * BUILD_COLS_PER_TASK;
            let nj = dst.len() / t_n;
            for r in 0..rows {
                let src = &tile[r * ncols + j0..r * ncols + j0 + nj];
                for (jj, &b) in src.iter().enumerate() {
                    dst[jj * t_n + t0 + r] = b;
                }
            }
        });
        self.tile_rows = 0;
    }

    pub fn finish(mut self) -> Result<Alignment> {
        if self.records.is_empty() {
            bail!("No records found in handle");
        }
        if self.records.len() != self.nrows {
            bail!("expected {} records, read {}", self.nrows, self.records.len());
        }
        self.flush();
        // The byte values present, in parallel over the matrix.
        let present = self
            .columns
            .par_chunks(1 << 20)
            .map(|chunk| {
                let mut p = [false; 256];
                for &b in chunk {
                    p[b as usize] = true;
                }
                p
            })
            .reduce(|| [false; 256], |mut a, b| {
                for (x, y) in a.iter_mut().zip(b) {
                    *x |= y;
                }
                a
            });
        Ok(Alignment {
            records: self.records,
            ncols: self.ncols.expect("at least one record"),
            columns: self.columns,
            present,
        })
    }
}

/// Python text-mode reading: strict UTF-8 with universal newlines.
pub fn read_text(path: &Path) -> Result<String> {
    let bytes = fs::read(path).with_context(|| format!("cannot read {}", path.display()))?;
    let text = String::from_utf8(bytes)
        .map_err(|e| anyhow::anyhow!("{} is not valid UTF-8: {e}", path.display()))?;
    Ok(match universal_newlines(&text) {
        Cow::Borrowed(_) => text,
        Cow::Owned(s) => s,
    })
}

/// Translate `\r\n` and lone `\r` to `\n`, as Python's text mode does.
pub fn universal_newlines(text: &str) -> Cow<'_, str> {
    if !text.contains('\r') {
        return Cow::Borrowed(text);
    }
    Cow::Owned(text.replace("\r\n", "\n").replace('\r', "\n"))
}

/// Read an alignment file (`Bio.AlignIO.read(file, format)` + upper-casing).
///
/// FASTA is streamed from the file (never held as text); NEXUS is parsed
/// from the whole text.
pub fn read_alignment(path: &Path, format: &str) -> Result<Alignment> {
    if format != format.to_lowercase() {
        bail!("Format string '{format}' should be lower case");
    }
    match format {
        "fasta" => fasta::read(path),
        _ => parse_alignment(&read_text(path)?, format),
    }
}

/// Parse alignment text that already has universal newlines applied.
pub fn parse_alignment(text: &str, format: &str) -> Result<Alignment> {
    let entries = match format {
        "fasta" => fasta::parse(text)?,
        "nexus" => nexus::parse(text)?
            .into_iter()
            .map(|(id, seq)| (Record { id, description: String::new() }, seq.into_bytes()))
            .collect(),
        other => bail!(
            "Unsupported alignment format '{other}': this port supports 'fasta' and 'nexus'"
        ),
    };
    Alignment::from_records(entries)
}

pub const DNA_BASES: &[u8] = b"ACGT";
pub const AMINO_ACIDS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

/// Outcome of `SequenceTypeDetector.detect`.
#[derive(Debug, Clone, PartialEq)]
pub enum SequenceType {
    Dna,
    Protein,
    /// Characters outside both alphabets: treated as amino acids.
    UnexpectedProtein(Vec<u8>),
}

impl SequenceType {
    pub fn char_set(&self) -> &'static [u8] {
        match self {
            SequenceType::Dna => DNA_BASES,
            _ => AMINO_ACIDS,
        }
    }

    /// The message Python prints for this outcome.
    pub fn message(&self) -> String {
        match self {
            SequenceType::Dna => "Detected DNA sequence.".to_string(),
            SequenceType::Protein => "Detected amino acid sequence.".to_string(),
            SequenceType::UnexpectedProtein(chars) => {
                let list: Vec<String> =
                    chars.iter().map(|&c| format!("'{}'", c as char)).collect();
                format!(
                    "⚠️  Sequence contains unexpected characters: [{}]\nDefaulting to amino acid interpretation.",
                    list.join(", ")
                )
            }
        }
    }
}

/// `SequenceTypeDetector.detect` on the upper-cased alignment.
pub fn detect_sequence_type(aln: &Alignment) -> Result<SequenceType> {
    let mut present = aln.present;
    for &b in b"N-?X." {
        present[b as usize] = false;
    }
    let chars: Vec<u8> = (0u8..=255).filter(|&b| present[b as usize]).collect();
    if chars.is_empty() {
        bail!("No valid sequence characters detected.");
    }
    if chars.iter().all(|c| DNA_BASES.contains(c)) {
        Ok(SequenceType::Dna)
    } else if chars.iter().all(|c| AMINO_ACIDS.contains(c)) {
        Ok(SequenceType::Protein)
    } else {
        Ok(SequenceType::UnexpectedProtein(chars))
    }
}

/// Biopython FASTA title for a record (`FastaWriter.write_record`).
pub fn fasta_title(record: &Record) -> String {
    let clean = |s: &str| s.replace(['\n', '\r'], " ");
    let id = clean(&record.id);
    let description = clean(&record.description);
    if !description.is_empty() && py_split(&description).next() == Some(id.as_str()) {
        description
    } else if !description.is_empty() {
        format!("{id} {description}")
    } else {
        id
    }
}

/// Render records as Biopython FASTA (60-column wrapping) into a byte buffer.
pub fn render_fasta<'a>(
    records: &[Record],
    rows: impl Iterator<Item = &'a [u8]>,
    out: &mut Vec<u8>,
) {
    for (record, seq) in records.iter().zip(rows) {
        out.push(b'>');
        out.extend_from_slice(fasta_title(record).as_bytes());
        out.push(b'\n');
        for chunk in seq.chunks(60) {
            out.extend_from_slice(chunk);
            out.push(b'\n');
        }
    }
}

/// Write the alignment's rows as FASTA to `path`.
pub fn write_fasta(path: &Path, aln: &Alignment) -> Result<()> {
    let rows: Vec<Vec<u8>> = (0..aln.nrows()).map(|t| aln.row(t)).collect();
    let mut out = Vec::with_capacity(aln.columns.len() + aln.columns.len() / 60 + aln.nrows() * 32);
    render_fasta(&aln.records, rows.iter().map(Vec::as_slice), &mut out);
    let mut f = fs::File::create(path).with_context(|| format!("cannot create {}", path.display()))?;
    f.write_all(&out)?;
    Ok(())
}

/// First word of a title as Biopython computes the record id.
pub fn first_word(title: &str) -> &str {
    py_split(title).next().unwrap_or("")
}

/// Strip trailing Python whitespace (re-exported for the FASTA parser).
pub fn rstrip(s: &str) -> &str {
    py_rstrip(s)
}
