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

use crate::pyfmt::{py_rstrip, py_split};
use crate::{fasta, nexus};

/// Identification of a sequence record, as Biopython's SeqRecord holds it.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub id: String,
    pub description: String,
}

/// A T x L alignment of upper-cased ASCII residues, stored row-major.
#[derive(Debug, Clone)]
pub struct Alignment {
    pub records: Vec<Record>,
    pub ncols: usize,
    /// Row-major residues: `data[t * ncols + j]`.
    pub data: Vec<u8>,
}

impl Alignment {
    pub fn nrows(&self) -> usize {
        self.records.len()
    }

    pub fn row(&self, t: usize) -> &[u8] {
        &self.data[t * self.ncols..(t + 1) * self.ncols]
    }

    /// Build from (record, raw sequence) pairs: checks equal lengths
    /// (`MultipleSeqAlignment`), ASCII-only residues, and upper-cases.
    pub fn from_records(entries: Vec<(Record, Vec<u8>)>) -> Result<Self> {
        if entries.is_empty() {
            bail!("No records found in handle");
        }
        let ncols = entries[0].1.len();
        if ncols == 0 {
            // numpy cannot upper-case the resulting float array (TypeError in Python).
            bail!("alignment has zero columns");
        }
        let mut data = Vec::with_capacity(entries.len() * ncols);
        let mut records = Vec::with_capacity(entries.len());
        for (record, seq) in entries {
            if seq.len() != ncols {
                bail!("Sequences must all be the same length");
            }
            if let Some(&b) = seq.iter().find(|b| !b.is_ascii()) {
                bail!(
                    "sequence {:?} contains a non-ASCII byte 0x{b:02x}; only ASCII alignments are supported",
                    record.id
                );
            }
            data.extend(seq.iter().map(u8::to_ascii_uppercase));
            records.push(record);
        }
        Ok(Alignment { records, ncols, data })
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
pub fn read_alignment(path: &Path, format: &str) -> Result<Alignment> {
    if format != format.to_lowercase() {
        bail!("Format string '{format}' should be lower case");
    }
    let text = read_text(path)?;
    parse_alignment(&text, format)
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
    let mut present = [false; 256];
    for &b in &aln.data {
        present[b as usize] = true;
    }
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
    let mut out = Vec::with_capacity(aln.data.len() + aln.data.len() / 60 + aln.nrows() * 32);
    render_fasta(&aln.records, (0..aln.nrows()).map(|t| aln.row(t)), &mut out);
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
