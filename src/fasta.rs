//! FASTA parsing with the semantics of Biopython 1.88 `FastaIO.FastaIterator`
//! (the `"fasta"` format used by `AlignIO.read`).

use anyhow::{Result, bail};

use crate::alignment::{Record, first_word, rstrip};

/// Parse FASTA text (universal newlines already applied) into records.
///
/// * the first line must start with `>` (otherwise Biopython reports comments);
/// * title = line after `>` with trailing whitespace removed; id = first word;
/// * sequence = concatenated lines with ` `, `\t`, `\r`, `\n` removed (only those).
pub fn parse(text: &str) -> Result<Vec<(Record, Vec<u8>)>> {
    let mut lines = text.split_inclusive('\n');
    let mut records = Vec::new();
    let Some(first) = lines.next() else {
        return Ok(records);
    };
    if !first.starts_with('>') {
        bail!(
            "This FASTA file contains comments at the beginning of the file, which are not \
             allowed by the 'fasta' parser."
        );
    }
    let mut title = rstrip(&first[1..]).to_string();
    let mut seq: Vec<u8> = Vec::new();
    for line in lines {
        if let Some(rest) = line.strip_prefix('>') {
            records.push(make_record(std::mem::take(&mut title), std::mem::take(&mut seq)));
            title = rstrip(rest).to_string();
            continue;
        }
        seq.extend(line.bytes().filter(|b| !matches!(b, b' ' | b'\t' | b'\r' | b'\n')));
    }
    records.push(make_record(title, seq));
    Ok(records)
}

fn make_record(title: String, seq: Vec<u8>) -> (Record, Vec<u8>) {
    let id = first_word(&title).to_string();
    (Record { id, description: title }, seq)
}
