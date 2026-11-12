//! FASTA parsing with the semantics of Biopython 1.88 `FastaIO.FastaIterator`
//! (the `"fasta"` format used by `AlignIO.read`).

use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use crate::alignment::{Alignment, AlignmentBuilder, Record, first_word, rstrip};

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

/// Read a FASTA file with the semantics of [`parse`] on its text (strict
/// UTF-8, Python's universal newlines), without holding the text or a list of
/// records: pass 1 validates UTF-8 and counts the records, pass 2 parses line
/// by line straight into an [`AlignmentBuilder`].
pub fn read(path: &Path) -> Result<Alignment> {
    let nrows = count_records(path)?;
    let file = File::open(path).with_context(|| format!("cannot read {}", path.display()))?;
    let mut reader = BufReader::with_capacity(1 << 22, file);
    let mut builder = AlignmentBuilder::new(nrows);
    let mut raw: Vec<u8> = Vec::new();
    let mut seq: Vec<u8> = Vec::new();
    let mut title: Option<String> = None;
    let mut first = true;
    let mut line_err: Option<anyhow::Error> = None;
    loop {
        raw.clear();
        if reader.read_until(b'\n', &mut raw)? == 0 {
            break;
        }
        // One physical '\n'-terminated chunk can hold several logical lines
        // (lone '\r' ends a line in universal-newline mode).
        for_each_logical_line(&raw, |line| {
            if line_err.is_some() {
                return;
            }
            if first {
                first = false;
                if !line.starts_with(b">") {
                    line_err = Some(anyhow::anyhow!(
                        "This FASTA file contains comments at the beginning of the file, which are not \
                         allowed by the 'fasta' parser."
                    ));
                    return;
                }
            }
            if let Some(rest) = line.strip_prefix(b">") {
                if let Some(t) = title.take() {
                    if let Err(e) = push_record(&mut builder, t, &seq) {
                        line_err = Some(e);
                        return;
                    }
                }
                seq.clear();
                // Pass 1 validated UTF-8, and lines never split a character.
                title = Some(rstrip(std::str::from_utf8(rest).expect("validated UTF-8")).to_string());
            } else if line.contains(&b' ') || line.contains(&b'\t') {
                seq.extend(line.iter().copied().filter(|b| !matches!(b, b' ' | b'\t' | b'\r' | b'\n')));
            } else {
                // Logical lines hold no '\r' or '\n'; nothing to remove.
                seq.extend_from_slice(line);
            }
        });
        if let Some(e) = line_err.take() {
            return Err(e);
        }
    }
    if let Some(t) = title {
        push_record(&mut builder, t, &seq)?;
    }
    builder.finish()
}

fn push_record(builder: &mut AlignmentBuilder, title: String, seq: &[u8]) -> Result<()> {
    let (record, _) = make_record(title, Vec::new());
    builder.push(record, seq)
}

/// Call `f` with each logical line of a '\n'-terminated (or final) chunk,
/// without its terminator, splitting at "\r\n", "\r" and "\n" as Python's
/// universal newlines do.
fn for_each_logical_line(raw: &[u8], mut f: impl FnMut(&[u8])) {
    let body = raw.strip_suffix(b"\n").unwrap_or(raw);
    if !body.contains(&b'\r') {
        f(body);
        return;
    }
    let body = body.strip_suffix(b"\r").unwrap_or(body);
    for line in body.split(|&b| b == b'\r') {
        f(line);
    }
}

/// Pass 1: strict UTF-8 check of the whole file (as Python's text-mode read
/// fails before parsing) and the number of records, i.e. logical lines that
/// start with '>' (at the start of the file or after '\n' or '\r').
fn count_records(path: &Path) -> Result<usize> {
    let mut file = File::open(path).with_context(|| format!("cannot read {}", path.display()))?;
    let mut buf = vec![0u8; 1 << 22];
    let mut carry = 0usize;
    let mut offset = 0usize;
    let mut prev = b'\n';
    let mut count = 0usize;
    loop {
        let n = file.read(&mut buf[carry..]).with_context(|| format!("cannot read {}", path.display()))?;
        let end = carry + n;
        if n == 0 {
            if carry > 0 {
                bail!("{} is not valid UTF-8: incomplete character at byte {offset}", path.display());
            }
            break;
        }
        let valid = match std::str::from_utf8(&buf[..end]) {
            Ok(_) => end,
            Err(e) if e.error_len().is_none() => e.valid_up_to(),
            Err(e) => bail!(
                "{} is not valid UTF-8: invalid byte at {}",
                path.display(),
                offset + e.valid_up_to()
            ),
        };
        // '>' at the start of a logical line: compare each byte with its
        // predecessor (vectorizable), carrying the last byte across reads.
        let new = &buf[carry..end];
        if let Some(&b0) = new.first() {
            count += (b0 == b'>' && (prev == b'\n' || prev == b'\r')) as usize;
            count += new[1..]
                .iter()
                .zip(new)
                .map(|(&b, &p)| (b == b'>' && (p == b'\n' || p == b'\r')) as usize)
                .sum::<usize>();
            prev = new[new.len() - 1];
        }
        // Keep an incomplete trailing character for the next read.
        buf.copy_within(valid..end, 0);
        offset += valid;
        carry = end - valid;
    }
    Ok(count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::Alignment;

    fn temp_file(name: &str, bytes: &[u8]) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!("witchi-fasta-{}-{name}", std::process::id()));
        std::fs::write(&p, bytes).unwrap();
        p
    }

    /// The streaming reader against the text pipeline it replaces
    /// (strict UTF-8 + universal newlines + `parse` + `from_records`).
    fn check(name: &str, bytes: &[u8]) {
        let path = temp_file(name, bytes);
        let streamed = read(&path);
        let via_text = std::str::from_utf8(bytes)
            .map_err(anyhow::Error::from)
            .and_then(|t| Alignment::from_records(parse(&crate::alignment::universal_newlines(t))?));
        match (streamed, via_text) {
            (Ok(a), Ok(b)) => {
                assert_eq!(a.records, b.records, "{name}");
                assert_eq!((a.ncols, &a.columns), (b.ncols, &b.columns), "{name}");
                assert_eq!(a.present, b.present, "{name}");
            }
            (Err(_), Err(_)) => {}
            (a, b) => panic!("{name}: streamed {:?} vs text {:?}", a.map(|x| x.ncols), b.map(|x| x.ncols)),
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn streaming_reader_matches_text_parser() {
        check("mixed_newlines", b">a one\r\nAC-G\rTT\n>b\rACGTTX\r\n");
        check("cr_then_header", b">a\rACGT\r>b x\rAC\rGT\r");
        check("trailing_cr", b">a\nAC\r\r");
        check("blank_and_ws", b">a  \t\n\nA C\tG\n\n>b\nACG\n");
        check("unicode_title", ">é ü\nACGT\n>ß\nACGT".as_bytes());
        check("first_line_blank", b"\n>a\nACGT\n");
        check("unequal", b">a\nACGT\n>b\nACG\n");
        check("lowercase", b">a\nacgtn\n>b\nACGTN\n");
        check("header_only", b">a\n");
        check("empty", b"");
        // A two-byte character straddling the 4 MiB read buffer of pass 1.
        let mut big = b">a ".to_vec();
        big.resize((1 << 22) - 1, b'x');
        big.extend("é\n".as_bytes());
        big.extend(b"ACGT\n>b\nACGT\n");
        check("utf8_across_buffer", &big);
    }

    #[test]
    fn streaming_reader_rejects_bad_utf8() {
        for (name, bytes) in [("invalid", &b">a\xff\nACGT\n"[..]), ("truncated", &b">a\nACGT\n>\xc3"[..])] {
            let path = temp_file(name, bytes);
            assert!(read(&path).is_err(), "{name}");
            std::fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn builder_transposes_across_tiles() {
        // More rows than one tile, and columns spanning several blocks.
        let (t_n, l) = (600usize, 150usize);
        let rows: Vec<Vec<u8>> = (0..t_n)
            .map(|t| (0..l).map(|j| b"ACDEFGHIKLMNPQRSTVWY-X"[(t * 7 + j * 13) % 22]).collect())
            .collect();
        let entries = rows
            .iter()
            .enumerate()
            .map(|(t, r)| (Record { id: format!("s{t}"), description: String::new() }, r.clone()))
            .collect();
        let aln = Alignment::from_records(entries).unwrap();
        for (t, r) in rows.iter().enumerate() {
            assert_eq!(&aln.row(t), r, "row {t}");
        }
        // Lossless encoding: decoding the codes gives the residues back.
        let mut aln = aln;
        let (enc, dec) = crate::chi2::Encoded::from_alignment(&mut aln, crate::alignment::AMINO_ACIDS);
        let lut = crate::chi2::lookup_table(crate::alignment::AMINO_ACIDS);
        for (t, r) in rows.iter().enumerate() {
            for j in 0..l {
                let code = enc.column(j)[t];
                assert_eq!(dec.byte(code), r[j]);
                assert_eq!(code >= 0, lut[r[j] as usize] >= 0);
                if code >= 0 {
                    assert_eq!(code, lut[r[j] as usize]);
                }
            }
        }
    }
}
