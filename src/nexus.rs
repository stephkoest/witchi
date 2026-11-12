//! NEXUS alignment parsing with the semantics of Biopython 1.88
//! `AlignIO.read(..., "nexus")` (`Bio.Nexus.Nexus` using the `cnexus` scanner).
//!
//! Supported: DATA / CHARACTERS / TAXA blocks with DIMENSIONS, FORMAT
//! (datatype dna/nucleotide/rna/protein, missing, gap, matchchar, interleave,
//! respectcase) and MATRIX (interleaved or not, quoted labels, duplicate labels,
//! parenthesised IUPAC ambiguities, matchchar). Commands that do not affect the
//! matrix are accepted and ignored; TREES / SETS / CODONS blocks are ignored.
//! Unsupported: `datatype=standard`, more than one MATRIX command.

use anyhow::{Result, anyhow, bail};
use std::collections::HashMap;

use crate::pyfmt::{py_split, py_strip};

const KNOWN_BLOCKS: &[&str] = &["trees", "data", "characters", "taxa", "sets", "codons"];
const MATRIX_BLOCKS: &[&str] = &["data", "characters", "taxa"];
const SPECIAL_COMMANDS: &[&str] = &[
    "charstatelabels", "charlabels", "taxlabels", "taxset", "charset", "charpartition",
    "taxpartition", "matrix", "tree", "utree", "translate", "codonposset", "title",
];
/// Commands Biopython implements that cannot change the parsed matrix.
const IGNORED_COMMANDS: &[&str] = &[
    "title", "link", "set", "options", "eliminate", "taxlabels", "charlabels",
    "charstatelabels", "statelabels", "translate", "tree", "utree", "taxset", "charset",
    "taxpartition", "codonposset", "codeset", "charpartition", "stateset", "changeset",
    "treeset", "treepartition",
];
const PUNCTUATION: &str = "()[]{}\\,;:=*\\'\"`+-<>";
const WHITESPACE: &str = " \t\n";

/// Parse NEXUS text (universal newlines applied) into (taxon id, sequence).
pub fn parse(text: &str) -> Result<Vec<(String, String)>> {
    let mut contents = py_strip(text);
    if let Some(rest) = contents.strip_prefix("#NEXUS") {
        contents = rest;
    }
    if contents.contains('\0') {
        bail!("embedded null character in NEXUS file");
    }
    let scanned = scanfile(contents)?;
    let mut commandlines = adjust_lines(scanned.split('\u{7}'));
    for cl in commandlines.iter_mut() {
        let head: String = cl.chars().take(6).collect();
        if head.to_uppercase() == "#NEXUS" {
            let rest: String = cl.chars().skip(6).collect();
            *cl = py_strip(&rest).to_string();
        }
    }

    let mut state = State::default();
    for (title, lines) in blocks(&commandlines)? {
        if !KNOWN_BLOCKS.contains(&title.as_str()) {
            continue;
        }
        let commands = lines
            .iter()
            .map(|l| Commandline::parse(l))
            .collect::<Result<Vec<_>>>()?;
        if !MATRIX_BLOCKS.contains(&title.as_str()) {
            continue;
        }
        for cmd in commands {
            state.execute(cmd)?;
        }
    }
    state.into_records()
}

/// Port of `cnexus.scanfile`: drop `[...]` comments (nested), keep `[&...]`
/// special comments, track quotes, and replace command-ending `;` by `\x07`.
fn scanfile(input: &str) -> Result<String> {
    let bytes = input.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut quotelevel = 0u8;
    let mut speciallevel = false;
    let mut commlevel: i64 = 0;
    for (i, &t) in bytes.iter().enumerate() {
        if commlevel == 0 && !speciallevel && quotelevel != 0 && t == quotelevel {
            quotelevel = 0;
        } else if quotelevel == 0 && commlevel == 0 && !speciallevel && (t == b'\'' || t == b'"') {
            quotelevel = t;
        } else if quotelevel == 0 && t == b'[' {
            if bytes.get(i + 1) == Some(&b'&') && commlevel == 0 && !speciallevel {
                speciallevel = true;
            } else {
                commlevel += 1;
            }
        } else if quotelevel == 0 && t == b']' {
            if speciallevel {
                speciallevel = false;
            } else {
                commlevel -= 1;
                if commlevel < 0 {
                    bail!("Unmatched ]");
                }
                continue;
            }
        }
        if commlevel == 0 {
            if t == b';' && quotelevel == 0 && !speciallevel {
                out.push(7);
            } else {
                out.push(t);
            }
        }
    }
    if commlevel > 0 {
        bail!("Unmatched [");
    }
    // Comment boundaries are ASCII, so whole UTF-8 sequences are kept or dropped.
    String::from_utf8(out).map_err(|e| anyhow!("NEXUS scan produced invalid UTF-8: {e}"))
}

/// Port of `_adjust_lines`.
fn adjust_lines<'a>(lines: impl Iterator<Item = &'a str>) -> Vec<String> {
    let mut formatted = Vec::new();
    for line in lines {
        let line = line.replace("\r\n", "\n").replace('\r', "\n");
        let line = py_strip(&line);
        if line.to_lowercase().starts_with("matrix") {
            formatted.push(line.to_string());
        } else {
            let line = line.replace('\n', " ");
            if !line.is_empty() {
                formatted.push(line);
            }
        }
    }
    formatted
}

/// Port of `Nexus._get_nexus_block`: (lower-cased block title, command lines).
fn blocks(commandlines: &[String]) -> Result<Vec<(String, Vec<String>)>> {
    let mut out = Vec::new();
    let mut current: Option<(String, Vec<String>)> = None;
    for cl in commandlines {
        let lower = cl.to_lowercase();
        if lower.starts_with("begin") {
            if let Some((title, _)) = &current {
                bail!("Illegal block nesting in block {title}");
            }
            let title = py_split(cl)
                .nth(1)
                .ok_or_else(|| anyhow!("NEXUS 'begin' without a block name"))?
                .to_lowercase();
            current = Some((title, Vec::new()));
        } else if lower.starts_with("end") {
            match current.take() {
                Some(block) => out.push(block),
                None => bail!("Unmatched 'end'."),
            }
        } else if let Some((_, lines)) = current.as_mut() {
            lines.push(cl.clone());
        }
    }
    Ok(out)
}

/// Port of `Commandline`.
struct Commandline {
    command: String,
    /// Raw option text for special commands.
    text: String,
    /// Parsed `key[=value]` options for the other commands.
    options: HashMap<String, Option<String>>,
}

impl Commandline {
    fn parse(line: &str) -> Result<Self> {
        let stripped = py_strip(line);
        let (command, options) = match stripped.split_once('\n') {
            Some((c, o)) => (c.to_string(), o.to_string()),
            None => {
                let words: Vec<&str> = py_split(line).collect();
                let first = words.first().ok_or_else(|| anyhow!("empty NEXUS command"))?;
                (first.to_string(), words[1..].join(" "))
            }
        };
        let command = py_strip(&command).to_lowercase();
        let mut cl = Commandline { command, text: String::new(), options: HashMap::new() };
        if SPECIAL_COMMANDS.contains(&cl.command.as_str()) {
            cl.text = py_strip(&options).to_string();
        } else if !options.is_empty() {
            let replaced = options.replace('=', " = ");
            let tokens: Vec<&str> = py_split(&replaced).collect();
            let mut in_pair = vec![false; tokens.len()];
            let mut pairs = Vec::new();
            for n in 1..tokens.len() {
                if tokens[n] == "=" {
                    if n + 1 >= tokens.len() {
                        bail!("Incorrect formatting in NEXUS line (dangling '='): {line}");
                    }
                    pairs.push((tokens[n - 1], tokens[n + 1]));
                    in_pair[n - 1] = true;
                    in_pair[n] = true;
                    in_pair[n + 1] = true;
                }
            }
            for (k, v) in pairs {
                cl.options.insert(k.to_lowercase(), Some(v.to_string()));
            }
            for (tok, used) in tokens.iter().zip(&in_pair) {
                if !used {
                    cl.options.insert(tok.to_lowercase(), None);
                }
            }
        }
        Ok(cl)
    }
}

#[derive(Default)]
struct State {
    ntax: i64,
    nchar: i64,
    datatype: Option<String>, // None = Biopython's default "dna"
    respectcase: bool,
    missing: Option<char>,
    gap: Option<char>,
    matchchar: Option<char>,
    interleave: bool,
    /// Set once a FORMAT command with a datatype has run (Biopython's
    /// `ambiguous_values` / `valid_characters` / `rev_ambiguous_values`).
    alphabet: Option<Alphabet>,
    matrix_seen: bool,
    /// Taxon label -> sequence, in insertion order.
    matrix: Vec<(String, String)>,
    taxlabels: Vec<String>,
}

struct Alphabet {
    valid_characters: Vec<char>,
    /// Sorted ambiguity set (e.g. "AC") -> IUPAC code (e.g. "M").
    rev_ambiguous: HashMap<String, String>,
}

impl State {
    fn missing(&self) -> char {
        self.missing.unwrap_or('?')
    }
    fn gap(&self) -> char {
        self.gap.unwrap_or('-')
    }
    fn datatype(&self) -> &str {
        self.datatype.as_deref().unwrap_or("dna")
    }

    fn execute(&mut self, cmd: Commandline) -> Result<()> {
        match cmd.command.as_str() {
            "dimensions" => {
                if let Some(v) = cmd.options.get("ntax") {
                    self.ntax = py_int(v.as_deref())?;
                }
                if let Some(v) = cmd.options.get("nchar") {
                    self.nchar = py_int(v.as_deref())?;
                }
                Ok(())
            }
            "format" => self.format(&cmd.options),
            "matrix" => self.parse_matrix(&cmd.text),
            c if IGNORED_COMMANDS.contains(&c) => Ok(()),
            c => bail!("Unknown command: {c} "),
        }
    }

    fn format(&mut self, opts: &HashMap<String, Option<String>>) -> Result<()> {
        let value = |key: &str| -> Result<Option<&str>> {
            match opts.get(key) {
                None => Ok(None),
                Some(None) => bail!("NEXUS FORMAT option '{key}' needs a value"),
                Some(Some(v)) => Ok(Some(v.as_str())),
            }
        };
        if opts.contains_key("respectcase") {
            self.respectcase = true;
        }
        value("symbols")?;
        if let Some(dt) = value("datatype")? {
            let dt = dt.to_lowercase();
            let (ambiguous, unambiguous): (Vec<(&str, &str)>, String) = match dt.as_str() {
                "dna" | "nucleotide" => (AMBIGUOUS_DNA.to_vec(), "GATC".to_string()),
                "rna" => (AMBIGUOUS_RNA.to_vec(), "GAUC".to_string()),
                "protein" => (
                    vec![("B", "DN"), ("Z", "EQ"), ("X", PROTEIN_LETTERS)],
                    format!("{PROTEIN_LETTERS}*"),
                ),
                "standard" => bail!("NEXUS datatype=standard is not supported"),
                other => bail!("Unsupported datatype: {other}"),
            };
            let mut valid: String = ambiguous.iter().map(|(k, _)| *k).collect();
            valid.push_str(&unambiguous);
            if !self.respectcase {
                valid = format!("{}{}", valid.to_lowercase(), valid.to_uppercase());
            }
            let mut rev: Vec<(&str, &str)> = Vec::new();
            for &(k, v) in &ambiguous {
                if k == "X" {
                    continue;
                }
                match rev.iter_mut().find(|(rv, _)| *rv == v) {
                    Some(entry) => entry.1 = k,
                    None => rev.push((v, k)),
                }
            }
            let rev_ambiguous = rev
                .into_iter()
                .map(|(v, k)| {
                    let mut chars: Vec<char> = v.chars().collect();
                    chars.sort_unstable();
                    (chars.into_iter().collect::<String>(), k.to_string())
                })
                .collect();
            self.alphabet = Some(Alphabet { valid_characters: valid.chars().collect(), rev_ambiguous });
            self.datatype = Some(dt);
        }
        if matches!(self.datatype(), "dna" | "nucleotide" | "rna") && self.alphabet.is_none() {
            // Biopython touches self.ambiguous_values here, which only exists
            // after a datatype was given: AttributeError -> "Unknown command".
            bail!("Unknown command: format ");
        }
        if let Some(v) = value("missing")? {
            self.missing = v.chars().next();
        }
        if let Some(v) = value("gap")? {
            self.gap = v.chars().next();
        }
        if let Some(v) = value("matchchar")? {
            self.matchchar = v.chars().next();
        }
        if let Some(v) = opts.get("interleave") {
            if v.as_deref().is_none_or(|s| s.to_lowercase() == "yes") {
                self.interleave = true;
            }
        }
        Ok(())
    }

    fn parse_matrix(&mut self, options: &str) -> Result<()> {
        if self.ntax == 0 || self.nchar == 0 {
            bail!("Dimensions must be specified before matrix!");
        }
        if self.ntax < 0 || self.nchar < 0 {
            bail!("negative NEXUS dimensions");
        }
        if self.matrix_seen {
            bail!("more than one NEXUS MATRIX command is not supported");
        }
        self.matrix_seen = true;
        let alphabet = self
            .alphabet
            .as_ref()
            .ok_or_else(|| anyhow!("Unknown command: matrix "))?;
        let ntax = self.ntax as usize;
        let nchar = self.nchar as usize;
        let (gap, missing) = (self.gap(), self.missing());

        let lines: Vec<&str> =
            options.split('\n').map(py_strip).filter(|l| !l.is_empty()).collect();
        let mut it = lines.into_iter();
        let mut taxcount = 0usize;
        let mut first_block = true;
        let mut refseq: Vec<char> = Vec::new();
        loop {
            let Some(line) = it.next() else {
                if taxcount < ntax {
                    bail!("Not enough taxa in matrix.");
                }
                break;
            };
            taxcount += 1;
            if taxcount > ntax {
                if !self.interleave {
                    bail!("Too many taxa in matrix - should matrix be interleaved?");
                }
                taxcount = 1;
                first_block = false;
            }
            let mut buf = CharBuffer::new(line);
            let id = quotestrip(buf.next_word()?).unwrap_or_default();
            let rest = buf.rest();
            let rest = py_strip(&rest);
            let mut chars: String;
            if self.interleave {
                if !rest.is_empty() {
                    chars = py_split(rest).collect();
                } else {
                    let next = it.next().ok_or_else(|| anyhow!("NEXUS matrix ended early"))?;
                    chars = py_split(next).collect();
                }
            } else {
                chars = py_split(rest).collect();
                while chars.chars().count() < nchar {
                    let next = it.next().ok_or_else(|| anyhow!("NEXUS matrix ended early"))?;
                    chars.extend(py_split(next));
                }
            }

            let mut seq: Vec<char> =
                replace_parenthesized_ambigs(&chars, &alphabet.rev_ambiguous)?.chars().collect();
            if taxcount == 1 {
                refseq = seq.clone();
            } else if let Some(m) = self.matchchar {
                for (p, c) in seq.iter_mut().enumerate() {
                    if *c == m {
                        let r = *refseq.get(p).ok_or_else(|| {
                            anyhow!("matchchar beyond the reference sequence for taxon {id}")
                        })?;
                        if r == m {
                            bail!("matchchar in the reference sequence for taxon {id}");
                        }
                        *c = r;
                    }
                }
            }
            for &c in &seq {
                if !alphabet.valid_characters.contains(&c) && c != gap && c != missing {
                    bail!(
                        "Taxon {id}: Illegal character {c} in sequence (check dimensions/interleaving)"
                    );
                }
            }
            let seq: String = seq.into_iter().collect();

            if first_block {
                let previous: Vec<&str> = self.matrix.iter().map(|(k, _)| k.as_str()).collect();
                let label = unique_label(&previous, &id)?;
                self.matrix.push((label.clone(), seq));
                self.taxlabels.push(label);
            } else {
                let previous: Vec<&str> =
                    self.taxlabels[..taxcount - 1].iter().map(String::as_str).collect();
                let label = unique_label(&previous, &id)?;
                let nexid = label.replace(' ', "_");
                let mut present: Option<&str> = None;
                for t in &self.taxlabels {
                    if t.replace(' ', "_") == nexid {
                        present = Some(t.as_str());
                    }
                }
                match present.filter(|t| !t.is_empty()) {
                    Some(t) => {
                        let t = t.to_string();
                        let entry = self.matrix.iter_mut().find(|(k, _)| *k == t).unwrap();
                        entry.1.push_str(&seq);
                    }
                    None => bail!(
                        "Taxon {label} not in first block of interleaved matrix. \
                         Check matrix dimensions and interleave."
                    ),
                }
            }
        }
        for (taxon, seq) in &self.matrix {
            let len = seq.chars().count();
            if len != nchar {
                bail!("Matrix Nchar {nchar} does not match data length ({len}) for taxon {taxon}");
            }
        }
        Ok(())
    }

    fn into_records(self) -> Result<Vec<(String, String)>> {
        if self.matrix.is_empty() {
            bail!("No records found in handle");
        }
        Ok(self.matrix)
    }
}

/// Python `int(s, 10)` for option values.
fn py_int(v: Option<&str>) -> Result<i64> {
    let s = v.ok_or_else(|| anyhow!("NEXUS dimension needs a value"))?;
    let t = py_strip(s);
    let (neg, digits) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    let valid = !digits.is_empty()
        && digits.bytes().all(|b| b.is_ascii_digit() || b == b'_')
        && !digits.starts_with('_')
        && !digits.ends_with('_')
        && !digits.contains("__");
    if !valid {
        bail!("invalid literal for int() with base 10: '{s}'");
    }
    let n: i64 = digits.replace('_', "").parse().map_err(|e| anyhow!("{e}"))?;
    Ok(if neg { -n } else { n })
}

/// Port of `quotestrip`.
fn quotestrip(word: Option<String>) -> Option<String> {
    let mut word = word?;
    if word.is_empty() {
        return None;
    }
    while word.len() >= 1
        && ((word.starts_with('\'') && word.ends_with('\''))
            || (word.starts_with('"') && word.ends_with('"')))
    {
        let n = word.chars().count();
        word = word.chars().skip(1).take(n.saturating_sub(2)).collect();
    }
    Some(word)
}

/// Port of `_unique_label`.
fn unique_label(previous: &[&str], label: &str) -> Result<String> {
    let mut label = label.to_string();
    while previous.contains(&label.as_str()) {
        let parts: Vec<&str> = label.split('.').collect();
        let last = *parts.last().unwrap();
        if last.starts_with("copy") {
            let copy_num = if last != "copy" {
                let n: Vec<char> = last.chars().skip(4).collect();
                py_int(Some(&n.into_iter().collect::<String>()))? + 1
            } else {
                1
            };
            label = format!("{}.copy{copy_num}", parts[..parts.len() - 1].join("."));
        } else {
            label.push_str(".copy");
        }
    }
    Ok(label)
}

/// Port of `_replace_parenthesized_ambigs`.
fn replace_parenthesized_ambigs(seq: &str, rev: &HashMap<String, String>) -> Result<String> {
    let mut seq: Vec<char> = seq.chars().collect();
    while let Some(opening) = seq.iter().position(|&c| c == '(') {
        let closing = seq
            .iter()
            .position(|&c| c == ')')
            .ok_or_else(|| anyhow!("Missing closing parenthesis in NEXUS sequence"))?;
        if closing < opening {
            bail!("Missing opening parenthesis in NEXUS sequence");
        }
        let mut inner: Vec<char> = seq[opening + 1..closing].to_vec();
        inner.sort_unstable();
        let ambig: String = inner.into_iter().collect();
        let mut code = rev
            .get(&ambig.to_uppercase())
            .ok_or_else(|| anyhow!("unknown ambiguity ({ambig}) in NEXUS sequence"))?
            .clone();
        if ambig != ambig.to_uppercase() {
            code = code.to_lowercase();
        }
        seq.splice(opening..=closing, code.chars());
    }
    Ok(seq.into_iter().collect())
}

/// Port of `CharBuffer` (only what `_matrix` needs).
struct CharBuffer {
    buf: Vec<char>,
    pos: usize,
}

impl CharBuffer {
    fn new(s: &str) -> Self {
        CharBuffer { buf: s.chars().collect(), pos: 0 }
    }
    fn peek(&self) -> Option<char> {
        self.buf.get(self.pos).copied()
    }
    fn next(&mut self) -> Option<char> {
        let c = self.peek();
        if c.is_some() {
            self.pos += 1;
        }
        c
    }
    fn next_nonwhitespace(&mut self) -> Option<char> {
        while let Some(c) = self.next() {
            if !WHITESPACE.contains(c) {
                return Some(c);
            }
        }
        None
    }
    fn next_word(&mut self) -> Result<Option<String>> {
        let Some(first) = self.next_nonwhitespace() else {
            return Ok(None);
        };
        let mut word = String::from(first);
        let quoted = match first {
            '\'' | '"' => Some(first),
            c if PUNCTUATION.contains(c) => return Ok(Some(word)),
            _ => None,
        };
        loop {
            let c = self.peek();
            if c.is_some() && c == quoted {
                word.push(self.next().unwrap());
                if self.peek() == quoted {
                    self.next();
                } else {
                    break;
                }
            } else if quoted.is_some() {
                match self.next() {
                    Some(ch) => word.push(ch),
                    // Biopython crashes (joins None) on an unterminated quote.
                    None => bail!("unterminated quote in NEXUS label {word}"),
                }
            } else {
                match c {
                    None => break,
                    Some(ch) if PUNCTUATION.contains(ch) || WHITESPACE.contains(ch) => break,
                    Some(ch) => {
                        word.push(ch);
                        self.next();
                    }
                }
            }
        }
        Ok(Some(word))
    }
    fn rest(&self) -> String {
        self.buf[self.pos..].iter().collect()
    }
}

const PROTEIN_LETTERS: &str = "ACDEFGHIKLMNPQRSTVWY";
const AMBIGUOUS_DNA: &[(&str, &str)] = &[
    ("A", "A"), ("C", "C"), ("G", "G"), ("T", "T"), ("M", "AC"), ("R", "AG"), ("W", "AT"),
    ("S", "CG"), ("Y", "CT"), ("K", "GT"), ("V", "ACG"), ("H", "ACT"), ("D", "AGT"),
    ("B", "CGT"), ("X", "GATC"), ("N", "GATC"),
];
const AMBIGUOUS_RNA: &[(&str, &str)] = &[
    ("A", "A"), ("C", "C"), ("G", "G"), ("U", "U"), ("M", "AC"), ("R", "AG"), ("W", "AU"),
    ("S", "CG"), ("Y", "CU"), ("K", "GU"), ("V", "ACG"), ("H", "ACU"), ("D", "AGU"),
    ("B", "CGU"), ("X", "GAUC"), ("N", "GAUC"),
];
