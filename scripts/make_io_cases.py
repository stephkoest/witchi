#!/usr/bin/env python3
"""Write the hand-crafted FASTA/NEXUS edge-case inputs to tests/io/inputs/.

Usage: python -I scripts/make_io_cases.py OUT_DIR
Files are written as raw bytes so CRLF/CR/odd whitespace survive exactly.
"""

import sys
from pathlib import Path

CASES = {
    # ---------------- FASTA ----------------
    "crlf.fa": b">a desc\r\nACGT\r\nAC\r\n>b\r\nTTGGCA\r\n",
    "cr_only.fa": b">a\rACGTAC\r>b\rTTGGCA\r",
    "leading_blank.fa": b"\n>a\nACGT\n",
    "comment_first.fa": b";comment\n>a\nACGT\n",
    "empty.fa": b"",
    "unequal.fa": b">a\nACGT\n>b\nACG\n",
    "empty_title.fa": b">\nACGTAA\n>b\nACGTTT\n",
    "spaces_tabs.fa": b">id1\tsome desc  \t\nAC GT\tAC \n  GT\n>id2 x\nACGTACGT\n",
    "vt_in_seq.fa": b">a\nAC\x0bGT\n>b\nACGTT\n",
    "lowercase_dna.fa": b">a\nacgtn\n>b\nACGGT\n>c\nacg-t\n",
    "all_ambiguous.fa": b">a\nNN-?\n>b\nX.nx\n",
    "semicolon_line.fa": b">a\nACDE\n;x\n>b\nACDEFG\n",
    "unicode_header.fa": ">séq desc　\nMKV\n>t2 ✓\nMKA\n".encode(),
    "no_trailing_newline.fa": b">a\nMKVL\n>b\nMKAL",
    "single.fa": b">only one\nACGTTGCA\n",
    "zero_length.fa": b">a\n>b\n",
    "blank_lines.fa": b">a\n\nAC\n\nGT\n\n>b\nACGA\n\n",
    "header_only_gt.fa": b">\n>\n",
    # ---------------- NEXUS ----------------
    "simple_protein.nex": b"""#NEXUS
begin data;
  dimensions ntax=3 nchar=12;
  format datatype=protein missing=? gap=-;
  matrix
  [ a comment [nested] here ]
  alpha MKVLAA
        GHIKLM
  beta  MKV-AAGHIKL?
  gamma mkvlaaghiklm
  ;
end;
""",
    "quoted_dup.nex": b"""#NEXUS
BEGIN DATA;
DIMENSIONS NTAX=6 NCHAR=8;
FORMAT DATATYPE=DNA GAP=- MISSING=?;
MATRIX
'Homo sapiens' ACGTACGT
'O''Brien'     ACGTACGA
"double q"     ACGTACGG
dup            ACGTACGC
dup            ACGTAAGC
dup            ACGTAAGG
;
END;
""",
    "ambig.nex": b"""#NEXUS
begin data; dimensions ntax=2 nchar=6; format datatype=dna;
matrix
t1 A(AC)GT(ag)T
t2 ACGT(CGT)T
;
end;
""",
    "rna.nex": b"""#NEXUS
begin data; dimensions ntax=2 nchar=5; format datatype=rna;
matrix
r1 ACGUU
r2 acguN
;
end;
""",
    "crlf.nex": b"#NEXUS\r\nbegin data;\r\ndimensions ntax=2 nchar=4;\r\nformat datatype=dna interleave;\r\nmatrix\r\nx AC\r\ny AG\r\nx GT\r\ny GA\r\n;\r\nend;\r\n",
    "interleave_split.nex": b"""#NEXUS
begin data;
dimensions ntax=2 nchar=8;
format datatype=dna interleave=yes matchchar=.;
matrix
first
ACGT
second
..GA

first ACCC
second .T..
;
end;
""",
    "other_blocks.nex": b"""#nexus
[ leading comment ]
begin taxa; dimensions ntax=2; taxlabels a b; end;
begin characters; dimensions nchar=5; format datatype=protein; matrix
a MKVLA
b MKVLG
; end;
begin trees; tree t1 = [&U] (a,b); end;
begin sets; charset first = 1-2; end;
begin assumptions; options deftype=unord; end;
""",
    "custom_gap.nex": b"""#NEXUS
begin data; dimensions ntax=2 nchar=4; format datatype=dna gap=~ missing=*;
matrix
a AC~*
b ACGT
;
end;
""",
    "missing_datatype.nex": b"#NEXUS\nbegin data; dimensions ntax=1 nchar=2; format gap=-; matrix\na AC\n; end;\n",
    "no_format.nex": b"#NEXUS\nbegin data; dimensions ntax=1 nchar=2; matrix\na AC\n; end;\n",
    "wrong_nchar.nex": b"#NEXUS\nbegin data; dimensions ntax=2 nchar=5; format datatype=dna; matrix\na ACGTA\nb ACGT\n; end;\n",
    "unknown_cmd.nex": b"#NEXUS\nbegin data; dimensions ntax=1 nchar=2; format datatype=dna; frobnicate x=1; matrix\na AC\n; end;\n",
    "illegal_char.nex": b"#NEXUS\nbegin data; dimensions ntax=1 nchar=2; format datatype=dna; matrix\na AJ\n; end;\n",
    "too_many_taxa.nex": b"#NEXUS\nbegin data; dimensions ntax=1 nchar=2; format datatype=dna; matrix\na AC\nb AC\n; end;\n",
    "unmatched_bracket.nex": b"#NEXUS\nbegin data; [ unclosed\n end;\n",
    "no_matrix.nex": b"#NEXUS\nbegin data; dimensions ntax=1 nchar=2; format datatype=dna; end;\n",
    "matrix_same_line.nex": b"#NEXUS\nbegin data; dimensions ntax=1 nchar=2; format datatype=dna; matrix a AC\n; end;\n",
    "standard.nex": b"#NEXUS\nbegin data; dimensions ntax=2 nchar=3; format datatype=standard symbols=\"01\"; matrix\na 010\nb 110\n; end;\n",
}

# Inputs where the Rust port intentionally errors although Python succeeds
# (PLAN.md section 6: unsupported NEXUS features / non-ASCII residues).
RUST_ONLY_ERRORS = set()


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    for name, data in CASES.items():
        (out / name).write_bytes(data)
    (out / "RUST_ONLY_ERRORS").write_text("".join(f"{n}\n" for n in sorted(RUST_ONLY_ERRORS)))


if __name__ == "__main__":
    main()
