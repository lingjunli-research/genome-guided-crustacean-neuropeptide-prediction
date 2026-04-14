# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# save as find_crz_candidates.py
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import regex as re

# Mature peptide core (no pQ or -NH2)
CORE = r"(?:Q|E)TFQYSRGWTN"
TAIL = r"G(?:K|R){1,2}"
PATTERN = re.compile(rf"(?b)({CORE}){{e<=1}}{TAIL}")  # allow <=1 edit in CORE

def six_frames(nt):
    s = str(nt).upper()
    fwd = [Seq(s).translate(), Seq(s[1:]).translate(), Seq(s[2:]).translate()]
    rc = str(Seq(s).reverse_complement())
    rev = [Seq(rc).translate(), Seq(rc[1:]).translate(), Seq(rc[2:]).translate()]
    return [str(x) for x in fwd + rev]

print("seq_id\tframe(0..5)\taa_start\taa_end\tcontext_±25aa\tmatched_pep")
for rec in SeqIO.parse(sys.argv[1], "fasta"):
    frames = six_frames(rec.seq)
    for fi, aa in enumerate(frames):
        for m in PATTERN.finditer(aa):
            a0, a1 = m.span()
            ctx = aa[max(0, a0-25):min(len(aa), a1+25)]
            print(f"{rec.id}\t{fi}\t{a0}\t{a1}\t{ctx}\t{m.group(0)}")


#python find_crz_candidates.py JC-genome_v3-230128\\JC-genome_v3-230128.fasta > crz_hits.tsv