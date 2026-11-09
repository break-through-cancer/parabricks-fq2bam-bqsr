#!/usr/bin/env python3
"""Writes the synthetic postalt fixtures: a three-contig reference, its .alt file, and reads covering each postalt path."""
import random

random.seed(7)
contigs = {"chrT": 700, "chrT1_alt": 200, "chrUn_x_decoy": 100, "chrT2_alt": 150}
seq = {c: "".join(random.choice("ACGT") for _ in range(n)) for c, n in {"chrT": 500, "chrT1_alt": 200, "chrUn_x_decoy": 100}.items()}
seq["chrT"] += "".join(random.choice("ACGT") for _ in range(200))
COMP = str.maketrans("ACGT", "TGCA")
revcomp = lambda s: s.translate(COMP)[::-1]
seq["chrT2_alt"] = revcomp(seq["chrT"][500:650])

with open("fixtures/ref.fa", "w") as f:
    for c, s in seq.items():
        f.write(f">{c}\n{s}\n")

header = ["@HD\tVN:1.6\tSO:coordinate"] + [f"@SQ\tSN:{c}\tLN:{n}" for c, n in contigs.items()] + ["@RG\tID:rg1\tSM:S1\tLB:L1\tPL:ILLUMINA"]
with open("fixtures/ref.alt", "w") as f:
    f.write("\n".join(header[1:5]) + "\n")
    f.write("chrT1_alt\t0\tchrT\t101\t60\t200M\t*\t0\t0\t*\t*\n")
    f.write("chrT2_alt\t16\tchrT\t501\t60\t150M\t*\t0\t0\t*\t*\n")
    f.write("chrUn_x_decoy\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n")

L = 50
def rec(name, flag, ctg, pos, mapq, tags, mate=("*", 0, 0), ref_ctg=None, qual=None):
    s = seq[ref_ctg or ctg][pos - 1:pos - 1 + L]
    return [name, flag, ctg, pos, mapq, f"{L}M", mate[0], mate[1], mate[2], s, qual or "I" * L] + ["RG:Z:rg1", "PG:Z:MarkDuplicates"] + tags

reads = [
    rec("pair_alt_hit", 99, "chrT", 150, 60, ["XA:Z:chrT1_alt,+50,50M,0;", "NM:i:0"], ("=", 400, 300)),
    rec("decoy_sa", 0, "chrT", 20, 60, ["pa:f:0.8", "SA:Z:chrUn_x_decoy,5,+,50M,60,0;", "NM:i:0"]),
    rec("repeat_tie", 0, "chrT", 200, 60, ["XA:Z:chrT1_alt,+100,50M,0;chrT,+450,50M,0;", "NM:i:0"]),
    rec("on_alt", 0, "chrT1_alt", 20, 30, ["XA:Z:chrT,+120,50M,0;", "NM:i:0"]),
    rec("plain_in_region", 0, "chrT", 250, 60, ["NM:i:0"]),
    rec("plain_outside", 0, "chrT", 330, 60, ["NM:i:0"]),
    rec("decoy_sa", 2048, "chrUn_x_decoy", 5, 60, ["NM:i:0"], ref_ctg="chrUn_x_decoy"),
    rec("pair_alt_hit", 147, "chrT", 400, 60, ["NM:i:0"], ("=", 150, -300)),
]
ramp = "".join(chr(33 + 10 + i) for i in range(L))
flip = [
    rec("flip_fwd", 0, "chrT", 520, 60, ["XA:Z:chrT2_alt,-82,50M,0;", "NM:i:0"], qual=ramp),
    rec("flip_rev", 16, "chrT", 580, 60, ["XA:Z:chrT2_alt,+22,50M,0;", "NM:i:0"], qual=ramp),
]
reads += flip
with open("fixtures/expected_flipped.tsv", "w") as f:
    for r in flip:
        f.write("\t".join([r[0], revcomp(r[9]), r[10][::-1]]) + "\n")
reads.sort(key=lambda r: (list(contigs).index(r[2]), r[3]))
with open("fixtures/reads.sam", "w") as f:
    f.write("\n".join(header) + "\n")
    for r in reads:
        f.write("\t".join(str(x) for x in r) + "\n")
