"""Compare postalt outputs record by record.

usage: compare_records.py INPUT.sam K8_OUT.sam RS_DEFAULT_OUT.sam
Prints one JSON object: how the corrected rs output differs from k8, and whether each engine's
generated ALT records carry a SEQ/QUAL that matches their strand flag.
"""
import json
import sys

COMP = str.maketrans("WSATUGCYRKMBDHVNwsatugcyrkmbdhvn", "WSTAACGRYMKVHDBNwstaacgrymkvhdbn")


def revcomp(s):
    return s.translate(COMP)[::-1]


def body(path):
    with open(path) as f:
        for line in f:
            if line[0] != "@":
                yield line.rstrip("\n")


inp = set()
for line in body(sys.argv[1]):
    t = line.split("\t", 4)
    inp.add((t[0], t[1], t[2], t[3]))

stat = dict(records=0, identical=0, seq_qual_only=0, seq_qual_reversed_pairs=0, other_diff=0,
            generated=0, generated_flipped=0, k8_flipped_ok=0, rs_flipped_ok=0, k8_same_ok=0, rs_same_ok=0,
            groups_without_source=0)
examples = []
group_key, group_orig, group_gen = None, [], []


def flush():
    global group_orig, group_gen
    for engine, recs in group_gen:
        pass
    group_orig, group_gen = [], []


def check(gen, orig, engine, stat):
    """gen/orig are field lists; correct when SEQ follows the strand flag relative to some original record."""
    g_rev = int(gen[1]) & 16
    for o in orig:
        same = (int(o[1]) & 16) == g_rev
        want_seq = o[9] if same else revcomp(o[9])
        want_qual = o[10] if same else o[10][::-1]
        if gen[9] == want_seq and gen[10] == want_qual:
            return True, same
    return False, None


# Walk the two outputs together; collect groups by (qname, read end) for the ground-truth check.
k8_groups, rs_groups = {}, {}
with open(sys.argv[2]) as a, open(sys.argv[3]) as b:
    for la, lb in zip(a, b):
        if la[0] == "@" or lb[0] == "@":
            if la != lb:
                stat["other_diff"] += 1
            continue
        la, lb = la.rstrip("\n"), lb.rstrip("\n")
        stat["records"] += 1
        fa, fb = la.split("\t"), lb.split("\t")
        key = (fa[0], int(fa[1]) & 0xc0)
        for fields, groups in ((fa, k8_groups), (fb, rs_groups)):
            groups.setdefault(key, []).append(fields)
        if la == lb:
            stat["identical"] += 1
        elif fa[:9] == fb[:9] and fa[11:] == fb[11:]:
            stat["seq_qual_only"] += 1
            if fb[9] == revcomp(fa[9]) and fb[10] == fa[10][::-1]:
                stat["seq_qual_reversed_pairs"] += 1
            elif len(examples) < 3:
                examples.append({"k8": la[:300], "rs": lb[:300]})
        else:
            stat["other_diff"] += 1
            if len(examples) < 3:
                examples.append({"k8": la[:300], "rs": lb[:300]})
    if a.readline() or b.readline():
        stat["other_diff"] += 1

for key, k8g in k8_groups.items():
    rsg = rs_groups.get(key, [])
    orig_k8 = [f for f in k8g if (f[0], f[1], f[2], f[3]) in inp]
    gen_k8 = [f for f in k8g if (f[0], f[1], f[2], f[3]) not in inp]
    gen_rs = [f for f in rsg if (f[0], f[1], f[2], f[3]) not in inp]
    if not gen_k8:
        continue
    if not orig_k8:
        stat["groups_without_source"] += 1
        continue
    for gk, gr in zip(gen_k8, gen_rs):
        stat["generated"] += 1
        ok_k, same = check(gk, orig_k8, "k8", stat)
        ok_r, same_r = check(gr, orig_k8, "rs", stat)
        flipped = False
        for o in orig_k8:
            if (int(o[1]) & 16) != (int(gk[1]) & 16):
                flipped = True
        # Classify by the first original record's strand relation (the reported hit is the XA-bearing primary).
        rel_same = (int(orig_k8[0][1]) & 16) == (int(gk[1]) & 16)
        if not rel_same:
            stat["generated_flipped"] += 1
            stat["k8_flipped_ok"] += ok_k
            stat["rs_flipped_ok"] += ok_r
        else:
            stat["k8_same_ok"] += ok_k
            stat["rs_same_ok"] += ok_r
stat["examples"] = examples
json.dump(stat, sys.stdout)
print()
