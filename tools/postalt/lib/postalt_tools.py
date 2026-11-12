#!/usr/bin/env python3
"""Helpers for postalt_apply: ALT-file checks, name-hash partitioning, postalt ordering, and statistics."""
import json
import re
import sys
import zlib

ALT_CLASSES = re.compile(r"(_alt$|^HLA-|_decoy$)")


def check_alt(alt_path, fai_path):
    fai = {line.split("\t")[0] for line in open(fai_path)}
    names = set()
    for line in open(alt_path):
        if line[0] != "@" and line.count("\t") >= 10:
            names.add(line.split("\t")[0])
    if not names:
        sys.exit(f"postalt: {alt_path} lists no contigs")
    unknown = sorted(n for n in names if not ALT_CLASSES.search(n))
    if unknown:
        sys.exit(f"postalt: {alt_path} lists contigs outside the _alt, HLA- and _decoy name classes "
                 f"({len(unknown)}, e.g. {', '.join(unknown[:3])}); the selector cannot recognize them")
    missing = sorted(n for n in names if n not in fai)
    if missing:
        sys.exit(f"postalt: {len(missing)} contigs in {alt_path} are absent from the FASTA index (e.g. {', '.join(missing[:3])})")


def partition(k, prefix):
    outs = [open(f"{prefix}.{i}.sam", "w") for i in range(k)]
    header = []
    for line in sys.stdin:
        if line[0] == "@":
            header.append(line)
            continue
        outs[zlib.crc32(line[: line.index("\t")].encode()) % k].write(line)
    with open(f"{prefix}.header.sam", "w") as h:
        h.writelines(header)
    for o in outs:
        o.close()


def order():
    rows = []
    for n, line in enumerate(sys.stdin):
        t = line.split("\t", 3)
        flag = int(t[1])
        rows.append((t[0], 0 if flag & 64 else 1, 1 if flag & 2048 else 0, n, line))
    rows.sort(key=lambda r: r[:4])
    sys.stdout.writelines(r[4] for r in rows)


def stats(selected_names, selected_in, paths):
    """Sums the per-partition counters that postalt_apply writes: output, om-tagged, decreased, increased, MAPQ0 before, MAPQ0 after."""
    keys = ["output_records", "with_om_tag", "mapq_decreased", "mapq_increased", "mapq0_before", "mapq0_after"]
    totals = [0] * len(keys)
    for path in paths:
        for i, v in enumerate(open(path).read().split()):
            totals[i] += int(v)
    out = {"selected_names": selected_names, "selected_records": selected_in, **dict(zip(keys, totals))}
    out["added_records"] = out["output_records"] - selected_in
    json.dump(out, sys.stdout, indent=2)
    print()


if __name__ == "__main__":
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "check-alt":
        check_alt(*args)
    elif cmd == "partition":
        partition(int(args[0]), args[1])
    elif cmd == "order":
        order()
    elif cmd == "stats":
        stats(int(args[0]), int(args[1]), args[2:])
    else:
        sys.exit(f"unknown command {cmd}")
