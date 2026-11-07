#!/usr/bin/env bash
# usage (inside postalt-compare image): compare_engines.sh IN.cram REF.fa ALT OUTDIR [THREADS] [RECORDS_PER_PARTITION] [PARALLEL_JOBS]
set -euo pipefail
IN=$1; REF=$2; ALT=$3; OUT=$4; T=${5:-8}; PER=${6:-400000}; PJ=${7:-3}
HERE=${COMPARE_HOME:-/opt/postalt}
TOOLS=/opt/postalt/postalt_tools.py; JS=/opt/postalt/bwa-postalt.js
mkdir -p "$OUT"; W=$OUT/work; rm -rf "$W"; mkdir -p "$W"
log() { echo "[$(date +%T)] $*" >&2; }

log "selecting names"
s=$(date +%s)
samtools view --no-PG -@ "$T" -T "$REF" -e '(rname =~ "_alt$|^HLA-|_decoy$") || ([XA] =~ "_alt,|_decoy,|(^|;)HLA-")' "$IN" | cut -f1 | LC_ALL=C sort -u -S 2G -T "$W" > "$W/names.txt"
echo "select_names_s $(( $(date +%s) - s ))" >> "$OUT/stage_times.txt"
log "extracting selected records"
s=$(date +%s)
samtools view --no-PG -@ "$T" -T "$REF" -N "$W/names.txt" --output-fmt bam,level=1 -o "$W/sel.bam" "$IN"
echo "extract_s $(( $(date +%s) - s ))" >> "$OUT/stage_times.txt"
N=$(samtools view -c -@ "$T" "$W/sel.bam"); K=$(( (N + PER - 1) / PER ))
echo "selected_records $N names $(wc -l < "$W/names.txt") partitions $K" | tee -a "$OUT/stage_times.txt" >&2
samtools view --no-PG -@ "$T" -h "$W/sel.bam" | python3 -I "$TOOLS" partition "$K" "$W/part"; rm -f "$W/sel.bam"

job() {
  i=$1
  { cat "$W/part.header.sam"; python3 -I "$TOOLS" order < "$W/part.$i.sam"; } > "$W/in.$i.sam"; rm -f "$W/part.$i.sam"
  /usr/bin/time -f "%e %U %M" -o "$W/k8.$i.time" k8 "$JS" "$ALT" "$W/in.$i.sam" > "$W/k8.$i.sam" 2> /dev/null
  /usr/bin/time -f "%e %U %M" -o "$W/rc.$i.time" postalt_rs --k8-1.2-compat "$ALT" "$W/in.$i.sam" > "$W/rc.$i.sam"
  /usr/bin/time -f "%e %U %M" -o "$W/rd.$i.time" postalt_rs "$ALT" "$W/in.$i.sam" > "$W/rd.$i.sam"
  if cmp -s "$W/k8.$i.sam" "$W/rc.$i.sam"; then compat=identical; else compat=DIFFERENT; fi
  python3 -I "$HERE/compare_records.py" "$W/in.$i.sam" "$W/k8.$i.sam" "$W/rd.$i.sam" > "$W/cmp.$i.json"
  printf '%s\t%s\t%s\t%s\t%s\n' "$i" "$compat" "$(cat "$W/k8.$i.time")" "$(cat "$W/rc.$i.time")" "$(cat "$W/rd.$i.time")" > "$W/row.$i.tsv"
  rm -f "$W/k8.$i.sam" "$W/rc.$i.sam" "$W/rd.$i.sam"
  [ "$i" -lt 6 ] || rm -f "$W/in.$i.sam"
}
export -f job; export W TOOLS JS ALT HERE
log "comparing $K partitions ($PJ at a time)"
seq 0 $((K - 1)) | xargs -P "$PJ" -I{} bash -c 'job {}'

log "parallel scaling on the first partitions"
S=$(( K < 6 ? K : 6 ))
{ cat "$W/part.header.sam"; for i in $(seq 0 $((S - 1))); do cat "$W/in.$i.sam" | grep -v '^@'; done; } > "$W/scale.sam"
start=$(date +%s.%N); /usr/bin/time -f "%e %U %M" -o "$W/sc_k8.time" bash -c "for i in \$(seq 0 $((S - 1))); do k8 $JS $ALT $W/in.\$i.sam > /dev/null & done; wait"
for j in 1 $S; do /usr/bin/time -f "%e %U %M" -o "$W/sc_rs$j.time" postalt_rs -j "$j" "$ALT" "$W/scale.sam" > /dev/null; done
rm -f "$W"/in.*.sam "$W/scale.sam"
python3 -I - "$W" "$OUT" "$K" "$S" <<'PY'
import glob, json, sys
w, out, k, s = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
rows = [open(p).read().rstrip("\n").split("\t") for p in sorted(glob.glob(f"{w}/row.*.tsv"))]
agg = {"partitions": k, "compat_identical": sum(r[1] == "identical" for r in rows), "compat_different": sum(r[1] != "identical" for r in rows)}
for name, col in (("k8", 2), ("rs_compat", 3), ("rs_default", 4)):
    t = [r[col].split() for r in rows]
    agg[f"{name}_single_thread"] = {"wall_s": round(sum(float(x[0]) for x in t), 1), "cpu_s": round(sum(float(x[1]) for x in t), 1), "max_rss_mb": round(max(int(x[2]) for x in t) / 1024, 1)}
cmp_keys = {}
examples = []
for p in sorted(glob.glob(f"{w}/cmp.*.json")):
    d = json.load(open(p)); examples += d.pop("examples", [])
    for key, v in d.items(): cmp_keys[key] = cmp_keys.get(key, 0) + v
agg["default_vs_k8"] = cmp_keys
agg["examples"] = examples[:5]
sc = {}
for name, f in (("k8_parallel", "sc_k8"), (f"rs_j1", "sc_rs1"), (f"rs_j{s}", f"sc_rs{s}")):
    x = open(f"{w}/{f}.time").read().split(); sc[name] = {"wall_s": float(x[0]), "cpu_s": float(x[1]), "max_rss_mb": round(int(x[2]) / 1024, 1)}
agg["scaling_first_partitions"] = sc
json.dump(agg, open(f"{out}/summary.json", "w"), indent=2)
print(json.dumps(agg, indent=2))
PY
