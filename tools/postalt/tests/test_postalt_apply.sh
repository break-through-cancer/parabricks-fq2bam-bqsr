#!/usr/bin/env bash
# Runs postalt_apply on the synthetic fixtures and compares with the checked-in expected output.
# Needs samtools and postalt_rs; run it inside the image: `make test`.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
cd "$W"
samtools view -b -o reads.bam "$HERE/fixtures/reads.sam"
samtools index reads.bam
postalt_apply --in reads.bam --ref "$HERE/fixtures/ref.fa" --alt "$HERE/fixtures/ref.alt" --out out --fmt bam --threads 2 --mem-gb 2
samtools quickcheck out.bam
samtools view out.bam | sort > got.sam
fail=0
if ! diff -u "$HERE/fixtures/expected.sam" got.sam; then echo "FAIL: output differs from expected.sam" >&2; fail=1; fi
[ "$(cat out.postalt_added.txt)" = "$(cat "$HERE/fixtures/expected_added.txt")" ] || { echo "FAIL: added record count" >&2; fail=1; }
samtools view -H out.bam | grep -q $'^@CO\tID:postalt' || { echo "FAIL: no @CO provenance line" >&2; fail=1; }
[ "$(samtools view -H out.bam | grep -c $'^@CO\tID:postalt')" = 1 ] || { echo "FAIL: @CO line is not unique" >&2; fail=1; }
[ "$(samtools view out.bam | grep -c 'PG:Z:MarkDuplicates')" = "$(samtools view reads.bam | grep -c 'PG:Z:MarkDuplicates')" ] || { echo "FAIL: PG tags were dropped" >&2; fail=1; }
# Generated records on the opposite strand must carry the reverse-complemented SEQ and reversed QUAL.
# (bwa-postalt.js under k8 >= 1.0 leaves them unchanged.)
while IFS=$'\t' read -r name seq qual; do
  got=$(samtools view out.bam chrT2_alt | awk -F'\t' -v n="$name" '$1 == n && int($2 / 2048) % 2 == 1 { print $10 "\t" $11 }')
  [ "$got" = "$(printf '%s\t%s' "$seq" "$qual")" ] || { echo "FAIL: $name generated record has wrong SEQ/QUAL: $got" >&2; fail=1; }
done < "$HERE/fixtures/expected_flipped.tsv"
# CRAM in and out: same records once CRAM's tag reordering and generated MD tags are ignored.
norm() { awk -F'\t' -v OFS='\t' '{n=0; delete t; for(i=12;i<=NF;i++) if($i !~ /^MD:Z:/) t[++n]=$i; for(a=2;a<=n;a++){v=t[a]; b=a-1; while(b>=1 && t[b]>v){t[b+1]=t[b]; b--} t[b+1]=v} s=""; for(i=1;i<=n;i++) s=s" "t[i]; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,s}' | sort; }
samtools view -C -T "$HERE/fixtures/ref.fa" -o reads.cram reads.bam
postalt_apply --in reads.cram --ref "$HERE/fixtures/ref.fa" --alt "$HERE/fixtures/ref.alt" --out outc --fmt cram --cram-version 3.0 --threads 2 --mem-gb 2
samtools quickcheck outc.cram
[ -s outc.cram.crai ] || { echo "FAIL: no CRAM index" >&2; fail=1; }
if ! diff <(norm < "$HERE/fixtures/expected.sam") <(samtools view -T "$HERE/fixtures/ref.fa" outc.cram | norm); then echo "FAIL: CRAM output differs from expected.sam" >&2; fail=1; fi
# An .alt contig class the selector does not know must stop the run.
printf 'chrT_hap1\t0\tchrT\t1\t60\t10M\t*\t0\t0\t*\t*\n' > bad.alt
if postalt_apply --in reads.bam --ref "$HERE/fixtures/ref.fa" --alt bad.alt --out bad --threads 2 --mem-gb 2 2> bad.err; then echo "FAIL: unknown ALT class accepted" >&2; fail=1; fi
grep -q 'outside the _alt, HLA- and _decoy name classes' bad.err || { echo "FAIL: unknown-class message missing" >&2; fail=1; }
[ "$fail" = 0 ] && echo "postalt_apply tests passed"
exit "$fail"
