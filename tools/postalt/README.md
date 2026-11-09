# postalt_apply

Runs `postalt_rs` (a Rust port of the BWA-kit alt-aware post-processing script `bwa-postalt.js`, r985) on a coordinate-sorted
BAM/CRAM without name-sorting the whole file: it selects the reads the script can change, runs the
script on those, and merges the result back.

```
postalt_apply --in FILE --ref FASTA --alt ALT --out PREFIX [--fmt bam|cram] [--cram-version V] [--threads N] [--mem-gb G]
```

Outputs: `PREFIX.<bam|cram>` and its index, `PREFIX.postalt.json` (counts), `PREFIX.postalt_added.txt`
(records added), `PREFIX.postalt_mqc.tsv` (MultiQC table).

## Selection

`postalt_rs`, like `bwa-postalt.js`, treats every contig named in the `.alt` file as ALT: `*_alt`, `HLA-*` and `*_decoy`.
It changes a read only when the read's group (name and read end) has a record on such a contig or an
`XA` hit on one. That test is the selector, applied to every name in one pass; all records of a
selected name go through the script, everything else passes through untouched. The tool stops if the
`.alt` file lists a contig outside those three name classes.

Selecting by genomic position instead is not exact: the script lowers the MAPQ of reads with a
decoy supplementary alignment and `pa < 1` far from any ALT locus.

## Merge details

- `samtools merge -c -p` keeps read-group and program IDs. Parabricks leaves `PG:Z:MarkDuplicates`
  tags on reads with no `@PG` line, and merge drops such tags, so the header gets a minimal `@PG`
  line for each missing ID found in the first 200,000 records.
- A single `@CO ID:postalt` line records the script version and the `.alt` checksum.
- CRAM re-encoding reorders tags and adds `MD` tags to the new ALT records; the content is the same.

## Pinned versions

samtools and htslib 1.22.1, checked against SHA-256 values in the Dockerfile, and `postalt_rs` built
locked from bsbludwig/postalt_rs commit `b48e8b3` (GPL-3.0). The image is built for linux/amd64.

## Licensing

`postalt_rs` is GPL-3.0. The image contains its binary, so publishing the image distributes that
program. The corresponding source is the pinned commit above, and the image carries the license and
notice files in `/usr/share/doc/postalt_rs/` with a source link in the
`io.github.break-through-cancer.postalt_rs.source` label. samtools and htslib are MIT-licensed, with
their license files in `/usr/local/share/doc/`.

The pipeline runs `postalt_rs` as a separate process through standard input, standard output, and
command-line arguments. `postalt_apply` and `lib/postalt_tools.py` call it that way and do not link
against it, so the GPL covers `postalt_rs` and not the pipeline code.

## Why not the original script

`bwa-postalt.js` run under k8 1.x does not reverse-complement SEQ or reverse QUAL on the opposite-strand
records it generates (indexed Bytes access was removed in k8 1.0). `postalt_rs` corrects this by default;
its `--k8-1.2-compat` mode reproduces the original output byte for byte, and was identical on every
partition of two whole-genome samples. The fixture test checks the corrected output.

## Tests

`make test` builds the image and runs `tests/test_postalt_apply.sh` inside it against synthetic
fixtures (`tests/make_fixtures.py` regenerates them): a read whose MAPQ drops through a decoy
supplementary, one that rises on an ALT contig, an ALT-versus-repeat tie, a pair that gains a
supplementary record, untouched reads, CRAM in and out, and the unknown-contig-class stop.
`tests/fixtures/expected.sam` is the output of the reference script on those fixtures.
