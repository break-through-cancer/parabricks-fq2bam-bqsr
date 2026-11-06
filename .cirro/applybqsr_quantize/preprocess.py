#!/usr/bin/env python3

import re
from collections import defaultdict
from typing import Iterable

import pandas as pd

PARABRICKS_FILE = re.compile(r"(?:^|/)preprocessing/parabricks/(?P<sample>[^/]+)/(?P<name>[^/]+)$")
MARKDUPLICATES_FILE = re.compile(r"(?:^|/)preprocessing/markduplicates/(?P<sample>[^/]+)/(?P<name>[^/]+)$")
RECAL_TABLE_FILE = re.compile(r"(?:^|/)preprocessing/recal_table/(?P<sample>[^/]+)/(?P<name>[^/]+)$")
INDEX_SUFFIX = {"bam": "bai", "cram": "crai"}
COLUMNS = ["sample", "alignment", "alignment_index", "recal_table"]
RECALIBRATED_FILE = re.compile(r"(?:^|/)preprocessing/recalibrated/(?P<sample>[^/]+)/(?P<name>[^/]+)$")
QUANTIZE_SOURCES = {"recalibrated": RECALIBRATED_FILE, "markduplicates": MARKDUPLICATES_FILE}
QUANTIZE_COLUMNS = ["sample", "alignment", "alignment_index"]
MODES = ("apply_bqsr_and_quantize", "apply_bqsr", "quantize")
FORM_ONLY = ("genome_source", "genome_index", "mode", "alignment_source")


def build_samplesheet(paths: Iterable[str]) -> pd.DataFrame:
    """Build the pipeline samplesheet from an upstream dataset's file paths.

    Accepts either sarek_align's layout (preprocessing/parabricks/<sample>/: the
    pre-BQSR fq2bam alignment, its index and the recalibration table, co-located) or
    this pipeline's own apply_bqsr=false output (preprocessing/markduplicates/<sample>/:
    alignment and index; preprocessing/recal_table/<sample>/: the table, separately).
    Alignments from other stages (for example recalibrated/) are never used, since
    applying the table again would recalibrate those reads twice.
    """
    paths = list(paths)
    present = set(paths)
    found = defaultdict(lambda: {"alignments": [], "tables": []})

    for path in paths:
        match = PARABRICKS_FILE.search(path) or MARKDUPLICATES_FILE.search(path)
        if match:
            name = match["name"]
            ext = name.rsplit(".", 1)[-1]
            if ext in INDEX_SUFFIX:
                found[match["sample"]]["alignments"].append(path)
            elif ext == "table":
                found[match["sample"]]["tables"].append(path)
            continue
        match = RECAL_TABLE_FILE.search(path)
        if match and match["name"].endswith(".table"):
            found[match["sample"]]["tables"].append(path)

    if not found:
        raise ValueError(
            "No files under preprocessing/parabricks/<sample>/ (sarek_align, Parabricks "
            "aligner, known sites supplied, save_mapped on, baserecalibrator skipped) or "
            "preprocessing/markduplicates/<sample>/ + preprocessing/recal_table/<sample>/ "
            "(this pipeline's own apply_bqsr=false output) in the input dataset(s). This "
            "pipeline needs the pre-BQSR alignment, its index, and the recalibration table."
        )

    rows, problems = [], []
    for sample in sorted(found):
        alignments, tables = found[sample]["alignments"], found[sample]["tables"]
        if len(alignments) != 1:
            problems.append(f"{sample}: expected one BAM/CRAM, found {len(alignments)}")
        if len(tables) != 1:
            problems.append(
                f"{sample}: no .table recalibration table" if not tables
                else f"{sample}: expected one .table, found {len(tables)}"
            )
        if len(alignments) != 1 or len(tables) != 1:
            continue
        alignment = alignments[0]
        index = f"{alignment}.{INDEX_SUFFIX[alignment.rsplit('.', 1)[-1]]}"
        if index not in present:
            problems.append(f"{sample}: no index for {alignment.rsplit('/', 1)[-1]}")
            continue
        rows.append(dict(sample=sample, alignment=alignment, alignment_index=index, recal_table=tables[0]))

    if problems:
        raise ValueError("Cannot build the samplesheet:\n  " + "\n  ".join(problems))
    return pd.DataFrame(rows, columns=COLUMNS)


def build_quantize_samplesheet(paths: Iterable[str], source: str = "auto") -> pd.DataFrame:
    """Build the quantize-only samplesheet (no recalibration table) from an input dataset.

    source picks the folder to quantize: recalibrated (BQSR applied, not yet quantized),
    markduplicates (the raw aligned alignment, for reads already recalibrated in the original
    BAM, so BQSR is never applied twice), or auto, which takes recalibrated per sample when
    present and markduplicates otherwise. Nothing checks that an alignment was not already
    quantized.
    """
    if source != "auto" and source not in QUANTIZE_SOURCES:
        raise ValueError(f"alignment_source must be auto, recalibrated or markduplicates; got {source!r}")
    paths = list(paths)
    present = set(paths)
    found = {name: defaultdict(list) for name in QUANTIZE_SOURCES}
    for path in paths:
        for name, pattern in QUANTIZE_SOURCES.items():
            match = pattern.search(path)
            if match and match["name"].rsplit(".", 1)[-1] in INDEX_SUFFIX:
                found[name][match["sample"]].append(path)

    if not any(found.values()):
        raise ValueError(
            "No alignments under preprocessing/recalibrated/<sample>/ or "
            "preprocessing/markduplicates/<sample>/ in the input dataset(s). Quantize only "
            "needs an existing alignment and its index."
        )

    rows, problems = [], []
    order = list(QUANTIZE_SOURCES) if source == "auto" else [source]
    for sample in sorted({s for name in order for s in found[name]}):
        name = next(n for n in order if sample in found[n])
        alignments = found[name][sample]
        if len(alignments) != 1:
            problems.append(f"{sample}: expected one BAM/CRAM under {name}, found {len(alignments)}")
            continue
        alignment = alignments[0]
        index = f"{alignment}.{INDEX_SUFFIX[alignment.rsplit('.', 1)[-1]]}"
        if index not in present:
            problems.append(f"{sample}: no index for {alignment.rsplit('/', 1)[-1]}")
            continue
        rows.append(dict(sample=sample, alignment=alignment, alignment_index=index))

    if source != "auto" and not rows and not problems:
        problems.append(f"no alignments under preprocessing/{source}/")
    if problems:
        raise ValueError("Cannot build the samplesheet:\n  " + "\n  ".join(problems))
    return pd.DataFrame(rows, columns=QUANTIZE_COLUMNS)


def build_for_mode(paths: Iterable[str], mode: str, source: str = "auto") -> tuple:
    """Return (samplesheet, quantize_quals_enabled) for a form mode."""
    if mode not in MODES:
        raise ValueError(f"mode must be one of {', '.join(MODES)}; got {mode!r}")
    if mode == "quantize":
        return build_quantize_samplesheet(paths, source), True
    return build_samplesheet(paths), mode == "apply_bqsr_and_quantize"


def apply_genome_params(params: dict) -> dict:
    """Translate Cirro form fields into pipeline genome parameters."""
    p = dict(params)
    if p.get("genome_source", "igenomes") == "dataset":
        index = p.get("genome_index")
        if not index:
            raise ValueError("Custom genome selected but no BWA genome index dataset was provided")
        p["genome"] = "null"
        p["bwa_index"] = index
        p["ref_fasta"] = f"{index}/genome.fasta"
        p["ref_fasta_fai"] = f"{index}/genome.fasta.fai"
    for k in FORM_ONLY:
        p.pop(k, None)
    return p


if __name__ == "__main__":
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()
    mode = ds.params.get("mode", "apply_bqsr_and_quantize")
    samplesheet, quantize = build_for_mode(ds.files["file"], mode, ds.params.get("alignment_source", "auto"))
    ds.add_param("quantize_quals_enabled", quantize, overwrite=True)
    ds.logger.info(samplesheet.to_csv(index=False))
    samplesheet.to_csv(ds.params["input"], index=False)
    ds.logger.info(f"Wrote {samplesheet.shape[0]} row(s) to {ds.params['input']}")
    resolved = apply_genome_params(ds.params)
    for key in FORM_ONLY:
        ds.remove_param(key, force=True)
    for key in ("genome", "bwa_index", "ref_fasta", "ref_fasta_fai"):
        if key in resolved:
            ds.add_param(key, resolved[key], overwrite=True)
