#!/usr/bin/env python3

import pandas as pd

INDEX_SUFFIX = {"bam": "bai", "cram": "crai"}
COLUMNS = ["sample", "alignment", "alignment_index", "recal_table"]
QUANTIZE_COLUMNS = ["sample", "alignment", "alignment_index"]
MODES = ("apply_bqsr_and_quantize", "apply_bqsr", "quantize")
FORM_ONLY = ("genome_source", "genome_index", "mode")


def build_samplesheet(table: pd.DataFrame, mode: str, log=None) -> pd.DataFrame:
    """Turn the table filled in on the form into the pipeline samplesheet. A row is one sample:
    its alignment (BAM or CRAM), the index (derived from the alignment when blank) and, for the
    modes that apply BQSR, the recalibration table. In quantize mode the table column is ignored.
    Under the apply modes an alignment from a recalibrated/ folder is rejected, since applying
    the table again would recalibrate those reads twice."""
    if mode not in MODES:
        raise ValueError(f"mode must be one of {', '.join(MODES)}; got {mode!r}")
    if table.empty:
        raise ValueError("The samplesheet has no rows")
    missing_columns = [c for c in ("sample", "alignment") if c not in table.columns]
    if missing_columns:
        raise ValueError("The samplesheet has no column " + ", ".join(missing_columns))
    t = table.reindex(columns=COLUMNS).fillna("").astype(str)
    t = t.apply(lambda column: column.str.strip())
    apply_bqsr = mode != "quantize"

    problems = []
    for column in ("sample", "alignment") + (("recal_table",) if apply_bqsr else ()):
        blank = t.index[t[column] == ""] + 1
        if len(blank):
            problems.append(f"{column} is blank in row(s) " + ", ".join(map(str, blank)))
    for sample, n in t.loc[t["sample"] != "", "sample"].value_counts().items():
        if n > 1:
            problems.append(f"{sample} is listed {n} times; give each sample one row")

    indexes = []
    for row_number, row in t.iterrows():
        alignment, index, label = row["alignment"], row["alignment_index"], f"{row['sample'] or 'row ' + str(row_number + 1)}"
        extension = alignment.rsplit(".", 1)[-1]
        if alignment == "":
            indexes.append("")
            continue
        if extension not in INDEX_SUFFIX:
            problems.append(f"{label}: alignment must be a .bam or .cram file")
            indexes.append("")
            continue
        if apply_bqsr and "recalibrated" in alignment.split("/")[:-1]:
            problems.append(f"{label}: {alignment.rsplit('/', 1)[-1]} is in a recalibrated folder; "
                            "applying the table again would recalibrate it twice. Choose the duplicate-marked alignment")
        if index == "":
            index = f"{alignment}.{INDEX_SUFFIX[extension]}"
        elif not index.endswith("." + INDEX_SUFFIX[extension]):
            problems.append(f"{label}: index for a .{extension} must end in .{INDEX_SUFFIX[extension]}")
        indexes.append(index)
    t["alignment_index"] = indexes
    if apply_bqsr:
        for sample, table_path in zip(t["sample"], t["recal_table"]):
            if table_path and not table_path.endswith(".table"):
                problems.append(f"{sample}: recalibration table must be a .table file")
    elif log is not None and (t["recal_table"] != "").any():
        log.warning("Quantize only: the recal_table column is ignored")

    files = pd.concat([t["alignment"], t["alignment_index"]] + ([t["recal_table"]] if apply_bqsr else []))
    for path, n in files[files != ""].value_counts().items():
        if n > 1:
            problems.append(f"{path} is listed {n} times")
    if problems:
        raise ValueError("Cannot build the samplesheet:\n  " + "\n  ".join(problems))
    return t[COLUMNS if apply_bqsr else QUANTIZE_COLUMNS].sort_values("sample").reset_index(drop=True)


def build_for_mode(table: pd.DataFrame, mode: str, log=None) -> tuple:
    """Return (samplesheet, quantize_quals_enabled) for a form mode."""
    return build_samplesheet(table, mode, log), mode != "apply_bqsr"


def apply_genome_params(params: dict) -> dict:
    """Translate Cirro form fields into pipeline genome parameters."""
    p = dict(params)
    if p.get("genome_source", "igenomes") == "dataset":
        index = p.get("genome_index")
        if not index:
            raise ValueError("Custom genome selected but no BWA genome index dataset was provided")
        index = index.rstrip("/")
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
    table = pd.read_csv(ds.params["input"], dtype=str, keep_default_na=False)
    samplesheet, quantize = build_for_mode(table, mode, ds.logger)
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
