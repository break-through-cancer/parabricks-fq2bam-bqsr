#!/usr/bin/env python3

import pandas as pd

COLUMNS = ["patient", "sample", "status", "lane", "fastq_1", "fastq_2"]
FORM_ONLY = ("genome_source", "genome_index", "dbsnp", "known_indels", "custom_intervals", "use_intervals", "intervals_mode")


def build_fastq_samplesheet(table: pd.DataFrame, log) -> pd.DataFrame:
    """Turn the table filled in on the form into the pipeline samplesheet. A row is one lane: a FASTQ
    pair, or one FASTQ with fastq_2 blank for a single-end (singleton) lane. Lanes are numbered per
    sample, paired lanes first, so a sample's read-group IDs do not depend on the other samples in the run."""
    if table.empty:
        raise ValueError("The samplesheet has no rows")
    missing_columns = [c for c in ("patient", "sample", "fastq_1") if c not in table.columns]
    if missing_columns:
        raise ValueError("The samplesheet has no column " + ", ".join(missing_columns))
    t = table.reindex(columns=["patient", "sample", "status", "fastq_1", "fastq_2"]).fillna("").astype(str)
    t = t.apply(lambda column: column.str.strip())

    problems = []
    for column in ("patient", "sample", "fastq_1"):
        blank = t.index[t[column] == ""] + 1
        if len(blank):
            problems.append(f"{column} is blank in row(s) " + ", ".join(map(str, blank)))
    bad_status = t.loc[~t["status"].isin(["", "0", "1"]), "status"].unique()
    if len(bad_status):
        problems.append("status must be 0 (normal) or 1 (tumor), got " + ", ".join(bad_status))
    files = pd.concat([t["fastq_1"], t["fastq_2"]])
    for path, n in files[files != ""].value_counts().items():
        if n > 1:
            problems.append(f"{path} is listed {n} times")
    for sample, group in t.groupby("sample", sort=True):
        for column in ("patient", "status"):
            values = sorted(set(group[column]) - {""}) if column == "status" else sorted(set(group[column]))
            if len(values) > 1:
                problems.append(f"{sample}: more than one {column} ({', '.join(values) or 'blank'})")
    if problems:
        raise ValueError("Cannot build the samplesheet:\n  " + "\n  ".join(problems))

    unset = int((t.groupby("sample")["status"].apply(lambda s: (s == "").all())).sum())
    if unset:
        log.warning(f"status not provided for {unset} sample(s), defaulting to 0 (normal)")
    # No fallback from patient to the sample id: it feeds SM:${patient}_${sample} in the CRAM's read
    # group (subworkflows/local/samplesheet_to_fastqs/main.nf), and a silent patient==sample fallback
    # produces a self-doubled SM that downstream consumers (e.g. ICONICC-NF's SUBSET_PAIR_VCF) cannot anticipate.
    status = t.groupby("sample")["status"].transform(lambda s: next((v for v in s if v != ""), "0")).astype(int)
    out = t.assign(status=status, single=t["fastq_2"] == "").sort_values(["sample", "single", "fastq_1"])
    out["lane"] = out.groupby("sample").cumcount().astype(str)
    return out.reset_index(drop=True)[COLUMNS]


def apply_genome_params(params: dict) -> dict:
    """Translate Cirro form fields into pipeline genome parameters."""
    p = dict(params)
    source = p.get("genome_source", "igenomes")
    if source == "dataset":
        index = p.get("genome_index")
        if not index:
            raise ValueError("Custom genome selected but no BWA genome index dataset was provided")
        p["genome"] = "null"
        # The folder picker stores a path that ends in a slash; joining onto it would give "//".
        index = index.rstrip("/")
        p["bwa_index"] = index
        p["ref_fasta"] = f"{index}/genome.fasta"
        p["ref_fasta_fai"] = f"{index}/genome.fasta.fai"
        sites = [p.get(k) for k in ("dbsnp", "known_indels") if p.get(k)]
        if sites:
            p["known_sites"] = ",".join(sites)
        if p.get("custom_intervals"):
            p["intervals"] = p["custom_intervals"]
    else:
        mode = p.get("intervals_mode") or ("none" if p.get("use_intervals") is False else "gatk_calling_regions")
        if mode == "none":
            p["no_intervals"] = True
        elif mode == "custom":
            if not p.get("custom_intervals"):
                raise ValueError("Intervals set to custom intervals but no intervals BED was selected")
            p["intervals"] = p["custom_intervals"]
    for k in FORM_ONLY:
        p.pop(k, None)
    return p


if __name__ == "__main__":
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()
    table = pd.read_csv(ds.params["input"], dtype=str, keep_default_na=False)
    sheet = build_fastq_samplesheet(table, ds.logger)
    sheet.to_csv(ds.params["input"], index=False)
    resolved = apply_genome_params(ds.params)
    for key in FORM_ONLY:
        ds.remove_param(key, force=True)
    for key in ("genome", "bwa_index", "ref_fasta", "ref_fasta_fai", "known_sites", "intervals", "no_intervals"):
        if key in resolved:
            ds.add_param(key, resolved[key], overwrite=True)
