#!/usr/bin/env python3

import re
from typing import Iterable

import pandas as pd

ALIGNMENT = re.compile(r"(?:^|/)(?P<name>[^/]+?)\.(?:md\.|recal\.)?(?P<ext>cram|bam)$")
INDEX_SUFFIX = {"bam": "bai", "cram": "crai"}


def build_samplesheet(paths: Iterable[str], samples: Iterable[str] = ()) -> pd.DataFrame:
    paths = list(paths)
    wanted = set(samples)
    present = set(paths)
    rows = []
    for path in sorted(paths):
        match = ALIGNMENT.search(path)
        if not match:
            continue
        if wanted and match["name"] not in wanted:
            continue
        index = f"{path}.{INDEX_SUFFIX[match['ext']]}"
        if index not in present:
            raise ValueError(f"No index found for {path}; expected {index}")
        rows.append({"sample": match["name"], "alignment": path, "alignment_index": index})
    if not rows:
        raise ValueError(f"No .cram or .bam files in the input dataset{' for the requested samples' if wanted else ''}.")
    return pd.DataFrame(rows, columns=["sample", "alignment", "alignment_index"])


def sources(ds, wanted):
    explicit = str(ds.params.get("alignment") or "").strip()
    if not explicit:
        return list(ds.files["file"]), wanted
    suffix = INDEX_SUFFIX[explicit.rsplit(".", 1)[-1]]
    return [explicit, f"{explicit}.{suffix}"], []


def main():
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()
    ds.logger.info(f"dataset files: {ds.files.to_dict(orient='records')}")
    wanted = [x.strip() for x in str(ds.params.get("samples") or "").split(",") if x.strip()]
    samplesheet = build_samplesheet(*sources(ds, wanted))
    for key in ("samples", "alignment"):
        ds.remove_param(key, force=True)
    ds.logger.info(samplesheet.to_csv(index=False))
    samplesheet.to_csv(ds.params["input"], index=False)


if __name__ == "__main__":
    main()
