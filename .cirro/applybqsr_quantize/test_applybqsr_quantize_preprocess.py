import importlib.util
import json
import pathlib
import re

import pandas as pd
import pytest

HERE = pathlib.Path(__file__).parent
_spec = importlib.util.spec_from_file_location("apply_bqsr_preprocess", HERE / "preprocess.py")
_preprocess = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_preprocess)
build_samplesheet = _preprocess.build_samplesheet
build_for_mode = _preprocess.build_for_mode
apply_genome_params = _preprocess.apply_genome_params

ROOT = "s3://bucket/datasets/ds1/data/preprocessing"
APPLY = "apply_bqsr_and_quantize"


class Log:
    def __init__(self):
        self.warnings = []

    def warning(self, message):
        self.warnings.append(message)


def table(*rows):
    return pd.DataFrame(rows, columns=["sample", "alignment", "alignment_index", "recal_table"])


def row(sample="S1", folder="parabricks", ext="bam", index=None, recal=True):
    path = f"{ROOT}/{folder}/{sample}/{sample}.{ext}"
    suffix = {"bam": "bai", "cram": "crai"}[ext]
    return (sample, path, f"{path}.{suffix}" if index is None else index, f"{ROOT}/recal_table/{sample}/{sample}.table" if recal else "")


def test_one_sample_bam():
    sheet = build_samplesheet(table(row()), APPLY)
    assert list(sheet.columns) == ["sample", "alignment", "alignment_index", "recal_table"]
    assert sheet.to_dict("records") == [{
        "sample": "S1",
        "alignment": f"{ROOT}/parabricks/S1/S1.bam",
        "alignment_index": f"{ROOT}/parabricks/S1/S1.bam.bai",
        "recal_table": f"{ROOT}/recal_table/S1/S1.table",
    }]


def test_cram_and_multiple_samples_sorted():
    sheet = build_samplesheet(table(row("S2", ext="cram"), row("S1")), APPLY)
    assert list(sheet["sample"]) == ["S1", "S2"]
    assert sheet.loc[1, "alignment_index"].endswith("S2.cram.crai")


def test_blank_index_is_derived_from_the_alignment():
    sheet = build_samplesheet(table(row("S1", index=""), row("S2", ext="cram", index="")), APPLY)
    assert sheet.loc[0, "alignment_index"].endswith("S1.bam.bai")
    assert sheet.loc[1, "alignment_index"].endswith("S2.cram.crai")


def test_an_index_with_the_other_formats_extension_is_an_error():
    with pytest.raises(ValueError, match="S1: index for a .bam must end in .bai"):
        build_samplesheet(table(row(index=f"{ROOT}/parabricks/S1/S1.cram.crai")), APPLY)


def test_apply_modes_reject_a_recalibrated_alignment():
    # Applying the table again to an already-recalibrated alignment would recalibrate it twice.
    for mode in ("apply_bqsr_and_quantize", "apply_bqsr"):
        with pytest.raises(ValueError, match="recalibrate it twice"):
            build_samplesheet(table(row(folder="recalibrated")), mode)


def test_apply_modes_require_the_recalibration_table():
    with pytest.raises(ValueError, match="recal_table is blank in row\\(s\\) 1"):
        build_samplesheet(table(row(recal=False)), "apply_bqsr")


def test_table_must_be_a_table_file():
    bad = (row()[0], row()[1], row()[2], f"{ROOT}/recal_table/S1/S1.txt")
    with pytest.raises(ValueError, match="S1: recalibration table must be a .table file"):
        build_samplesheet(table(bad), APPLY)


def test_alignment_must_be_bam_or_cram():
    with pytest.raises(ValueError, match="S1: alignment must be a .bam or .cram file"):
        build_samplesheet(table(("S1", f"{ROOT}/S1.sam", "", f"{ROOT}/S1.table")), APPLY)


def test_blank_cells_and_duplicates_are_reported_together():
    with pytest.raises(ValueError) as err:
        build_samplesheet(table(row("S1"), row("S1"), ("", "", "", "")), APPLY)
    msg = str(err.value)
    assert "sample is blank in row(s) 3" in msg
    assert "alignment is blank in row(s) 3" in msg
    assert "S1 is listed 2 times" in msg
    assert f"{ROOT}/parabricks/S1/S1.bam is listed 2 times" in msg


def test_the_same_file_for_two_samples_is_an_error():
    a, b = row("S1"), row("S2")
    with pytest.raises(ValueError, match="S1.bam is listed 2 times"):
        build_samplesheet(table(a, (b[0], a[1], a[2], b[3])), APPLY)


def test_empty_table_and_missing_columns():
    with pytest.raises(ValueError, match="no rows"):
        build_samplesheet(table(), APPLY)
    with pytest.raises(ValueError, match="no column alignment"):
        build_samplesheet(pd.DataFrame({"sample": ["S1"]}), APPLY)


def test_quantize_needs_only_the_alignment_and_ignores_the_table():
    log = Log()
    quant = build_samplesheet(table(row(folder="recalibrated", recal=False)), "quantize", log)
    assert list(quant.columns) == ["sample", "alignment", "alignment_index"]
    assert quant.loc[0, "alignment"].endswith("recalibrated/S1/S1.bam") and not log.warnings
    build_samplesheet(table(row()), "quantize", log)
    assert log.warnings == ["Quantize only: the recal_table column is ignored"]


def test_mode_decides_quantization():
    bqsr = table(row())
    assert build_for_mode(bqsr, "apply_bqsr_and_quantize")[1] is True
    assert build_for_mode(bqsr, "apply_bqsr")[1] is False
    quant, flag = build_for_mode(table(row(recal=False)), "quantize")
    assert flag is True and list(quant.columns) == ["sample", "alignment", "alignment_index"]


def test_unknown_mode_is_an_error():
    with pytest.raises(ValueError, match="mode must be one of"):
        build_for_mode(table(row()), "everything")


def test_igenomes_is_the_default_genome():
    p = apply_genome_params({"genome": "GATK.GRCh38", "igenomes_base": "s3://refs/igenomes/"})
    assert p["genome"] == "GATK.GRCh38" and "ref_fasta" not in p
    p = apply_genome_params({"genome_source": "igenomes", "genome": "GATK.GRCh38"})
    assert p["genome"] == "GATK.GRCh38" and "genome_source" not in p


def test_custom_genome_resolves_dataset_files_and_strips_the_trailing_slash():
    for index in ("s3://b/ds/data", "s3://b/ds/data/"):
        p = apply_genome_params({"genome_source": "dataset", "genome_index": index})
        assert p["genome"] == "null"
        assert p["bwa_index"] == "s3://b/ds/data"
        assert p["ref_fasta"] == "s3://b/ds/data/genome.fasta" and p["ref_fasta_fai"] == "s3://b/ds/data/genome.fasta.fai"
        assert not {"genome_source", "genome_index"} & p.keys()


def test_custom_genome_without_index_is_an_error():
    with pytest.raises(ValueError, match="BWA genome index"):
        apply_genome_params({"genome_source": "dataset"})


def test_form_only_fields_are_dropped_from_the_pipeline_params():
    out = apply_genome_params({"mode": "quantize", "genome_source": "igenomes", "output_fmt": "bam"})
    assert out == {"output_fmt": "bam"}


def form():
    return json.loads((HERE / "process-form.json").read_text())


def test_bam_is_the_first_output_format_and_cram_version_shows_only_for_cram():
    group = form()["form"]["properties"]["output"]
    assert group["properties"]["output_fmt"]["enum"] == ["bam", "cram"]
    assert group["properties"]["output_fmt"]["default"] == "bam"
    assert "cram_version" not in group["properties"]
    branches = {b["properties"]["output_fmt"]["enum"][0]: b for b in group["dependencies"]["output_fmt"]["oneOf"]}
    assert "cram_version" not in branches["bam"]["properties"]
    assert branches["cram"]["properties"]["cram_version"]["default"] == "3.0"


def test_form_and_input_expose_the_genome_selection():
    genome = form()["form"]["properties"]["genome_selection"]
    assert genome["properties"]["genome_source"]["enum"] == ["igenomes", "dataset"]
    custom = genome["dependencies"]["genome_source"]["oneOf"][1]
    assert custom["required"] == ["genome_index"]
    index = custom["properties"]["genome_index"]
    assert index["format"] == "directory-path" and "pathType" not in index
    assert index["process"] == ["process-cirro-genome-index-bwa-1-0", "genome_bwa_index"]
    mapping = json.loads((HERE / "process-input.json").read_text())
    for key in ("genome_source", "genome", "genome_index"):
        assert mapping[key].startswith("$.dataset.params.genome_selection.")
    assert "|" not in mapping["genome_index"]


def test_form_and_input_expose_the_mode_selection():
    props = form()["form"]["properties"]
    mode = props["mode_selection"]["properties"]["mode"]
    assert mode["enum"] == ["apply_bqsr_and_quantize", "apply_bqsr", "quantize"]
    assert mode["default"] == "apply_bqsr_and_quantize"
    branches = {b["properties"]["mode"]["enum"][0]: b["properties"] for b in props["mode_selection"]["dependencies"]["mode"]["oneOf"]}
    assert set(branches) == set(mode["enum"])
    assert "applybqsr_gpus" in branches["apply_bqsr_and_quantize"] and "applybqsr_gpus" in branches["apply_bqsr"]
    assert "applybqsr_gpus" not in branches["quantize"] and "alignment_source" not in branches["quantize"]
    assert "quantize_quals_enabled" not in props["quantization"]["properties"]
    mapping = json.loads((HERE / "process-input.json").read_text())
    assert mapping["mode"] == "$.dataset.params.mode_selection.mode"
    assert "alignment_source" not in mapping and "quantize_quals_enabled" not in mapping


def test_samplesheet_is_a_structured_input_and_replaces_the_dataset_selector():
    spec = form()
    assert spec["form"]["properties"]["input"] == {"type": "string", "schema": "pipeline_samplesheet"}
    assert spec["form"]["required"] == ["input"]
    assert json.loads((HERE / "process-input.json").read_text())["input"] == "$.structuredInputs.pipeline_samplesheet.path"
    items = spec["structuredInputs"]["pipeline_samplesheet"]["schema"]["items"]
    assert items["required"] == ["sample", "alignment"]
    assert set(items["properties"]) == {"sample", "alignment", "alignment_index", "recal_table"}
    for column in ("alignment", "alignment_index", "recal_table"):
        assert items["properties"][column]["format"] == "file-path"


def test_column_patterns_accept_the_expected_files_and_reject_others():
    props = form()["structuredInputs"]["pipeline_samplesheet"]["schema"]["items"]["properties"]
    ok = lambda column, value: re.search(props[column]["pattern"], value) is not None
    assert ok("alignment", "s3://b/ds/data/preprocessing/parabricks/S1/S1.bam") and ok("alignment", "S1.cram")
    assert not ok("alignment", "S1.bam.bai") and not ok("alignment", "s3://b/my file.bam") and not ok("alignment", "")
    assert ok("alignment_index", "S1.bam.bai") and ok("alignment_index", "S1.cram.crai") and ok("alignment_index", "")
    assert not ok("alignment_index", "S1.bam")
    assert ok("recal_table", "s3://b/S1.table") and ok("recal_table", "") and not ok("recal_table", "S1.tsv")
