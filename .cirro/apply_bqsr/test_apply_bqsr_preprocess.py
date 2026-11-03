import json
import pytest

import importlib.util
import pathlib

_spec = importlib.util.spec_from_file_location("apply_bqsr_preprocess", pathlib.Path(__file__).with_name("preprocess.py"))
_preprocess = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_preprocess)
build_samplesheet = _preprocess.build_samplesheet
apply_genome_params = _preprocess.apply_genome_params

ROOT = "s3://bucket/datasets/ds1/data/preprocessing"


def files(*names):
    return [f"{ROOT}/{n}" for n in names]


def test_one_sample_bam():
    sheet = build_samplesheet(files(
        "parabricks/S1/S1.bam", "parabricks/S1/S1.bam.bai", "parabricks/S1/S1.table",
    ))
    assert list(sheet.columns) == ["sample", "alignment", "alignment_index", "recal_table"]
    assert sheet.to_dict("records") == [{
        "sample": "S1",
        "alignment": f"{ROOT}/parabricks/S1/S1.bam",
        "alignment_index": f"{ROOT}/parabricks/S1/S1.bam.bai",
        "recal_table": f"{ROOT}/parabricks/S1/S1.table",
    }]


def test_cram_and_multiple_samples_sorted():
    sheet = build_samplesheet(files(
        "parabricks/S2/S2.cram", "parabricks/S2/S2.cram.crai", "parabricks/S2/S2.table",
        "parabricks/S1/S1.bam", "parabricks/S1/S1.bam.bai", "parabricks/S1/S1.table",
    ))
    assert list(sheet["sample"]) == ["S1", "S2"]
    assert sheet.loc[1, "alignment_index"].endswith("S2.cram.crai")


def test_recalibrated_stage_is_ignored():
    # Applying the table again to an already-recalibrated alignment would recalibrate it twice.
    sheet = build_samplesheet(files(
        "parabricks/S1/S1.bam", "parabricks/S1/S1.bam.bai", "parabricks/S1/S1.table",
        "recalibrated/S1/S1.recal.bam", "recalibrated/S1/S1.recal.bam.bai",
    ))
    assert sheet["alignment"].tolist() == [f"{ROOT}/parabricks/S1/S1.bam"]


def test_our_own_markduplicates_layout_bam():
    # This pipeline's own apply_bqsr=false output: alignment+index under markduplicates/,
    # the table separately under recal_table/ -- not co-located like sarek_align's layout.
    sheet = build_samplesheet(files(
        "markduplicates/S1/S1.md.bam", "markduplicates/S1/S1.md.bam.bai", "recal_table/S1/S1.table",
    ))
    assert sheet.to_dict("records") == [{
        "sample": "S1",
        "alignment": f"{ROOT}/markduplicates/S1/S1.md.bam",
        "alignment_index": f"{ROOT}/markduplicates/S1/S1.md.bam.bai",
        "recal_table": f"{ROOT}/recal_table/S1/S1.table",
    }]


def test_our_own_markduplicates_layout_cram():
    sheet = build_samplesheet(files(
        "markduplicates/S1/S1.md.cram", "markduplicates/S1/S1.md.cram.crai", "recal_table/S1/S1.table",
    ))
    assert sheet.loc[0, "alignment"].endswith("S1.md.cram")
    assert sheet.loc[0, "alignment_index"].endswith("S1.md.cram.crai")


def test_sarek_and_our_own_layout_in_the_same_batch():
    sheet = build_samplesheet(files(
        "parabricks/S1/S1.bam", "parabricks/S1/S1.bam.bai", "parabricks/S1/S1.table",
        "markduplicates/S2/S2.md.bam", "markduplicates/S2/S2.md.bam.bai", "recal_table/S2/S2.table",
    ))
    assert list(sheet["sample"]) == ["S1", "S2"]


def test_no_recognized_files():
    with pytest.raises(ValueError, match="preprocessing/parabricks"):
        build_samplesheet(files("recalibrated/S1/S1.recal.bam", "recalibrated/S1/S1.recal.bam.bai"))


def test_missing_table_and_index_are_reported_per_sample():
    with pytest.raises(ValueError) as err:
        build_samplesheet(files(
            "parabricks/S1/S1.bam", "parabricks/S1/S1.bam.bai",
            "parabricks/S2/S2.cram", "parabricks/S2/S2.table",
        ))
    msg = str(err.value)
    assert "S1: no .table recalibration table" in msg
    assert "S2: no index for" in msg


def test_two_alignments_for_one_sample_is_an_error():
    with pytest.raises(ValueError, match="S1: expected one BAM/CRAM, found 2"):
        build_samplesheet(files(
            "parabricks/S1/S1.bam", "parabricks/S1/S1.bam.bai",
            "parabricks/S1/S1.cram", "parabricks/S1/S1.cram.crai",
            "parabricks/S1/S1.table",
        ))


def test_same_sample_in_two_input_datasets_is_an_error():
    other = "s3://bucket/datasets/ds2/data/preprocessing"
    with pytest.raises(ValueError, match="S1: expected one BAM/CRAM, found 2"):
        build_samplesheet(files(
            "parabricks/S1/S1.bam", "parabricks/S1/S1.bam.bai", "parabricks/S1/S1.table",
        ) + [f"{other}/parabricks/S1/S1.bam", f"{other}/parabricks/S1/S1.bam.bai", f"{other}/parabricks/S1/S1.table"])


def test_bam_is_the_first_output_format_and_cram_version_shows_only_for_cram():
    import json
    group = json.loads((pathlib.Path(__file__).parent / "process-form.json").read_text())["form"]["properties"]["output"]
    assert group["properties"]["output_fmt"]["enum"] == ["bam", "cram"]
    assert group["properties"]["output_fmt"]["default"] == "bam"
    assert "cram_version" not in group["properties"]
    branches = {b["properties"]["output_fmt"]["enum"][0]: b for b in group["dependencies"]["output_fmt"]["oneOf"]}
    assert "cram_version" not in branches["bam"]["properties"]
    assert branches["cram"]["properties"]["cram_version"]["default"] == "3.0"


def test_igenomes_is_the_default_genome():
    p = apply_genome_params({"genome": "GATK.GRCh38", "igenomes_base": "s3://refs/igenomes/"})
    assert p["genome"] == "GATK.GRCh38" and "ref_fasta" not in p
    p = apply_genome_params({"genome_source": "igenomes", "genome": "GATK.GRCh38"})
    assert p["genome"] == "GATK.GRCh38" and "genome_source" not in p


def test_custom_genome_resolves_dataset_files():
    p = apply_genome_params({"genome_source": "dataset", "genome_index": "s3://b/ds/data"})
    assert p["genome"] == "null"
    assert p["bwa_index"] == "s3://b/ds/data"
    assert p["ref_fasta"] == "s3://b/ds/data/genome.fasta" and p["ref_fasta_fai"] == "s3://b/ds/data/genome.fasta.fai"
    assert not {"genome_source", "genome_index"} & p.keys()


def test_custom_genome_without_index_is_an_error():
    with pytest.raises(ValueError, match="BWA genome index"):
        apply_genome_params({"genome_source": "dataset"})


def test_form_and_input_expose_the_genome_selection():
    here = pathlib.Path(__file__).parent
    form = json.loads((here / "process-form.json").read_text())["form"]["properties"]["genome_selection"]
    assert form["properties"]["genome_source"]["enum"] == ["igenomes", "dataset"]
    custom = form["dependencies"]["genome_source"]["oneOf"][1]
    assert custom["required"] == ["genome_index"]
    assert custom["properties"]["genome_index"]["process"] == ["process-cirro-genome-index-bwa-1-0", "genome_bwa_index"]
    mapping = json.loads((here / "process-input.json").read_text())
    for key in ("genome_source", "genome", "genome_index"):
        assert mapping[key].startswith("$.dataset.params.genome_selection.")
    assert mapping["genome_index"].endswith("|/data")
