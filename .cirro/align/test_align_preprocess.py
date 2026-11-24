import logging

import pandas as pd
import pytest

import importlib.util
import pathlib

_spec = importlib.util.spec_from_file_location("align_preprocess", pathlib.Path(__file__).with_name("preprocess.py"))
_preprocess = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_preprocess)
apply_genome_params = _preprocess.apply_genome_params
build_fastq_samplesheet = _preprocess.build_fastq_samplesheet

LOG = logging.getLogger("test")


def table(rows):
    columns = ["patient", "sample", "status", "fastq_1", "fastq_2", "fastq_singleton"]
    return pd.DataFrame([tuple(r) + ("",) * (len(columns) - len(r)) for r in rows], columns=columns)


def test_rows_become_lanes_numbered_per_sample_with_status_default():
    sheet = build_fastq_samplesheet(table([
        ("P1", "S1", "1", "s3://b/b_R1.fq.gz", "s3://b/b_R2.fq.gz"), ("P1", "S1", "1", "s3://b/a_R1.fq.gz", "s3://b/a_R2.fq.gz"),
        ("P2", "S2", "", "", "", "s3://b/c.fq.gz"),
    ]), LOG)
    assert list(sheet.columns) == ["patient", "sample", "status", "lane", "fastq_1", "fastq_2"]
    assert sheet[sheet["sample"] == "S1"][["lane", "fastq_1"]].values.tolist() == [["0", "s3://b/a_R1.fq.gz"], ["1", "s3://b/b_R1.fq.gz"]]
    s2 = sheet[sheet["sample"] == "S2"].iloc[0]
    assert s2["patient"] == "P2" and s2["status"] == 0 and s2["fastq_1"] == "s3://b/c.fq.gz" and s2["fastq_2"] == "" and s2["lane"] == "0"
    assert sheet[sheet["sample"] == "S1"]["status"].tolist() == [1, 1]


def test_lane_numbers_do_not_depend_on_the_other_samples_in_the_run():
    s1 = [("P1", "S1", "", "s3://b/a_R1.fq.gz", "s3://b/a_R2.fq.gz"), ("P1", "S1", "", "s3://b/b_R1.fq.gz", "s3://b/b_R2.fq.gz")]
    s2 = [("P2", "S2", "", "s3://b/c_R1.fq.gz", "s3://b/c_R2.fq.gz")]
    alone = build_fastq_samplesheet(table(s2), LOG)
    together = build_fastq_samplesheet(table(s1 + s2), LOG)
    assert alone[["sample", "lane"]].values.tolist() == [["S2", "0"]]
    assert together[together["sample"] == "S2"][["sample", "lane"]].values.tolist() == [["S2", "0"]]
    assert together[together["sample"] == "S1"]["lane"].tolist() == ["0", "1"]


def test_a_lane_with_a_pair_and_singletons_becomes_a_paired_row_and_a_singleton_row_after_the_paired_lanes():
    sheet = build_fastq_samplesheet(table([
        ("P1", "S1", "", "s3://b/S1_R1.fq.gz", "s3://b/S1_R2.fq.gz", "s3://b/S1__singletons.fastq.gz"),
        ("P1", "S1", "", "s3://b/S1b_R1.fq.gz", "s3://b/S1b_R2.fq.gz"),
    ]), LOG)
    assert sheet[["lane", "fastq_1", "fastq_2"]].values.tolist() == [
        ["0", "s3://b/S1_R1.fq.gz", "s3://b/S1_R2.fq.gz"], ["1", "s3://b/S1b_R1.fq.gz", "s3://b/S1b_R2.fq.gz"],
        ["2", "s3://b/S1__singletons.fastq.gz", ""]]


def test_a_singleton_only_row_is_a_singleton_lane():
    sheet = build_fastq_samplesheet(table([("P1", "S1", "", "", "", "s3://b/S1_se.fastq.gz")]), LOG)
    assert sheet[["lane", "fastq_1", "fastq_2"]].values.tolist() == [["0", "s3://b/S1_se.fastq.gz", ""]]


def test_a_table_without_the_singleton_column_still_builds():
    sheet = build_fastq_samplesheet(pd.DataFrame([("P1", "S1", "", "s3://b/a_R1.fq.gz", "s3://b/a_R2.fq.gz")],
                                                 columns=["patient", "sample", "status", "fastq_1", "fastq_2"]), LOG)
    assert sheet[["lane", "fastq_1", "fastq_2"]].values.tolist() == [["0", "s3://b/a_R1.fq.gz", "s3://b/a_R2.fq.gz"]]


def test_r1_without_r2_is_an_error_that_points_to_the_singleton_column():
    with pytest.raises(ValueError, match=r"both be set or both be blank in row\(s\) 1, 2.*fastq_singleton"):
        build_fastq_samplesheet(table([("P1", "S1", "", "s3://b/a.fq.gz", ""), ("P1", "S1", "", "", "s3://b/b_R2.fq.gz")]), LOG)


def test_a_row_without_any_fastq_is_an_error():
    with pytest.raises(ValueError, match=r"no FASTQ in row\(s\) 2"):
        build_fastq_samplesheet(table([("P1", "S1", "", "s3://b/a_R1.fq.gz", "s3://b/a_R2.fq.gz"), ("P1", "S1", "", "", "")]), LOG)


def test_blank_patient_is_an_error_with_no_fallback_to_the_sample():
    with pytest.raises(ValueError, match="patient is blank in row"):
        build_fastq_samplesheet(table([("", "S1", "", "s3://b/a_R1.fq.gz", "s3://b/a_R2.fq.gz")]), LOG)


def test_one_sample_cannot_have_two_patients_or_two_statuses():
    rows = [("P1", "S1", "0", "s3://b/a_R1.fq.gz", "s3://b/a_R2.fq.gz"), ("P2", "S1", "1", "s3://b/b_R1.fq.gz", "s3://b/b_R2.fq.gz")]
    with pytest.raises(ValueError) as error:
        build_fastq_samplesheet(table(rows), LOG)
    assert "S1: more than one patient" in str(error.value) and "S1: more than one status" in str(error.value)


def test_a_file_listed_twice_is_an_error_not_a_silent_duplicate():
    rows = [("P1", "S1", "", "s3://b/a_R1.fq.gz", "s3://b/a_R2.fq.gz"), ("P1", "S1", "", "s3://b/a_R1.fq.gz", "s3://b/b_R2.fq.gz")]
    with pytest.raises(ValueError, match="a_R1.fq.gz is listed 2 times"):
        build_fastq_samplesheet(table(rows), LOG)
    rows = [("P1", "S1", "", "s3://b/a_R1.fq.gz", "s3://b/a_R2.fq.gz", "s3://b/a_R2.fq.gz")]
    with pytest.raises(ValueError, match="a_R2.fq.gz is listed 2 times"):
        build_fastq_samplesheet(table(rows), LOG)


def test_status_other_than_0_or_1_is_an_error():
    with pytest.raises(ValueError, match="status must be 0"):
        build_fastq_samplesheet(table([("P1", "S1", "2", "s3://b/a_R1.fq.gz", "s3://b/a_R2.fq.gz")]), LOG)


def test_blank_sample_is_an_error():
    with pytest.raises(ValueError, match=r"sample is blank in row\(s\) 1"):
        build_fastq_samplesheet(table([("P1", "", "", "s3://b/a_R1.fq.gz", "s3://b/a_R2.fq.gz")]), LOG)


def test_empty_table_is_an_error():
    with pytest.raises(ValueError, match="no rows"):
        build_fastq_samplesheet(table([]), LOG)


def test_whitespace_around_cells_is_ignored():
    sheet = build_fastq_samplesheet(table([(" P1 ", " S1 ", " 1 ", " s3://b/a_R1.fq.gz ", " s3://b/a_R2.fq.gz ")]), LOG)
    assert sheet.iloc[0].tolist() == ["P1", "S1", 1, "0", "s3://b/a_R1.fq.gz", "s3://b/a_R2.fq.gz"]


def test_igenomes_with_intervals_unticked():
    p = apply_genome_params({"genome_source": "igenomes", "genome": "GATK.GRCh38", "use_intervals": False,
                             "igenomes_base": "s3://refs/igenomes/"})
    assert p["no_intervals"] is True and "genome_source" not in p and "use_intervals" not in p
    assert p["genome"] == "GATK.GRCh38"


def test_custom_genome_resolves_dataset_files_and_known_sites():
    p = apply_genome_params({"genome_source": "dataset", "genome_index": "s3://b/ds/data/",
                             "dbsnp": "s3://r/a/germline_resource.vcf.gz", "known_indels": "s3://r/b/germline_resource.vcf.gz",
                             "custom_intervals": "s3://r/c/regions.bed", "use_intervals": True})
    assert p["genome"] == "null"
    assert p["bwa_index"] == "s3://b/ds/data"
    assert p["ref_fasta"] == "s3://b/ds/data/genome.fasta" and p["ref_fasta_fai"] == "s3://b/ds/data/genome.fasta.fai"
    assert p["known_sites"] == "s3://r/a/germline_resource.vcf.gz,s3://r/b/germline_resource.vcf.gz"
    assert p["intervals"] == "s3://r/c/regions.bed"
    assert not {"genome_source", "genome_index", "dbsnp", "known_indels", "custom_intervals", "use_intervals"} & p.keys()


def test_custom_genome_without_index_is_an_error():
    with pytest.raises(ValueError, match="BWA genome index"):
        apply_genome_params({"genome_source": "dataset"})


def _fields(node, key=None):
    if isinstance(node, dict):
        yield key, node
        for k, value in node.items():
            yield from _fields(value, k)
    elif isinstance(node, list):
        for value in node:
            yield from _fields(value, key)


def test_fastq_columns_accept_reads_and_a_blank_cell():
    import json
    import re
    form = json.loads(pathlib.Path(__file__).with_name("process-form.json").read_text())
    columns = form["structuredInputs"]["pipeline_samplesheet"]["schema"]["items"]["properties"]
    patterns = [re.compile(columns[k]["pattern"]) for k in ("fastq_1", "fastq_2", "fastq_singleton")]
    for name in ("s3://b/d/S1_R1.fastq.gz", "s3://b/d/S1_2.fq", "/data/S1_R2.fq.gz", ""):
        assert all(r.search(name) for r in patterns), name
    for name in ("s3://b/d/S1_I1_001.fastq.gz", "s3://b/d/S1_I2_001.fq.gz", "s3://b/d/S 1_R1.fastq.gz", "s3://b/d/S1.bam", "s3://b/d/"):
        assert not any(r.search(name) for r in patterns), name


def test_form_is_in_path_mode_with_a_structured_samplesheet_and_no_legacy_pickers():
    import json
    form = json.loads(pathlib.Path(__file__).with_name("process-form.json").read_text())
    legacy = [k for k, n in _fields(form["form"]) if "pathType" in n or "file" in n and isinstance(n["file"], str) or "useS3Path" in n]
    assert not legacy, legacy
    assert form["form"]["required"] == ["input"]
    assert form["form"]["properties"]["input"] == {"type": "string", "schema": "pipeline_samplesheet"}
    columns = form["structuredInputs"]["pipeline_samplesheet"]["schema"]["items"]["properties"]
    assert list(columns) == ["patient", "sample", "status", "fastq_1", "fastq_2", "fastq_singleton"]
    assert columns["fastq_1"]["format"] == columns["fastq_2"]["format"] == columns["fastq_singleton"]["format"] == "file-path"
    assert form["structuredInputs"]["pipeline_samplesheet"]["schema"]["items"]["required"] == ["patient", "sample"]
    inputs = json.loads(pathlib.Path(__file__).with_name("process-input.json").read_text())
    assert inputs["input"] == "$.structuredInputs.pipeline_samplesheet.path"
    assert not [k for k, v in inputs.items() if v.endswith(("|/data", "|/config/manifest.csv"))]


def test_every_reference_picker_is_an_end_anchored_file_path():
    import json
    form = json.loads(pathlib.Path(__file__).with_name("process-form.json").read_text())
    refs = {k: n for k, n in _fields(form["form"]) if k in ("dbsnp", "known_indels", "custom_intervals") and "format" in n}
    assert len(refs) == 3
    for name, n in refs.items():
        assert n["format"] == "file-path" and n["pattern"].endswith("$") and not n["pattern"].startswith("^"), name


def test_every_dataset_folder_picker_is_filtered_by_process():
    import json
    form = json.loads(pathlib.Path(__file__).with_name("process-form.json").read_text())
    folders = {k: n for k, n in _fields(form["form"]) if n.get("format") == "directory-path"}
    assert folders, "expected at least one folder picker"
    assert all(f.get("process") for f in folders.values()), [k for k, f in folders.items() if not f.get("process")]
    assert folders["genome_index"]["process"] == ["process-cirro-genome-index-bwa-1-0", "genome_bwa_index"]


def test_form_offers_low_memory_off_and_gpuwrite_on_by_default():
    # gpuwrite stays on (no runtime difference), but a deep exome failed in PBCuBGZFWriter with
    # "Index buffer is too small..." under both low-memory settings, so it must be switchable.
    import json
    here = pathlib.Path(__file__).parent
    form = (here / "process-form.json").read_text()
    resources = json.loads(form)["form"]["properties"]["alignment_resources"]["properties"]
    assert resources["fq2bam_low_memory"]["default"] is False
    assert resources["fq2bam_gpuwrite"]["default"] is True
    inputs = json.loads((here / "process-input.json").read_text())
    assert inputs["fq2bam_gpuwrite"] == "$.dataset.params.alignment_resources.fq2bam_gpuwrite"


def test_igenomes_with_custom_intervals_uses_them():
    p = apply_genome_params({"genome_source": "igenomes", "genome": "GATK.GRCh38",
                             "intervals_mode": "custom", "custom_intervals": "s3://r/c/regions.bed"})
    assert p["intervals"] == "s3://r/c/regions.bed" and "no_intervals" not in p
    assert "intervals_mode" not in p and "custom_intervals" not in p


def test_igenomes_intervals_none_sets_no_intervals():
    p = apply_genome_params({"genome_source": "igenomes", "genome": "GATK.GRCh38", "intervals_mode": "none"})
    assert p["no_intervals"] is True and "intervals" not in p


def test_igenomes_default_intervals_leaves_the_genome_default():
    p = apply_genome_params({"genome_source": "igenomes", "genome": "GATK.GRCh38", "intervals_mode": "gatk_calling_regions"})
    assert "intervals" not in p and "no_intervals" not in p


def test_igenomes_custom_intervals_without_a_file_is_an_error():
    with pytest.raises(ValueError, match="custom intervals"):
        apply_genome_params({"genome_source": "igenomes", "genome": "GATK.GRCh38", "intervals_mode": "custom"})


def test_gpu_field_defaults_to_2_and_caps_at_4():
    import json
    adv = json.loads((pathlib.Path(__file__).parent / "process-form.json").read_text())["form"]["properties"]["alignment_resources"]["properties"]
    assert adv["fq2bam_gpus"]["default"] == 2 and adv["fq2bam_gpus"]["maximum"] == 4


def test_igenomes_branch_offers_three_interval_modes_and_a_custom_bed():
    import json
    form = json.loads((pathlib.Path(__file__).parent / "process-form.json").read_text())
    igenomes = form["form"]["properties"]["genome_selection"]["dependencies"]["genome_source"]["oneOf"][0]
    mode = igenomes["properties"]["intervals_mode"]
    assert mode["enum"] == ["gatk_calling_regions", "none", "custom"] and mode["default"] == "gatk_calling_regions"
    custom = [b for b in igenomes["dependencies"]["intervals_mode"]["oneOf"] if b["properties"]["intervals_mode"]["enum"] == ["custom"]][0]
    bed = custom["properties"]["custom_intervals"]
    assert bed["format"] == "file-path" and custom["required"] == ["custom_intervals"]


def test_form_does_not_expose_the_fq2bam_intermediate_format():
    here = pathlib.Path(__file__).parent
    assert "fq2bam_intermediate_fmt" not in (here / "process-form.json").read_text()
    assert "fq2bam_intermediate_fmt" not in (here / "process-input.json").read_text()


def test_form_offers_an_optional_fq2bam_memory_override_and_maps_it():
    import json
    here = pathlib.Path(__file__).parent
    advanced = json.loads((here / "process-form.json").read_text())["form"]["properties"]["alignment_resources"]
    field = advanced["properties"]["fq2bam_memory_gb"]
    assert field["type"] == "integer" and field["minimum"] == 16 and field["maximum"] == 768
    assert "default" not in field and "fq2bam_memory_gb" not in advanced.get("required", [])
    mapping = json.loads((here / "process-input.json").read_text())
    assert mapping["fq2bam_memory_gb"] == "$.dataset.params.alignment_resources.fq2bam_memory_gb"


def test_form_offers_poly_g_trimming_auto_by_default_and_maps_it():
    import json
    here = pathlib.Path(__file__).parent
    trimming = json.loads((here / "process-form.json").read_text())["form"]["properties"]["read_trimming"]["properties"]
    assert trimming["poly_g_trimming"]["default"] == "auto"
    assert trimming["poly_g_trimming"]["enum"] == ["auto", "on", "off"]
    assert "trim_nextseq" not in trimming
    mapping = json.loads((here / "process-input.json").read_text())
    assert mapping["poly_g_trimming"] == "$.dataset.params.read_trimming.poly_g_trimming"
    assert "trim_nextseq" not in mapping


def test_form_offers_publish_markduplicates_off_by_default_and_maps_it():
    import json
    here = pathlib.Path(__file__).parent
    props = json.loads((here / "process-form.json").read_text())["form"]["properties"]
    field = props["publish_markduplicates"]
    assert field["type"] == "boolean" and field["default"] is False
    mapping = json.loads((here / "process-input.json").read_text())
    assert mapping["publish_markduplicates"] == "$.dataset.params.publish_markduplicates"


def test_form_offers_postalt_off_by_default_in_advanced_and_maps_it():
    import json
    here = pathlib.Path(__file__).parent
    props = json.loads((here / "process-form.json").read_text())["form"]["properties"]
    field = props["advanced"]["properties"]["postalt"]
    assert field["type"] == "boolean" and field["default"] is False
    mapping = json.loads((here / "process-input.json").read_text())
    assert mapping["postalt"] == "$.dataset.params.advanced.postalt"


def test_form_offers_apply_bqsr_on_by_default_and_maps_it():
    import json
    here = pathlib.Path(__file__).parent
    props = json.loads((here / "process-form.json").read_text())["form"]["properties"]
    field = props["bqsr"]["properties"]["apply_bqsr"]
    assert field["type"] == "boolean" and field["default"] is True
    mapping = json.loads((here / "process-input.json").read_text())
    assert mapping["apply_bqsr"] == "$.dataset.params.bqsr.apply_bqsr"


def test_quantization_section_only_offered_when_apply_bqsr_is_on():
    import json
    group = json.loads((pathlib.Path(__file__).parent / "process-form.json").read_text())["form"]["properties"]["bqsr"]
    assert "quantization" not in group["properties"]
    branches = group["dependencies"]["apply_bqsr"]["oneOf"]
    on_branch = [b for b in branches if b["properties"]["apply_bqsr"]["enum"] == [True]][0]
    off_branch = [b for b in branches if b["properties"]["apply_bqsr"]["enum"] == [False]][0]
    assert on_branch["properties"]["quantization"]["properties"]["quantize_quals_enabled"]["default"] is True
    assert "quantization" not in off_branch["properties"]


def test_bam_is_the_first_output_format_and_cram_version_shows_only_for_cram():
    import json
    group = json.loads((pathlib.Path(__file__).parent / "process-form.json").read_text())["form"]["properties"]["output"]
    assert group["properties"]["output_fmt"]["enum"] == ["bam", "cram"]
    assert group["properties"]["output_fmt"]["default"] == "bam"
    assert "cram_version" not in group["properties"]
    branches = {b["properties"]["output_fmt"]["enum"][0]: b for b in group["dependencies"]["output_fmt"]["oneOf"]}
    assert "cram_version" not in branches["bam"]["properties"]
    assert branches["cram"]["properties"]["cram_version"]["default"] == "3.0"


def test_alignment_resources_are_their_own_block_above_advanced():
    import json
    props = json.loads((pathlib.Path(__file__).parent / "process-form.json").read_text())["form"]["properties"]
    keys = list(props)
    assert keys.index("alignment_resources") == keys.index("advanced") - 1
    assert list(props["alignment_resources"]["properties"]) == ["fq2bam_gpus", "fq2bam_low_memory", "fq2bam_gpuwrite", "fq2bam_memory_gb"]
    assert "fq2bam_gpus" not in props["advanced"]["properties"]


def test_clip_settings_only_offered_when_trimming_is_on():
    import json
    here = pathlib.Path(__file__).parent
    trim = json.loads((here / "process-form.json").read_text())["form"]["properties"]["read_trimming"]["properties"]["trim"]
    assert list(trim["properties"]) == ["trim_fastq"]
    branches = {b["properties"]["trim_fastq"]["enum"][0]: b for b in trim["dependencies"]["trim_fastq"]["oneOf"]}
    assert set(branches[True]["properties"]) == {"trim_fastq", "clip_r1", "clip_r2", "three_prime_clip_r1", "three_prime_clip_r2"}
    assert set(branches[False]["properties"]) == {"trim_fastq"}
    mapping = json.loads((here / "process-input.json").read_text())
    assert mapping["clip_r1"] == "$.dataset.params.read_trimming.trim.clip_r1"
