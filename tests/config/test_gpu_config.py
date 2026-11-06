import os
import pathlib
import re
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[2]


def flat_config(*extra):
    out = subprocess.run(["nextflow", *extra, "config", "-flat"], cwd=ROOT, capture_output=True, text=True,
                         env={**os.environ, "PW_ONDEMAND_JOB_QUEUE": "q"}, check=True)
    return out.stdout


def parabricks_processes():
    names = []
    for nf in (ROOT / "modules" / "local").glob("*.nf"):
        text = nf.read_text()
        if "clara-parabricks" in text:
            names += re.findall(r"^process\s+(\w+)", text, re.M)
    return sorted(names)


def test_every_parabricks_process_requests_a_gpu():
    procs = parabricks_processes()
    assert procs == ["PARABRICKS_APPLYBQSR", "PARABRICKS_FQ2BAM", "PARABRICKS_FQ2BAM_PART", "PARABRICKS_MARKDUP"]
    for cfg in (flat_config(), flat_config("-c", ".cirro/align/process-compute.config")):
        for p in procs:
            assert re.search(rf"withName:{p}'\.accelerator = ", cfg), f"{p} has no accelerator"


def test_applybqsr_requests_at_most_two_gpus():
    for cfg in (flat_config(), flat_config("-c", ".cirro/align/process-compute.config"),
                flat_config("-c", ".cirro/applybqsr_quantize/process-compute.config")):
        assert re.search(r"withName:PARABRICKS_APPLYBQSR'\.accelerator = ", cfg)
    assert "applybqsr_gpus = 1" in flat_config()
    assert "in 1..2" in (ROOT / "main.nf").read_text()


def test_applybqsr_gpu_count_is_a_parameter_limited_to_what_pbrun_accepts():
    import json
    import pathlib

    root = pathlib.Path(__file__).resolve().parents[2]
    assert "accelerator      = { params.applybqsr_gpus }" in (root / "conf" / "modules.config").read_text()
    assert "accelerator" not in (root / ".cirro" / "applybqsr_quantize" / "process-compute.config").read_text()
    form = json.loads((root / ".cirro" / "applybqsr_quantize" / "process-form.json").read_text())["form"]["properties"]
    gpus = form["mode_selection"]["dependencies"]["mode"]["oneOf"][0]["properties"]["applybqsr_gpus"]
    assert (gpus["default"], gpus["minimum"], gpus["maximum"]) == (1, 1, 2)
    mapping = json.loads((root / ".cirro" / "applybqsr_quantize" / "process-input.json").read_text())
    assert mapping["applybqsr_gpus"] == "$.dataset.params.mode_selection.applybqsr_gpus"
