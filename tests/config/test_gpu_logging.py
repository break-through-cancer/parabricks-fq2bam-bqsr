import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]


def test_parabricks_tasks_log_gpu_model_and_memory_before_pbrun():
    # GPU type and memory decide whether --low-memory is needed and make A/B timings comparable.
    for nf in sorted((ROOT / "modules" / "local").glob("*.nf")):
        text = nf.read_text()
        if "clara-parabricks" not in text:
            continue
        script = text.split("script:", 1)[1].split("stub:", 1)[0]
        gpu = script.find("nvidia-smi --query-gpu=index,name,memory.total")
        pbrun = script.find("pbrun ")
        assert gpu != -1, f"{nf.name} does not log the GPU"
        assert gpu < pbrun, f"{nf.name} logs the GPU after pbrun"


def test_parabricks_tasks_log_the_driver_version():
    # CUDA 13 images need a recent host driver, so a failed pull-and-run should show which driver it met.
    for nf in sorted((ROOT / "modules" / "local").glob("*.nf")):
        text = nf.read_text()
        if "clara-parabricks" in text:
            assert "memory.total,driver_version" in text, nf.name


def test_every_parabricks_module_uses_the_same_image_tag():
    tags = set()
    for nf in sorted((ROOT / "modules" / "local").glob("*.nf")):
        tags.update(re.findall(r"clara-parabricks:([\w.\-]+)", nf.read_text()))
    assert tags == {"4.8.0-1"}
