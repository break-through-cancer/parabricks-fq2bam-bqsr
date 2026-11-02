import pathlib
import re
import shutil
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
MODULES = ROOT / "modules" / "local"


def helper():
    text = (MODULES / "pbrun_wrapper.nf").read_text()
    return re.search(r"def dedupeHdFunction\(\) \{\s*'''(.*?)'''", text, re.S).group(1)


def test_every_stage_that_can_carry_a_queryname_header_repairs_it():
    # Merged split-path parts leave a stale SO:queryname @HD after the SO:coordinate one; smoove and GRIDSS read the last.
    for name, anchor in (("parabricks_markdup.nf", "dedupe_hd ${out}"), ("parabricks_applybqsr.nf", "dedupe_hd ${meta.sample}.recal.bam")):
        text = (MODULES / name).read_text()
        assert "dedupeHdFunction()" in text and anchor in text, name


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not installed")
def test_helper_keeps_only_the_coordinate_header_and_drops_stale_indexes(tmp_path):
    (tmp_path / "in.sam").write_text(
        "@HD\tVN:1.5\tSO:coordinate\n@HD\tVN:1.5\tSO:queryname\n@SQ\tSN:c1\tLN:28\n"
        "r1\t0\tc1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n"
    )
    script = f"""set -euo pipefail
{helper()}
samtools view -b -o t.bam in.sam && samtools index t.bam
dedupe_hd t.bam
samtools view -H t.bam | grep -c '^@HD'
samtools view -H t.bam | grep '^@HD'
test ! -e t.bam.bai
samtools view -c t.bam
dedupe_hd t.bam
"""
    out = subprocess.run(["bash", "-c", script], cwd=tmp_path, capture_output=True, text=True, check=True).stdout.split("\n")
    assert out[:3] == ["1", "@HD\tVN:1.5\tSO:coordinate", "1"]
