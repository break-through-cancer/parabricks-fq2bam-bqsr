import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]


def test_postalt_is_off_by_default_and_the_container_is_pinned():
    config = (ROOT / "nextflow.config").read_text()
    assert re.search(r"^params\.postalt\s*=\s*false$", config, re.M)
    image = re.search(r"^params\.postalt_container\s*=\s*'([^']+)'$", config, re.M).group(1)
    assert re.search(r":\d+\.\d+\.\d+$", image), "the postalt image must carry an explicit version tag"
    dockerfile = (ROOT / "tools" / "postalt" / "Dockerfile").read_text()
    assert f'org.opencontainers.image.version="{image.rsplit(":", 1)[1]}"' in dockerfile
    assert f"IMAGE ?= {image}" in (ROOT / "tools" / "postalt" / "Makefile").read_text()


def test_every_download_in_the_postalt_image_is_checksummed():
    dockerfile = (ROOT / "tools" / "postalt" / "Dockerfile").read_text()
    downloads = len(re.findall(r"curl -fsSL", dockerfile))
    checks = len(re.findall(r"sha256sum -c", dockerfile))
    assert downloads == checks == 2
    assert re.search(r"ARG POSTALT_RS_COMMIT=[0-9a-f]{40}\n", dockerfile)
    assert "cargo build --release --locked" in dockerfile


def test_postalt_image_ships_the_gpl_notices_and_names_its_source():
    dockerfile = (ROOT / "tools" / "postalt" / "Dockerfile").read_text()
    assert "/src/LICENSE /src/NOTICE.md /usr/share/doc/postalt_rs/" in dockerfile
    assert 'org.opencontainers.image.licenses="GPL-3.0 AND MIT"' in dockerfile
    assert "postalt_rs.source=" in dockerfile and "${POSTALT_RS_COMMIT}" in dockerfile


def test_fq2bam_and_markdup_stop_publishing_md_directly_when_postalt_runs():
    # With postalt on, the published .md alignment is POSTALT's output, not fq2bam's.
    text = (ROOT / "conf" / "modules.config").read_text()
    for name in ("PARABRICKS_FQ2BAM", "PARABRICKS_MARKDUP"):
        block = text[text.index(f"withName: '{name}'"):]
        block = block[:block.index("withName:", 20)] if "withName:" in block[20:] else block
        assert "meta.publish_directly == true && !params.postalt" in block, name
    postalt = text[text.index("withName: 'POSTALT'"):]
    postalt = postalt[:postalt.index("withName:", 20)]
    assert "meta.publish_directly == true ?" in postalt
    assert "reports/postalt/" in postalt


def test_postalt_is_wired_between_alignment_and_bqsr_and_only_for_the_fastq_entry():
    main = (ROOT / "main.nf").read_text()
    assert main.index("POSTALT(") < main.index("PARABRICKS_APPLYBQSR(\n                aligned_cram.join")
    assert "--postalt applies to the FASTQ entry only" in main
    assert "postalt_added" in main


def test_cirro_compute_config_retries_postalt_with_more_memory():
    config = (ROOT / ".cirro" / "align" / "process-compute.config").read_text()
    block = config[config.index("withName: 'POSTALT'"):]
    assert "memory = { 32.GB * task.attempt }" in block[:block.index("}") + 1]
