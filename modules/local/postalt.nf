process POSTALT {
    tag "${meta.sample}"
    container params.postalt_container
    label 'process_high'

    input:
    tuple val(meta), path(alignment, stageAs: 'input/*'), path(index, stageAs: 'input/*')
    tuple path(ref_fasta), path(ref_fasta_fai)
    path alt_file
    val output_fmt

    output:
    tuple val(meta), path("${meta.sample}.md.${output_fmt}"), path("${meta.sample}.md.${output_fmt}.${output_fmt == 'cram' ? 'crai' : 'bai'}"), emit: alignment
    tuple val(meta), path("${meta.sample}.md.postalt.json"), emit: json
    tuple val(meta), path("${meta.sample}.md.postalt_added.txt"), emit: added
    path "${meta.sample}.md.postalt_mqc.tsv", emit: mqc

    script:
    """
    set -euo pipefail
    postalt_apply \\
        --in ${alignment} \\
        --ref ${ref_fasta} \\
        --alt ${alt_file} \\
        --out ${meta.sample}.md \\
        --fmt ${output_fmt} \\
        --cram-version ${params.cram_version} \\
        --threads ${task.cpus} \\
        --mem-gb ${task.memory.toGiga()}
    """

    stub:
    def idx = output_fmt == 'cram' ? 'crai' : 'bai'
    """
    touch ${meta.sample}.md.${output_fmt} ${meta.sample}.md.${output_fmt}.${idx}
    printf '{"added_records": 0}\\n' > ${meta.sample}.md.postalt.json
    echo 0 > ${meta.sample}.md.postalt_added.txt
    printf '# id: postalt\\n# section_name: Alt-aware postalt\\n# plot_type: table\\nSample\\tadded_records\\n${meta.sample}\\t0\\n' > ${meta.sample}.md.postalt_mqc.tsv
    """
}
