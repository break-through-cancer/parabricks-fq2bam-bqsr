nextflow.enable.dsl = 2

params.input         = null
params.outdir        = './results'
params.igenomes_base = 's3://ngi-igenomes/igenomes/'
params.container     = 'ghcr.io/break-through-cancer/postalt-compare:0.1.0'
params.records_per_partition = 400000

process COMPARE_ENGINES {
    tag { sample }
    container params.container
    cpus 16
    memory { 48.GB * task.attempt }
    publishDir { "${params.outdir}/reports/postalt_compare/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(alignment), path(alignment_index)
    tuple path(ref_fasta), path(ref_fai), path(alt_file)

    output:
    tuple val(sample), path("${sample}/summary.json"), path("${sample}/stage_times.txt")

    script:
    """
    set -euo pipefail
    compare_engines.sh ${alignment} ${ref_fasta} ${alt_file} ${sample} ${task.cpus} ${params.records_per_partition} 8
    rm -rf ${sample}/work
    """
}

workflow {
    if (!params.input) error "--input samplesheet (sample,alignment,alignment_index) is required"
    def base = params.igenomes_base.toString().replaceAll('/+\$', '')
    def fasta = "${base}/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
    def alt   = "${base}/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.alt"
    ref_ch = channel.value([file(fasta, checkIfExists: true), file("${fasta}.fai", checkIfExists: true), file(alt, checkIfExists: true)])
    samples = channel.fromPath(params.input).splitCsv(header: true).map { row ->
        tuple(row.sample, file(row.alignment, checkIfExists: true), file(row.alignment_index, checkIfExists: true))
    }
    COMPARE_ENGINES(samples, ref_ch)
}
