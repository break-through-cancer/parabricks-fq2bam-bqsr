nextflow.enable.dsl = 2

params.input         = null
params.outdir        = './results'
params.igenomes_base = 's3://ngi-igenomes/igenomes/'
params.container     = 'ghcr.io/break-through-cancer/parabricks-fq2bam-bqsr-postalt:0.2.0'
params.postalt       = false
params.fmt           = 'cram'
params.region        = ''

process MAKE_ALIGNMENT {
    tag { sample }
    container params.container
    cpus 16
    memory { 64.GB * task.attempt }
    publishDir { "${params.outdir}/preprocessing/recalibrated/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(alignment), path(alignment_index)
    tuple path(ref_fasta), path(ref_fai), path(alt_file)

    output:
    tuple val(sample), path("${sample}.recal.${params.fmt}"), path("${sample}.recal.${params.fmt}.${params.fmt == 'cram' ? 'crai' : 'bai'}")

    script:
    def idx    = params.fmt == 'cram' ? 'crai' : 'bai'
    def ofmt   = params.fmt == 'cram' ? 'cram,version=3.0' : 'bam'
    def source = params.postalt ? "full.cram" : "${alignment}"
    def run_postalt = params.postalt ? """
        postalt_apply --in ${alignment} --ref ${ref_fasta} --alt ${alt_file} --out full \\
            --fmt cram --cram-version 3.0 --threads ${task.cpus} --mem-gb ${task.memory.toGiga()}
    """ : ""
    """
    set -euo pipefail
    ${run_postalt}
    regions=""
    for c in ${params.region}; do
        regions="\$regions \$c \$(cut -f1 ${ref_fai} | grep -E "^\${c}_.*_alt\$" | tr '\\n' ' ')"
    done
    samtools view --no-PG -@ ${task.cpus} -T ${ref_fasta} -O ${ofmt} --write-index \\
        -o ${sample}.recal.${params.fmt}##idx##${sample}.recal.${params.fmt}.${idx} ${source} \$regions
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
    MAKE_ALIGNMENT(samples, ref_ch)
}
