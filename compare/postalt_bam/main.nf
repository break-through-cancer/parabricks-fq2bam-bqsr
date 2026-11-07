nextflow.enable.dsl = 2

params.input         = null
params.outdir        = './results'
params.igenomes_base = 's3://ngi-igenomes/igenomes/'
params.container     = 'ghcr.io/break-through-cancer/parabricks-fq2bam-bqsr-postalt:0.2.0'
params.postalt       = false

process MAKE_BAM {
    tag { sample }
    container params.container
    cpus 16
    memory { 64.GB * task.attempt }
    publishDir { "${params.outdir}/preprocessing/recalibrated/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(alignment), path(alignment_index)
    tuple path(ref_fasta), path(ref_fai), path(alt_file)

    output:
    tuple val(sample), path("${sample}.recal.bam"), path("${sample}.recal.bam.bai")

    script:
    if (params.postalt)
        """
        set -euo pipefail
        postalt_apply --in ${alignment} --ref ${ref_fasta} --alt ${alt_file} --out ${sample}.recal \\
            --fmt bam --threads ${task.cpus} --mem-gb ${task.memory.toGiga()}
        """
    else
        """
        set -euo pipefail
        samtools view -b --no-PG -@ ${task.cpus} -T ${ref_fasta} -o ${sample}.recal.bam ${alignment}
        samtools index -@ ${task.cpus} ${sample}.recal.bam
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
    MAKE_BAM(samples, ref_ch)
}
