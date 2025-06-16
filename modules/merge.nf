process MERGE {
    publishDir "results/merge/${sample}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.fasta"), emit: sequence

    script:
    """
    vsearch --fastq_mergepairs ${fastq[0]} --reverse ${fastq[1]} \\
    --threads ${task.cpus} --fastqout ${sample}.fastq

    seqtk seq ${sample}.fastq -A > ${sample}.fasta
    """
}