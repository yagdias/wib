process ASSEMBLY {
    publishDir "results/megahit/${sample}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.contigs.fa"), emit: contigs

    script:
    """
    megahit -1 ${fastq[0]} -2 ${fastq[1]} -o ${sample} -t ${task.cpus}
    mv ${sample}/final.contigs.fa ${sample}.contigs.fa
    """
}