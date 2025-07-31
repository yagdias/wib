process QUALITY_CONTROL {
    publishDir { "${projectDir}/results/fastp/${sample}" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path('*.fq.gz'), emit: sample
    tuple val(sample), path('*.json'), emit: json
    tuple val(sample), path('*.html'), emit: html



    script: // colocar o argumento de qualidade
    """
    fastp \\
    -i ${fastq[0]} \\
    -I ${fastq[1]} \\
    -w ${task.cpus} \\
    -o ${sample}.R1.fq.gz \\
    -O ${sample}.R2.fq.gz \\
    -h ${sample}.fastp.html \\
    -j ${sample}.fastp.json \\
    --qualified_quality_phred ${params.phred_score} \\
    """
}
