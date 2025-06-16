process ALIGN {
    publishDir "./results/diamond/${sample}", mode: "copy", overwrite: true

    input:
    tuple val(sample), path(sequence)
    file db

    output:
    tuple val(sample), path("${sample}.dblastx"), emit: blast_result

    script:
    """
    diamond blastx \\
            --out ${sample}.dblastx \\
            --outfmt 6\\
            --query ${sequence} \\
            --db ${db} \\
            --evalue 0.001 \\
            --threads ${task.cpus} \\
            --ultra-sensitive
    """
}