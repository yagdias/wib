process ALIGN {
    publishDir "./results/blast/${sample}", mode: "copy", overwrite: true

    input:
    tuple val(sample), path(sequence)
    tuple val(db), file("${db}*")

    output:
    tuple val(sample), path("${sample}.blastn"), emit: blast_result

    script:
    """
    blastn \\
            -query ${sequence} \\
            -db ${db} \\
            -outfmt "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms stitle salltitles nident qcovhsp"  \\
            -num_threads ${task.cpus} \\
            -out ${sample}.blastn
    """
}