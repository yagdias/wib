process TAX {
    input:
    tuple val(sample), path(blast_output)

    output:
    tuple val(sample), path("${sample}.fasta"), emit: sequence

    script:
    """
    tax.py \\
    --blast_file $blast_output \\
    --sample_name $sample \\
    --pident ${params.min_pident} \\
    --taxonomy_db_path ${params.taxonomy_db_path} --cut_evalue ${params.cut_evalue} --update_db true \\
    """

}