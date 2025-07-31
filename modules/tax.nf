process TAX {
    publishDir { "${projectDir}/results/tax/${sample}" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(blast_output)

    output:
    tuple val(sample), path("${sample}.csv"), emit: tax

    script:
    """
    tax.py \\
    --blast_file $blast_output \\
    --sample_name $sample \\
    --pident ${params.min_pident} \\
    --taxonomy_db_path ${params.taxonomy_db_path} --cut_evalue ${params.cut_evalue} --update_db false \\
    """

}
