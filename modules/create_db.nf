process CREATE_DB {
    publishDir { "${projectDir}/results/tmp/${db}" }, mode: "copy", overwrite: true

    input:
    file db

    output:
    tuple path(db), path("${db}.*")

    script:
    """
    makeblastdb -dbtype nucl -in ${db} -taxid_map ${params.taxon_names} -parse_seqids
    """
}
