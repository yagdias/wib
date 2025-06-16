process CREATE_DB {
    publishDir "./results/tmp/${db}", mode: "copy", overwrite: true

    input:
    file db

    output:
    file "${db}.*"

    script:
    """
    diamond makedb --db ${db} --in ${db} --threads ${task.cpus}
    """
}