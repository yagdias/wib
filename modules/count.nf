
process COUNT {
    input:
    tuple val(sample), path(sequences), path(blast_result)

    output:
    path("*.txt"), emit: counts

    script:
    """
    total_seqs=\$(grep -c "^>" ${sequences})
    # Count unique query IDs in blast result (column 1) to get number of classified sequences
    total_hits=\$(cut -f1 ${blast_result} | sort -u | wc -l)
    
    echo "${sample} \${total_seqs} \${total_hits}" > ${sample}_counts.txt
    """
}

process MERGE_COUNTS {
    publishDir "${projectDir}/results", mode: 'copy', overwrite: true

    input:
    path(counts)

    output:
    path("assembly_classification_counts.txt")

    script:
    """
    echo "montados" > assembly_classification_counts.txt
    cat ${counts} | awk '{print \$1, \$2}' >> assembly_classification_counts.txt
    
    echo "classificados" >> assembly_classification_counts.txt
    cat ${counts} | awk '{print \$3, \$1}' >> assembly_classification_counts.txt
    """
}
