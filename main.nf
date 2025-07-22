include { QUALITY_CONTROL } from './modules/quality_control.nf'
include { ASSEMBLY } from './modules/assembly.nf'
include { MERGE } from './modules/merge.nf'
include { CREATE_DB } from './modules/create_db.nf'
include { ALIGN } from './modules/align.nf'
include {TAX} from './modules/tax.nf'

workflow {
    input_ch = channel.fromPath(params.input_csv)
                .splitCsv(header:true)
                .map {row -> 
                tuple(row.sample, [row.r1, row.r2])}
    db_ch = file(params.db)
    CREATE_DB(db_ch)

    
    
    QUALITY_CONTROL(input_ch)
    if (params.mode == 1){
        ASSEMBLY(QUALITY_CONTROL.out.sample)
        sequences = ASSEMBLY.out.contigs
    }
    else{
        MERGE(QUALITY_CONTROL.out.sample)
        sequences = MERGE.out.sequence
    }

    ALIGN(sequences, CREATE_DB.out)

    TAX(ALIGN.out.blast_result)

}
