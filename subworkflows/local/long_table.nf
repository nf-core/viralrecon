//
// Variant calling QC
//
include { BCFTOOLS_QUERY } from '../../modules/nf-core/modules/bcftools/query/main'
include { CREATE_LONG_TABLE } from '../../modules/local/create_long_table'

workflow LONG_TABLE {
    take:
    vcf           // channel: [ val(meta), [ vcf ] ]
    tbi           // channel: [ val(meta), [ tbi ] ]
    pangolin      // channel: [ val(meta), [ ] ]
    snpsift       // channel: [ val(meta), [ ] ]

    main:

    ch_versions = Channel.empty()

    if (!params.skip_long_table) {
        BCFTOOLS_QUERY (
            vcf.join(tbi, by: [0]),
            "",
            "",
            ""
        )
        ch_query_table   = BCFTOOLS_QUERY.out.vcf
        ch_versions     = ch_versions.mix(BCFTOOLS_QUERY.out.versions)

        CREATE_LONG_TABLE (
            ch_query_table.collect{it[1]},
            snpsift.collect{it[1]},
            pangolin.collect{it[1]}
        )
        ch_long_table = CREATE_LONG_TABLE.out.csv_variants
    }

    emit:
    longtable_csv      = ch_long_table      // channel: [ val(meta), [ vcf.gz ] ]
    query_table        = ch_query_table     // channel: [ val(meta), [ txt ] ]

    versions           = ch_versions        // channel: [ versions.yml ]
}
