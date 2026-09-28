//
// Create a long table with variant information including AA changes and lineage info
//

include { BCFTOOLS_QUERY           } from '../../../modules/nf-core/bcftools/query/main'
include { MAKE_VARIANTS_LONG_TABLE } from '../../../modules/local/make_variants_long_table'

workflow VARIANTS_LONG_TABLE {
    take:
    vcf      // channel: [ val(meta), [ vcf ] ]
    tbi      // channel: [ val(meta), [ tbi ] ]
    snpsift  // channel: [ val(meta), [ txt ] ]
    pangolin // channel: [ val(meta), [ csv ] ]

    main:

    BCFTOOLS_QUERY (
        vcf.join(tbi, by: [0]),
        [],
        [],
        []
    )

    MAKE_VARIANTS_LONG_TABLE (
        BCFTOOLS_QUERY.out.output.collect{_meta, output -> output},
        snpsift.collect{_meta, txt -> txt}.ifEmpty([]),
        pangolin.collect{_meta, csv -> csv}.ifEmpty([])
    )

    emit:
    query_table = BCFTOOLS_QUERY.out.output        // channel: [ val(meta), [ txt ] ]
    long_table  = MAKE_VARIANTS_LONG_TABLE.out.csv // channel: [ val(meta), [ csv ] ]
}
