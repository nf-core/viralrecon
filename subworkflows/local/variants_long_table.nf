//
// Create a long table with variant information including AA changes and lineage info
//

include { BCFTOOLS_QUERY           } from '../../modules/nf-core/bcftools/query/main'
include { MAKE_VARIANTS_LONG_TABLE } from '../../modules/local/make_variants_long_table'

workflow VARIANTS_LONG_TABLE {
    take:
    vcf      // channel: [ val(meta), [ vcf ] ]
    tbi      // channel: [ val(meta), [ tbi ] ]
    snpsift  // channel: [ val(meta), [ txt ] ]
    pangolin // channel: [ val(meta), [ csv ] ]

    main:

    ch_versions = Channel.empty()

    BCFTOOLS_QUERY (
        vcf.join(tbi, by: [0]),
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())

    MAKE_VARIANTS_LONG_TABLE (
        BCFTOOLS_QUERY.out.txt.collect{it[1]},
        snpsift.collect{it[1]}.ifEmpty([]),
        pangolin.collect{it[1]}.ifEmpty([])
    )
    ch_versions = ch_versions.mix(MAKE_VARIANTS_LONG_TABLE.out.versions)

    emit:
    query_table = BCFTOOLS_QUERY.out.txt           // channel: [ val(meta), [ txt ] ]
    long_table  = MAKE_VARIANTS_LONG_TABLE.out.csv // channel: [ val(meta), [ csv ] ]

    versions    = ch_versions    // channel: [ versions.yml ]
}
