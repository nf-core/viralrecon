//
// Run BCFTools bgzip, tabix and stats commands
//

include { TABIX_BGZIP     } from '../../modules/nf-core/modules/tabix/bgzip/main'
include { VCF_TABIX_STATS } from './vcf_tabix_stats'

workflow VCF_BGZIP_TABIX_STATS {
    take:
    vcf // channel: [ val(meta), [ vcf ] ]

    main:

    ch_versions = Channel.empty()

    TABIX_BGZIP (
        vcf
    )
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())

    VCF_TABIX_STATS (
        TABIX_BGZIP.out.output
    )
    ch_versions = ch_versions.mix(VCF_TABIX_STATS.out.versions)

    emit:
    vcf      = TABIX_BGZIP.out.output    // channel: [ val(meta), [ vcf.gz ] ]
    tbi      = VCF_TABIX_STATS.out.tbi   // channel: [ val(meta), [ tbi ] ]
    stats    = VCF_TABIX_STATS.out.stats // channel: [ val(meta), [ txt ] ]

    versions = ch_versions               // channel: [ versions.yml ]
}
