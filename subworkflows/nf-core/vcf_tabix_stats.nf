//
// Run BCFTools tabix and stats commands
//

include { TABIX_TABIX    } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_STATS } from '../../modules/nf-core/bcftools/stats/main'

workflow VCF_TABIX_STATS {
    take:
    vcf     // channel: [ val(meta), [ vcf ] ]
    regions //    file: regions.txt
    targets //    file: targets.txt
    samples //    file: samples.txt

    main:

    ch_versions = Channel.empty()

    TABIX_TABIX (
        vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    BCFTOOLS_STATS (
        vcf.join(TABIX_TABIX.out.tbi, by: [0]),
        regions,
        targets,
        samples
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())

    emit:
    tbi      = TABIX_TABIX.out.tbi      // channel: [ val(meta), [ tbi ] ]
    csi      = TABIX_TABIX.out.csi      // channel: [ val(meta), [ csi ] ]

    stats    = BCFTOOLS_STATS.out.stats // channel: [ val(meta), [ txt ] ]

    versions = ch_versions              // channel: [ versions.yml ]

}
