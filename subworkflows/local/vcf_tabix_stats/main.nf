//
// Run BCFTools tabix and stats commands
//

include { TABIX_TABIX    } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_STATS } from '../../../modules/nf-core/bcftools/stats/main'

workflow VCF_TABIX_STATS {
    take:
    vcf     // channel: [ val(meta), [ vcf ] ]
    regions //    file: regions.txt
    targets //    file: targets.txt
    samples //    file: samples.txt

    main:

    ch_vcf_for_tabix = vcf.map { meta, vcf_file -> [ meta, vcf_file, [], [] ] }

    TABIX_TABIX (
        ch_vcf_for_tabix
    )

    BCFTOOLS_STATS (
        vcf.join(TABIX_TABIX.out.index, by: [0]),
        regions,
        targets,
        samples,
        [ [:], [] ],
        [ [:], [] ]
    )

    emit:
    tbi      = TABIX_TABIX.out.index    // channel: [ val(meta), [ tbi ] ]

    stats    = BCFTOOLS_STATS.out.stats // channel: [ val(meta), [ txt ] ]

}
