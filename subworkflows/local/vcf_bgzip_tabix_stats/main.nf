//
// Run BCFTools bgzip, tabix and stats commands
//

include { TABIX_BGZIP     } from '../../../modules/nf-core/tabix/bgzip/main'
include { VCF_TABIX_STATS } from '../vcf_tabix_stats'

workflow VCF_BGZIP_TABIX_STATS {
    take:
    vcf     // channel: [ val(meta), [ vcf ] ]
    regions //    file: regions.txt
    targets //    file: targets.txt
    samples //    file: samples.txt

    main:

    TABIX_BGZIP (
        vcf
    )

    VCF_TABIX_STATS (
        TABIX_BGZIP.out.output,
        regions,
        targets,
        samples
    )

    emit:
    vcf      = TABIX_BGZIP.out.output    // channel: [ val(meta), [ vcf.gz ] ]

    tbi      = VCF_TABIX_STATS.out.tbi   // channel: [ val(meta), [ tbi ] ]
    stats    = VCF_TABIX_STATS.out.stats // channel: [ val(meta), [ txt ] ]
}
