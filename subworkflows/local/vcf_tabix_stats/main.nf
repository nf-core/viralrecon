//
// Run BCFTools tabix and stats commands
//

include { HTSLIB_BGZIPTABIX } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { BCFTOOLS_STATS    } from '../../../modules/nf-core/bcftools/stats/main'

workflow VCF_TABIX_STATS {
    take:
    vcf     // channel: [ val(meta), [ vcf ] ]
    regions //    file: regions.txt
    targets //    file: targets.txt
    samples //    file: samples.txt

    main:

    HTSLIB_BGZIPTABIX (
        vcf.map { meta, vcf_file -> [ meta, vcf_file, [], [] ] },
        'compress',
        true,
        'vcf'
    )

    BCFTOOLS_STATS (
        vcf.join(HTSLIB_BGZIPTABIX.out.index, by: [0]),
        regions,
        targets,
        samples,
        [ [:], [] ],
        [ [:], [] ]
    )

    emit:
    tbi      = HTSLIB_BGZIPTABIX.out.index    // channel: [ val(meta), [ tbi ] ]

    stats    = BCFTOOLS_STATS.out.stats // channel: [ val(meta), [ txt ] ]

}
