//
// Run snpEff, bgzip, tabix, stats and SnpSift commands
//

include { SNPEFF_ANN            } from '../../../modules/local/snpeff/ann'
include { SNPSIFT_EXTRACTFIELDS } from '../../../modules/local/snpsift/extractfields'

include { VCF_BGZIP_TABIX_STATS } from '../vcf_bgzip_tabix_stats'

workflow SNPEFF_SNPSIFT {
    take:
    vcf    // channel: [ val(meta), [ vcf ] ]
    db     // path   : snpEff database
    config // path   : snpEff config
    fasta  // path   : genome.fasta

    main:

    SNPEFF_ANN (
        vcf,
        db,
        config,
        fasta
    )

    VCF_BGZIP_TABIX_STATS (
        SNPEFF_ANN.out.vcf,
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ]
    )

    SNPSIFT_EXTRACTFIELDS (
        VCF_BGZIP_TABIX_STATS.out.vcf
    )

    emit:
    csv         = SNPEFF_ANN.out.csv              // channel: [ val(meta), [ csv ] ]
    txt         = SNPEFF_ANN.out.txt              // channel: [ val(meta), [ txt ] ]
    html        = SNPEFF_ANN.out.html             // channel: [ val(meta), [ html ] ]

    vcf         = VCF_BGZIP_TABIX_STATS.out.vcf   // channel: [ val(meta), [ vcf.gz ] ]
    tbi         = VCF_BGZIP_TABIX_STATS.out.tbi   // channel: [ val(meta), [ tbi ] ]
    stats       = VCF_BGZIP_TABIX_STATS.out.stats // channel: [ val(meta), [ txt ] ]

    snpsift_txt = SNPSIFT_EXTRACTFIELDS.out.txt   // channel: [ val(meta), [ txt ] ]
}
