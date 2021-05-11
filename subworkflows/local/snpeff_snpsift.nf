//
// Run snpEff, bgzip, tabix, stats and SnpSift commands
//

params.snpeff_options  = [:]
params.bgzip_options   = [:]
params.tabix_options   = [:]
params.stats_options   = [:]
params.snpsift_options = [:]

include { SNPEFF_ANN            } from '../../modules/local/snpeff_ann'            addParams( options: params.snpeff_options  )
include { SNPSIFT_EXTRACTFIELDS } from '../../modules/local/snpsift_extractfields' addParams( options: params.snpsift_options )
include { VCF_BGZIP_TABIX_STATS } from '../nf-core/vcf_bgzip_tabix_stats'          addParams( bgzip_options: params.bgzip_options, tabix_options: params.tabix_options, stats_options: params.stats_options )

workflow SNPEFF_SNPSIFT {
    take:
    vcf    // channel: [ val(meta), [ vcf ] ]
    db     // path   : snpEff database
    config // path   : snpEff config
    fasta  // path   : genome.fasta

    main:
    
    SNPEFF_ANN ( vcf, db, config, fasta )

    VCF_BGZIP_TABIX_STATS ( SNPEFF_ANN.out.vcf )

    SNPSIFT_EXTRACTFIELDS ( VCF_BGZIP_TABIX_STATS.out.vcf )

    emit:
    csv              = SNPEFF_ANN.out.csv                         // channel: [ val(meta), [ csv ] ]
    txt              = SNPEFF_ANN.out.txt                         // channel: [ val(meta), [ txt ] ]
    html             = SNPEFF_ANN.out.html                        // channel: [ val(meta), [ html ] ]
    snpeff_version   = SNPEFF_ANN.out.version                     //    path: *.version.txt

    vcf              = VCF_BGZIP_TABIX_STATS.out.vcf              // channel: [ val(meta), [ vcf.gz ] ]
    tbi              = VCF_BGZIP_TABIX_STATS.out.tbi              // channel: [ val(meta), [ tbi ] ]
    stats            = VCF_BGZIP_TABIX_STATS.out.stats            // channel: [ val(meta), [ txt ] ]
    tabix_version    = VCF_BGZIP_TABIX_STATS.out.tabix_version    //    path: *.version.txt
    bcftools_version = VCF_BGZIP_TABIX_STATS.out.bcftools_version //    path: *.version.txt

    snpsift_txt      = SNPSIFT_EXTRACTFIELDS.out.txt              // channel: [ val(meta), [ txt ] ]
    snpsift_version  = SNPSIFT_EXTRACTFIELDS.out.version          //    path: *.version.txt
}
