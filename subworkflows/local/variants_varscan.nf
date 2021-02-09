/*
 * Variant calling and downstream processing for Varscan
 */

params.varscan_mpileup2cns_options = [:]
params.bcftools_bgzip_options      = [:]
params.bcftools_tabix_options      = [:]
params.bcftools_stats_options      = [:]
params.consensus_genomecov_options = [:]
params.consensus_merge_options     = [:]
params.consensus_mask_options      = [:]
params.consensus_maskfasta_options = [:]
params.consensus_bcftools_options  = [:]
params.consensus_plot_options      = [:]
params.quast_options               = [:]
params.snpeff_options              = [:]
params.snpsift_options             = [:]
params.snpeff_bgzip_options        = [:]
params.snpeff_tabix_options        = [:]
params.snpeff_stats_options        = [:]

include { VARSCAN_MPILEUP2CNS   } from '../../modules/local/varscan_mpileup2cns'   addParams( options: params.varscan_mpileup2cns_options )
include { QUAST                 } from '../../modules/nf-core/software/quast/main' addParams( options: params.quast_options               )
include { VCF_BGZIP_TABIX_STATS } from '../nf-core/vcf_bgzip_tabix_stats'          addParams( bgzip_options: params.bcftools_bgzip_options, tabix_options: params.bcftools_tabix_options, stats_options: params.bcftools_stats_options )
include { MAKE_CONSENSUS        } from './make_consensus'                          addParams( genomecov_options: params.consensus_genomecov_options, merge_options: params.consensus_merge_options, mask_options: params.consensus_mask_options, maskfasta_options: params.consensus_maskfasta_options, bcftools_options: params.consensus_bcftools_options, plot_bases_options: params.consensus_plot_options )
include { SNPEFF_SNPSIFT        } from './snpeff_snpsift'                          addParams( snpeff_options: params.snpeff_options, snpsift_options: params.snpsift_options, bgzip_options: params.snpeff_bgzip_options, tabix_options: params.snpeff_tabix_options, stats_options:  params.snpeff_stats_options )

workflow VARIANTS_VARSCAN {
    take:
    mpileup       // channel: [ val(meta), [ mpileup ] ]
    bam           // channel: [ val(meta), [ bam ] ]
    fasta         // channel: /path/to/genome.fasta
    gff           // channel: /path/to/genome.gff
    snpeff_db     // channel: /path/to/snpeff_db/
    snpeff_config // channel: /path/to/snpeff.config
    
    main:
    /*
     * Call variants
     */
    VARSCAN_MPILEUP2CNS ( mpileup, fasta )

    /*
     * Zip, index VCF output
     */
    VCF_BGZIP_TABIX_STATS ( VARSCAN_MPILEUP2CNS.out.vcf )

    /*
     * Create genome consensus using variants in VCF
     */
    if (!params.skip_consensus) {
        MAKE_CONSENSUS ( bam.join(VCF_BGZIP_TABIX_STATS.out.vcf, by: [0]).join(VCF_BGZIP_TABIX_STATS.out.tbi, by: [0]), fasta )

        if (!params.skip_variants_quast) {
            QUAST ( MAKE_CONSENSUS.out.fasta.collect{ it[1] }, fasta, gff, true, params.gff )
        }
    }

    /*
     * Annotate variants
     */
    if (params.gff && !params.skip_variants_snpeff) {
        SNPEFF_SNPSIFT ( VCF_BGZIP_TABIX_STATS.out.vcf, snpeff_db, snpeff_config, fasta )
    }

    emit:
    vcf_orig         = VARSCAN_MPILEUP2CNS.out.vcf        // channel: [ val(meta), [ vcf ] ]
    log_out          = VARSCAN_MPILEUP2CNS.out.log        // channel: [ val(meta), [ log ] ]
    varscan_version  = VARSCAN_MPILEUP2CNS.out.version    //    path: *.version.txt

    vcf              = VCF_BGZIP_TABIX_STATS.out.vcf      // channel: [ val(meta), [ vcf ] ]
    tbi              = VCF_BGZIP_TABIX_STATS.out.tbi      // channel: [ val(meta), [ tbi ] ]
    stats            = VCF_BGZIP_TABIX_STATS.out.stats    // channel: [ val(meta), [ txt ] ]
    bcftools_version = VCF_BGZIP_TABIX_STATS.out.version  //    path: *.version.txt

    consensus        = MAKE_CONSENSUS.out.fasta            // channel: [ val(meta), [ fasta ] ]
    bases_tsv        = MAKE_CONSENSUS.out.tsv              // channel: [ val(meta), [ tsv ] ]
    bases_pdf        = MAKE_CONSENSUS.out.pdf              // channel: [ val(meta), [ pdf ] ]
    bedtools_version = MAKE_CONSENSUS.out.bedtools_version //    path: *.version.txt 
    
    quast_results    = QUAST.out.results                   // channel: [ val(meta), [ results ] ]
    quast_tsv        = QUAST.out.tsv                       // channel: [ val(meta), [ tsv ] ]
    quast_version    = QUAST.out.version                   //    path: *.version.txt

    snpeff_vcf       = SNPEFF_SNPSIFT.out.vcf              // channel: [ val(meta), [ vcf.gz ] ]
    snpeff_tbi       = SNPEFF_SNPSIFT.out.tbi              // channel: [ val(meta), [ tbi ] ]
    snpeff_stats     = SNPEFF_SNPSIFT.out.stats            // channel: [ val(meta), [ txt ] ]
    snpeff_csv       = SNPEFF_SNPSIFT.out.csv              // channel: [ val(meta), [ csv ] ]
    snpeff_txt       = SNPEFF_SNPSIFT.out.txt              // channel: [ val(meta), [ txt ] ]
    snpeff_html      = SNPEFF_SNPSIFT.out.html             // channel: [ val(meta), [ html ] ]
    snpsift_txt      = SNPEFF_SNPSIFT.out.snpsift_txt      // channel: [ val(meta), [ txt ] ]
    snpeff_version   = SNPEFF_SNPSIFT.out.snpeff_version   //    path: *.version.txt
    snpsift_version  = SNPEFF_SNPSIFT.out.snpsift_version  //    path: *.version.txt
}