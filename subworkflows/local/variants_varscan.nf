/*
 * Variant calling and downstream processing for Varscan
 */

params.varscan_mpileup2cns_options   = [:]
params.quast_options                 = [:]
params.bcftools_filter_options       = [:]
params.consensus_genomecov_options   = [:]
params.consensus_merge_options       = [:]
params.consensus_mask_options        = [:]
params.consensus_maskfasta_options   = [:]
params.consensus_bcftools_options    = [:]
params.consensus_plot_options        = [:]
params.bgzip_options                 = [:]
params.tabix_options                 = [:]
params.stats_options                 = [:]
params.bcftools_filter_tabix_options = [:]
params.bcftools_filter_stats_options = [:]
params.snpeff_lowfreq_options        = [:]
params.snpsift_lowfreq_options       = [:]
params.snpeff_lowfreq_bgzip_options  = [:]
params.snpeff_lowfreq_tabix_options  = [:]
params.snpeff_lowfreq_stats_options  = [:]
params.snpeff_highfreq_options       = [:]
params.snpsift_highfreq_options      = [:]
params.snpeff_highfreq_bgzip_options = [:]
params.snpeff_highfreq_tabix_options = [:]
params.snpeff_highfreq_stats_options = [:]

include { VARSCAN_MPILEUP2CNS                       } from '../../modules/local/varscan_mpileup2cns' addParams( options: params.varscan_mpileup2cns_options )
include { BCFTOOLS_FILTER                           } from '../../modules/local/bcftools_filter'     addParams( options: params.bcftools_filter_options     )
include { MAKE_CONSENSUS                            } from './make_consensus'                        addParams( genomecov_options: params.consensus_genomecov_options, merge_options: params.consensus_merge_options, mask_options: params.consensus_mask_options, maskfasta_options: params.consensus_maskfasta_options, bcftools_options: params.consensus_bcftools_options, plot_bases_options: params.consensus_plot_options )
include { VCF_BGZIP_TABIX_STATS                     } from './vcf_bgzip_tabix_stats'                 addParams( bgzip_options: params.bgzip_options, tabix_options: params.tabix_options, stats_options: params.stats_options )
include { VCF_TABIX_STATS                           } from './vcf_tabix_stats'                       addParams( tabix_options: params.bcftools_filter_tabix_options, stats_options: params.bcftools_filter_stats_options      )
include { SNPEFF_SNPSIFT as SNPEFF_SNPSIFT_LOWFREQ  } from './snpeff_snpsift'                        addParams( snpeff_options: params.snpeff_lowfreq_options, snpsift_options: params.snpsift_lowfreq_options, bgzip_options: params.snpeff_lowfreq_bgzip_options, tabix_options: params.snpeff_lowfreq_tabix_options, stats_options: params.snpeff_lowfreq_stats_options      )
include { SNPEFF_SNPSIFT as SNPEFF_SNPSIFT_HIGHFREQ } from './snpeff_snpsift'                        addParams( snpeff_options: params.snpeff_highfreq_options, snpsift_options: params.snpsift_highfreq_options, bgzip_options: params.snpeff_highfreq_bgzip_options, tabix_options: params.snpeff_highfreq_tabix_options, stats_options: params.snpeff_highfreq_stats_options )
include { QUAST                                     } from '../../modules/nf-core/quast/main'        addParams( options: params.quast_options )

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
     * Filter VCF for high frequency variants
     */
    BCFTOOLS_FILTER ( VCF_BGZIP_TABIX_STATS.out.vcf )

    /*
     * Index VCF for high frequency variants
     */
    VCF_TABIX_STATS ( BCFTOOLS_FILTER.out.vcf )

    /*
     * Create genome consensus using variants in VCF
     */
    if (!params.skip_consensus) {
        MAKE_CONSENSUS ( bam.join(BCFTOOLS_FILTER.out.vcf, by: [0]).join(VCF_TABIX_STATS.out.tbi, by: [0]), fasta )

        if (!params.skip_variants_quast) {
            QUAST ( MAKE_CONSENSUS.out.fasta.collect{ it[1] }, fasta, gff, true, if params.gff )
        }
    }

    /*
     * Annotate low/high frequency variants
     */
    if (params.gff && !params.skip_variants_snpeff) {
        SNPEFF_SNPSIFT_LOWFREQ ( VCF_BGZIP_TABIX_STATS.out.vcf, snpeff_db, snpeff_config, fasta )

        SNPEFF_SNPSIFT_HIGHFREQ ( BCFTOOLS_FILTER.out.vcf, snpeff_db, snpeff_config, fasta )
    }

    // emit:
    // scaffolds         = SPADES.out.scaffolds              // channel: [ val(meta), [ scaffolds ] ]
    // spades_version    = SPADES.out.version                //    path: *.version.txt

}