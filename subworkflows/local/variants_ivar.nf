/*
 * Variant calling and downstream processing for IVar
 */

params.ivar_variants_options        = [:]
params.ivar_variants_to_vcf_options = [:]
params.tabix_bgzip_options          = [:]
params.tabix_tabix_options          = [:]
params.bcftools_stats_options       = [:]
params.ivar_consensus_options       = [:]
params.consensus_plot_options       = [:]
params.quast_options                = [:]
params.snpeff_options               = [:]
params.snpsift_options              = [:]
params.snpeff_bgzip_options         = [:]
params.snpeff_tabix_options         = [:]
params.snpeff_stats_options         = [:]
params.pangolin_options             = [:]

include { IVAR_VARIANTS_TO_VCF  } from '../../modules/local/ivar_variants_to_vcf'           addParams( options: params.ivar_variants_to_vcf_options )
include { PLOT_BASE_DENSITY     } from '../../modules/local/plot_base_density'              addParams( options: params.consensus_plot_options       )
include { IVAR_VARIANTS         } from '../../modules/nf-core/software/ivar/variants/main'  addParams( options: params.ivar_variants_options        )
include { IVAR_CONSENSUS        } from '../../modules/nf-core/software/ivar/consensus/main' addParams( options: params.ivar_consensus_options       )
include { QUAST                 } from '../../modules/nf-core/software/quast/main'          addParams( options: params.quast_options                )
include { PANGOLIN              } from '../../modules/nf-core/software/pangolin/main'       addParams( options: params.pangolin_options             )
include { VCF_BGZIP_TABIX_STATS } from '../nf-core/vcf_bgzip_tabix_stats'                   addParams( bgzip_options: params.tabix_bgzip_options, tabix_options: params.tabix_tabix_options, stats_options: params.bcftools_stats_options )
include { SNPEFF_SNPSIFT        } from './snpeff_snpsift'                                   addParams( snpeff_options: params.snpeff_options, snpsift_options: params.snpsift_options, bgzip_options: params.snpeff_bgzip_options, tabix_options: params.snpeff_tabix_options, stats_options:  params.snpeff_stats_options )

workflow VARIANTS_IVAR {
    take:
    bam                 // channel: [ val(meta), [ bam ] ]
    fasta               // channel: /path/to/genome.fasta
    gff                 // channel: /path/to/genome.gff
    snpeff_db           // channel: /path/to/snpeff_db/
    snpeff_config       // channel: /path/to/snpeff.config
    ivar_multiqc_header // channel: /path/to/multiqc_header for ivar variants
    
    main:
    /*
     * Call variants
     */
    IVAR_VARIANTS ( bam, fasta, gff )

    /*
     * Convert original iVar output to VCF, zip and index
     */
    IVAR_VARIANTS_TO_VCF ( IVAR_VARIANTS.out.tsv, ivar_multiqc_header )

    VCF_BGZIP_TABIX_STATS ( IVAR_VARIANTS_TO_VCF.out.vcf )

    /*
     * Create genome consensus
     */
    ch_consensus        = Channel.empty()
    ch_consensus_qual   = Channel.empty()
    ch_bases_tsv        = Channel.empty()
    ch_bases_pdf        = Channel.empty()
    ch_quast_results    = Channel.empty()
    ch_quast_tsv        = Channel.empty()
    ch_quast_version    = Channel.empty()
    ch_pangolin_report  = Channel.empty()
    ch_pangolin_version = Channel.empty()
    if (!params.skip_consensus) {
        IVAR_CONSENSUS ( bam, fasta )
        ch_consensus       = IVAR_CONSENSUS.out.fasta
        ch_consensus_qual  = IVAR_CONSENSUS.out.qual
    
        PLOT_BASE_DENSITY ( ch_consensus )
        ch_bases_tsv = PLOT_BASE_DENSITY.out.tsv
        ch_bases_pdf = PLOT_BASE_DENSITY.out.pdf

        if (!params.skip_variants_quast) {
            QUAST ( ch_consensus.collect{ it[1] }, fasta, gff, true, params.gff )
            ch_quast_results = QUAST.out.results
            ch_quast_tsv     = QUAST.out.tsv
            ch_quast_version = QUAST.out.version
        }

        if (!params.skip_pangolin) {
            PANGOLIN ( ch_consensus )
            ch_pangolin_report  = PANGOLIN.out.report
            ch_pangolin_version = PANGOLIN.out.version
        }
    }

    /*
     * Annotate variants
     */
    ch_snpeff_vcf      = Channel.empty()
    ch_snpeff_tbi      = Channel.empty()
    ch_snpeff_stats    = Channel.empty()
    ch_snpeff_csv      = Channel.empty()
    ch_snpeff_txt      = Channel.empty()
    ch_snpeff_html     = Channel.empty()
    ch_snpsift_txt     = Channel.empty()
    ch_snpeff_version  = Channel.empty()
    ch_snpsift_version = Channel.empty()
    if (params.gff && !params.skip_snpeff) {
        SNPEFF_SNPSIFT ( VCF_BGZIP_TABIX_STATS.out.vcf, snpeff_db, snpeff_config, fasta )
        ch_snpeff_vcf      = SNPEFF_SNPSIFT.out.vcf
        ch_snpeff_tbi      = SNPEFF_SNPSIFT.out.tbi
        ch_snpeff_stats    = SNPEFF_SNPSIFT.out.stats
        ch_snpeff_csv      = SNPEFF_SNPSIFT.out.csv
        ch_snpeff_txt      = SNPEFF_SNPSIFT.out.txt
        ch_snpeff_html     = SNPEFF_SNPSIFT.out.html
        ch_snpsift_txt     = SNPEFF_SNPSIFT.out.snpsift_txt
        ch_snpeff_version  = SNPEFF_SNPSIFT.out.snpeff_version
        ch_snpsift_version = SNPEFF_SNPSIFT.out.snpsift_version
    }

    emit:
    tsv              = IVAR_VARIANTS.out.tsv                      // channel: [ val(meta), [ tsv ] ]
    ivar_version     = IVAR_VARIANTS.out.version                  //    path: *.version.txt

    vcf_orig         = IVAR_VARIANTS_TO_VCF.out.vcf               // channel: [ val(meta), [ vcf ] ]
    log_out          = IVAR_VARIANTS_TO_VCF.out.log               // channel: [ val(meta), [ log ] ]
    multiqc_tsv      = IVAR_VARIANTS_TO_VCF.out.tsv               // channel: [ val(meta), [ tsv ] ]

    vcf              = VCF_BGZIP_TABIX_STATS.out.vcf              // channel: [ val(meta), [ vcf ] ]
    tbi              = VCF_BGZIP_TABIX_STATS.out.tbi              // channel: [ val(meta), [ tbi ] ]
    stats            = VCF_BGZIP_TABIX_STATS.out.stats            // channel: [ val(meta), [ txt ] ]
    tabix_version    = VCF_BGZIP_TABIX_STATS.out.tabix_version    //    path: *.version.txt
    bcftools_version = VCF_BGZIP_TABIX_STATS.out.bcftools_version //    path: *.version.txt
    
    consensus        = ch_consensus                               // channel: [ val(meta), [ fasta ] ]
    consensus_qual   = ch_consensus_qual                          // channel: [ val(meta), [ fasta ] ]
    bases_tsv        = ch_bases_tsv                               // channel: [ val(meta), [ tsv ] ]
    bases_pdf        = ch_bases_pdf                               // channel: [ val(meta), [ pdf ] ]
    
    quast_results    = ch_quast_results                           // channel: [ val(meta), [ results ] ]
    quast_tsv        = ch_quast_tsv                               // channel: [ val(meta), [ tsv ] ]
    quast_version    = ch_quast_version                           //    path: *.version.txt

    snpeff_vcf       = ch_snpeff_vcf                              // channel: [ val(meta), [ vcf.gz ] ]
    snpeff_tbi       = ch_snpeff_tbi                              // channel: [ val(meta), [ tbi ] ]
    snpeff_stats     = ch_snpeff_stats                            // channel: [ val(meta), [ txt ] ]
    snpeff_csv       = ch_snpeff_csv                              // channel: [ val(meta), [ csv ] ]
    snpeff_txt       = ch_snpeff_txt                              // channel: [ val(meta), [ txt ] ]
    snpeff_html      = ch_snpeff_html                             // channel: [ val(meta), [ html ] ]
    snpsift_txt      = ch_snpsift_txt                             // channel: [ val(meta), [ txt ] ]
    snpeff_version   = ch_snpeff_version                          //    path: *.version.txt
    snpsift_version  = ch_snpsift_version                         //    path: *.version.txt

    pangolin_report  = ch_pangolin_report                         // channel: [ val(meta), [ csv ] ]
    pangolin_version = ch_pangolin_version                        //    path: *.version.txt
}

