/*
 * Variant calling and downstream processing for IVar
 */

params.ivar_variants_options        = [:]
params.ivar_variants_to_vcf_options = [:]
params.bcftools_bgzip_options       = [:]
params.bcftools_tabix_options       = [:]
params.bcftools_stats_options       = [:]
params.ivar_consensus_options       = [:]
params.consensus_plot_options       = [:]
params.quast_options                = [:]
params.snpeff_options               = [:]
params.snpsift_options              = [:]
params.snpeff_bgzip_options         = [:]
params.snpeff_tabix_options         = [:]
params.snpeff_stats_options         = [:]

include { IVAR_VARIANTS         } from '../../modules/local/ivar_variants'         addParams( options: params.ivar_variants_options        )
include { IVAR_CONSENSUS        } from '../../modules/local/ivar_consensus'        addParams( options: params.ivar_consensus_options       )
include { IVAR_VARIANTS_TO_VCF  } from '../../modules/local/ivar_variants_to_vcf'  addParams( options: params.ivar_variants_to_vcf_options )
include { PLOT_BASE_DENSITY     } from '../../modules/local/plot_base_density'     addParams( options: params.consensus_plot_options       )
include { QUAST                 } from '../../modules/nf-core/software/quast/main' addParams( options: params.quast_options                )
include { VCF_BGZIP_TABIX_STATS } from './vcf_bgzip_tabix_stats'                   addParams( bgzip_options: params.bcftools_bgzip_options, tabix_options: params.bcftools_tabix_options, stats_options: params.bcftools_stats_options )
include { SNPEFF_SNPSIFT        } from './snpeff_snpsift'                          addParams( snpeff_options: params.snpeff_options, snpsift_options: params.snpsift_options, bgzip_options: params.snpeff_bgzip_options, tabix_options: params.snpeff_tabix_options, stats_options:  params.snpeff_stats_options )

workflow VARIANTS_IVAR {
    take:
    mpileup             // channel: [ val(meta), [ mpileup ] ]
    fasta               // channel: /path/to/genome.fasta
    gff                 // channel: /path/to/genome.gff
    snpeff_db           // channel: /path/to/snpeff_db/
    snpeff_config       // channel: /path/to/snpeff.config
    ivar_multiqc_header // channel: /path/to/multiqc_header for ivar variants
    
    main:
    /*
     * Call variants
     */
    IVAR_VARIANTS ( mpileup, fasta, gff )

    /*
     * Convert original iVar output to VCF, zip and index
     */
    IVAR_VARIANTS_TO_VCF ( IVAR_VARIANTS.out.tsv, ivar_multiqc_header )

    VCF_BGZIP_TABIX_STATS ( IVAR_VARIANTS_TO_VCF.out.vcf )

    /*
     * Create genome consensus
     */
    if (!params.skip_consensus) {
        IVAR_CONSENSUS ( mpileup )

        PLOT_BASE_DENSITY ( IVAR_CONSENSUS.out.fasta )

        if (!params.skip_variants_quast) {
            QUAST ( IVAR_CONSENSUS.out.fasta.collect{ it[1] }, fasta, gff, true, params.gff )
        }
    }

    /*
     * Annotate variants
     */
    if (params.gff && !params.skip_variants_snpeff) {
        SNPEFF_SNPSIFT ( VCF_BGZIP_TABIX_STATS.out.vcf, snpeff_db, snpeff_config, fasta )
    }

    emit:
    tsv              = IVAR_VARIANTS.out.tsv              // channel: [ val(meta), [ tsv ] ]
    ivar_version     = IVAR_VARIANTS.out.version          //    path: *.version.txt

    vcf_orig         = IVAR_VARIANTS_TO_VCF.out.vcf       // channel: [ val(meta), [ vcf ] ]
    log_out          = IVAR_VARIANTS_TO_VCF.out.log       // channel: [ val(meta), [ log ] ]
    multiqc_tsv      = IVAR_VARIANTS_TO_VCF.out.tsv       // channel: [ val(meta), [ tsv ] ]

    vcf              = VCF_BGZIP_TABIX_STATS.out.vcf      // channel: [ val(meta), [ vcf ] ]
    tbi              = VCF_BGZIP_TABIX_STATS.out.tbi      // channel: [ val(meta), [ tbi ] ]
    stats            = VCF_BGZIP_TABIX_STATS.out.stats    // channel: [ val(meta), [ txt ] ]
    bcftools_version = VCF_BGZIP_TABIX_STATS.out.version  //    path: *.version.txt
    
    consensus        = IVAR_CONSENSUS.out.fasta           // channel: [ val(meta), [ fasta ] ]
    consensus_qual   = IVAR_CONSENSUS.out.txt             // channel: [ val(meta), [ txt ] ]
    
    bases_tsv        = PLOT_BASE_DENSITY.out.tsv          // channel: [ val(meta), [ tsv ] ]
    bases_pdf        = PLOT_BASE_DENSITY.out.pdf          // channel: [ val(meta), [ pdf ] ]
    
    quast_results    = QUAST.out.results                  // channel: [ val(meta), [ results ] ]
    quast_tsv        = QUAST.out.tsv                      // channel: [ val(meta), [ tsv ] ]
    quast_version    = QUAST.out.version                  //    path: *.version.txt

    snpeff_vcf       = SNPEFF_SNPSIFT.out.vcf             // channel: [ val(meta), [ vcf.gz ] ]
    snpeff_tbi       = SNPEFF_SNPSIFT.out.tbi             // channel: [ val(meta), [ tbi ] ]
    snpeff_stats     = SNPEFF_SNPSIFT.out.stats           // channel: [ val(meta), [ txt ] ]
    snpeff_csv       = SNPEFF_SNPSIFT.out.csv             // channel: [ val(meta), [ csv ] ]
    snpeff_txt       = SNPEFF_SNPSIFT.out.txt             // channel: [ val(meta), [ txt ] ]
    snpeff_html      = SNPEFF_SNPSIFT.out.html            // channel: [ val(meta), [ html ] ]
    snpsift_txt      = SNPEFF_SNPSIFT.out.snpsift_txt     // channel: [ val(meta), [ txt ] ]
    snpeff_version   = SNPEFF_SNPSIFT.out.snpeff_version  //    path: *.version.txt
    snpsift_version  = SNPEFF_SNPSIFT.out.snpsift_version //    path: *.version.txt
}

