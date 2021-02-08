/*
 * Variant calling and downstream processing for BCFTools
 */

params.bcftools_mpileup_options    = [:]
params.quast_options               = [:]
params.consensus_genomecov_options = [:]
params.consensus_merge_options     = [:]
params.consensus_mask_options      = [:]
params.consensus_maskfasta_options = [:]
params.consensus_bcftools_options  = [:]
params.consensus_plot_options      = [:]
params.snpeff_options              = [:]
params.snpsift_options             = [:]
params.snpeff_bgzip_options        = [:]
params.snpeff_tabix_options        = [:]
params.snpeff_stats_options        = [:]

include { BCFTOOLS_MPILEUP } from '../../modules/local/bcftools_mpileup'      addParams( options: params.bcftools_mpileup_options ) 
include { QUAST            } from '../../modules/nf-core/software/quast/main' addParams( options: params.quast_options            )
include { MAKE_CONSENSUS   } from './make_consensus'                          addParams( genomecov_options: params.consensus_genomecov_options, merge_options: params.consensus_merge_options, mask_options: params.consensus_mask_options, maskfasta_options: params.consensus_maskfasta_options, bcftools_options: params.consensus_bcftools_options, plot_bases_options: params.consensus_plot_options )
include { SNPEFF_SNPSIFT   } from './snpeff_snpsift'                          addParams( snpeff_options: params.snpeff_options, snpsift_options: params.snpsift_options, bgzip_options: params.snpeff_bgzip_options, tabix_options: params.snpeff_tabix_options, stats_options:  params.snpeff_stats_options )

workflow VARIANTS_BCFTOOLS {
    take:
    bam           // channel: [ val(meta), [ bam ] ]
    fasta         // channel: /path/to/genome.fasta
    gff           // channel: /path/to/genome.gff
    snpeff_db     // channel: /path/to/snpeff_db/
    snpeff_config // channel: /path/to/snpeff.config
    
    main:
    /*
     * Call variants
     */
    BCFTOOLS_MPILEUP ( bam, fasta )

    /*
     * Create genome consensus using variants in VCF
     */
    if (!params.skip_consensus) {
        MAKE_CONSENSUS ( bam.join(BCFTOOLS_MPILEUP.out.vcf, by: [0]).join(BCFTOOLS_MPILEUP.out.tbi, by: [0]), fasta )

        if (!params.skip_variants_quast) {
            QUAST ( MAKE_CONSENSUS.out.fasta.collect{ it[1] }, fasta, gff, true, params.gff )
        }
    }

    /*
     * Annotate variants
     */
    if (params.gff && !params.skip_variants_snpeff) {
        SNPEFF_SNPSIFT ( BCFTOOLS_MPILEUP.out.vcf, snpeff_db, snpeff_config, fasta )
    }

    emit:
    vcf              = BCFTOOLS_MPILEUP.out.vcf            // channel: [ val(meta), [ vcf ] ]
    tbi              = BCFTOOLS_MPILEUP.out.tbi            // channel: [ val(meta), [ tbi ] ]
    stats            = BCFTOOLS_MPILEUP.out.txt            // channel: [ val(meta), [ txt ] ]
    bcftools_version = BCFTOOLS_MPILEUP.out.version        //    path: *.version.txt
    
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