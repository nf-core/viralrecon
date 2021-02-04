/*
 * Variant calling and downstream processing for IVar
 */

params.ivar_variants_options                 = [:]
params.ivar_consensus_options                = [:]
params.quast_options                         = [:]
params.ivar_variants_to_vcf_lowfreq_options  = [:]
params.ivar_variants_to_vcf_highfreq_options = [:]
params.ivar_bgzip_lowfreq_options            = [:]
params.ivar_tabix_lowfreq_options            = [:]
params.ivar_stats_lowfreq_options            = [:]
params.ivar_bgzip_highfreq_options           = [:]
params.ivar_tabix_highfreq_options           = [:]
params.ivar_stats_highfreq_options           = [:]
params.snpeff_lowfreq_options                = [:]
params.snpsift_lowfreq_options               = [:]
params.snpeff_lowfreq_bgzip_options          = [:]
params.snpeff_lowfreq_tabix_options          = [:]
params.snpeff_lowfreq_stats_options          = [:]
params.snpeff_highfreq_options               = [:]
params.snpsift_highfreq_options              = [:]
params.snpeff_highfreq_bgzip_options         = [:]
params.snpeff_highfreq_tabix_options         = [:]
params.snpeff_highfreq_stats_options         = [:]

include { IVAR_VARIANTS                                           } from '../../modules/local/ivar_variants'        addParams( options: params.ivar_variants_options                 )
include { IVAR_CONSENSUS                                          } from '../../modules/local/ivar_consensus'       addParams( options: params.ivar_consensus_options                )
include { IVAR_VARIANTS_TO_VCF as IVAR_VARIANTS_TO_VCF_LOWFREQ    } from '../../modules/local/ivar_variants_to_vcf' addParams( options: params.ivar_variants_to_vcf_lowfreq_options  )
include { IVAR_VARIANTS_TO_VCF as IVAR_VARIANTS_TO_VCF_HIGHFREQ   } from '../../modules/local/ivar_variants_to_vcf' addParams( options: params.ivar_variants_to_vcf_highfreq_options )
include { VCF_BGZIP_TABIX_STATS as VCF_BGZIP_TABIX_STATS_LOWFREQ  } from './vcf_bgzip_tabix_stats'                  addParams( bgzip_options: params.ivar_bgzip_lowfreq_options, tabix_options: params.ivar_tabix_lowfreq_options, stats_options: params.ivar_stats_lowfreq_options    )
include { VCF_BGZIP_TABIX_STATS as VCF_BGZIP_TABIX_STATS_HIGHFREQ } from './vcf_bgzip_tabix_stats'                  addParams( bgzip_options: params.ivar_bgzip_highfreq_options, tabix_options: params.ivar_tabix_highfreq_options, stats_options: params.ivar_stats_highfreq_options )
include { SNPEFF_SNPSIFT as SNPEFF_SNPSIFT_LOWFREQ                } from './snpeff_snpsift'                         addParams( snpeff_options: params.snpeff_lowfreq_options, snpsift_options: params.snpsift_lowfreq_options, bgzip_options: params.snpeff_lowfreq_bgzip_options, tabix_options: params.snpeff_lowfreq_tabix_options, stats_options: params.snpeff_lowfreq_stats_options      )
include { SNPEFF_SNPSIFT as SNPEFF_SNPSIFT_HIGHFREQ               } from './snpeff_snpsift'                         addParams( snpeff_options: params.snpeff_highfreq_options, snpsift_options: params.snpsift_highfreq_options, bgzip_options: params.snpeff_highfreq_bgzip_options, tabix_options: params.snpeff_highfreq_tabix_options, stats_options: params.snpeff_highfreq_stats_options )
include { QUAST                                                   } from '../../modules/nf-core/quast/main'         addParams( options: params.quast_options )

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
     * Convert original ivar output to VCF, zip and index
     */
    IVAR_VARIANTS_TO_VCF_LOWFREQ ( IVAR_VARIANTS.out.tsv, ivar_multiqc_header )

    VCF_BGZIP_TABIX_STATS_LOWFREQ ( IVAR_VARIANTS_TO_VCF_LOWFREQ.out.vcf )

    /*
     * Convert original ivar output to VCF, filter for high frequency variants, zip and index
     */
    IVAR_VARIANTS_TO_VCF_HIGHFREQ ( IVAR_VARIANTS.out.tsv, ivar_multiqc_header )

    VCF_BGZIP_TABIX_STATS_HIGHFREQ ( IVAR_VARIANTS_TO_VCF_HIGHFREQ.out.vcf )

    /*
     * Create genome consensus
     */
    if (!params.skip_consensus) {
        IVAR_CONSENSUS ( mpileup )

        if (!params.skip_variants_quast) {
            QUAST ( IVAR_CONSENSUS.out.fasta.collect{ it[1] }, fasta, gff, true, if params.gff )
        }
    }

    /*
     * Annotate low/high frequency variants
     */
    if (params.gff && !params.skip_variants_snpeff) {
        SNPEFF_SNPSIFT_LOWFREQ ( VCF_BGZIP_TABIX_STATS_LOWFREQ.out.vcf, snpeff_db, snpeff_config, fasta )

        SNPEFF_SNPSIFT_HIGHFREQ ( VCF_BGZIP_TABIX_STATS_HIGHFREQ.out.vcf, snpeff_db, snpeff_config, fasta )
    }

    // emit:
    // scaffolds         = SPADES.out.scaffolds              // channel: [ val(meta), [ scaffolds ] ]
    // spades_version    = SPADES.out.version                //    path: *.version.txt

}
