//
// Variant calling and downstream processing for IVar
//

include { IVAR_VARIANTS_TO_VCF  } from '../../modules/local/ivar_variants_to_vcf'
include { PLOT_BASE_DENSITY     } from '../../modules/local/plot_base_density'
include { IVAR_VARIANTS         } from '../../modules/nf-core/modules/ivar/variants/main'
include { IVAR_CONSENSUS        } from '../../modules/nf-core/modules/ivar/consensus/main'
include { QUAST                 } from '../../modules/nf-core/modules/quast/main'
include { PANGOLIN              } from '../../modules/nf-core/modules/pangolin/main'
include { NEXTCLADE_RUN         } from '../../modules/nf-core/modules/nextclade/run/main'
include { ASCIIGENOME           } from '../../modules/local/asciigenome'

include { VCF_BGZIP_TABIX_STATS } from '../nf-core/vcf_bgzip_tabix_stats'
include { SNPEFF_SNPSIFT        } from './snpeff_snpsift'

workflow VARIANTS_IVAR {
    take:
    bam                 // channel: [ val(meta), [ bam ] ]
    fasta               // channel: /path/to/genome.fasta
    sizes               // channel: /path/to/genome.sizes
    gff                 // channel: /path/to/genome.gff
    bed                 // channel: /path/to/primers.bed
    nextclade_db        // channel: /path/to/nextclade_db/
    snpeff_db           // channel: /path/to/snpeff_db/
    snpeff_config       // channel: /path/to/snpeff.config
    ivar_multiqc_header // channel: /path/to/multiqc_header for ivar variants

    main:

    ch_versions = Channel.empty()

    //
    // Call variants
    //
    IVAR_VARIANTS (
        bam,
        fasta,
        gff,
        params.save_mpileup
    )
    ch_versions = ch_versions.mix(IVAR_VARIANTS.out.versions.first())

    //
    // Convert original iVar output to VCF, zip and index
    //
    IVAR_VARIANTS_TO_VCF (
        IVAR_VARIANTS.out.tsv,
        ivar_multiqc_header
    )
    ch_versions = ch_versions.mix(IVAR_VARIANTS_TO_VCF.out.versions.first())

    VCF_BGZIP_TABIX_STATS (
        IVAR_VARIANTS_TO_VCF.out.vcf
    )
    ch_versions = ch_versions.mix(VCF_BGZIP_TABIX_STATS.out.versions)

    //
    // Create genome consensus
    //
    ch_consensus        = Channel.empty()
    ch_consensus_qual   = Channel.empty()
    ch_bases_tsv        = Channel.empty()
    ch_bases_pdf        = Channel.empty()
    ch_quast_results    = Channel.empty()
    ch_quast_tsv        = Channel.empty()
    ch_pangolin_report  = Channel.empty()
    ch_nextclade_report = Channel.empty()
    if (!params.skip_consensus) {
        IVAR_CONSENSUS (
            bam,
            fasta,
            params.save_mpileup
        )
        ch_consensus      = IVAR_CONSENSUS.out.fasta
        ch_consensus_qual = IVAR_CONSENSUS.out.qual
        ch_versions       = ch_versions.mix(IVAR_CONSENSUS.out.versions.first())

        PLOT_BASE_DENSITY (
            ch_consensus
        )
        ch_bases_tsv = PLOT_BASE_DENSITY.out.tsv
        ch_bases_pdf = PLOT_BASE_DENSITY.out.pdf
        ch_versions  = ch_versions.mix(PLOT_BASE_DENSITY.out.versions.first())

        if (!params.skip_variants_quast) {
            QUAST (
                ch_consensus.collect{ it[1] },
                fasta,
                gff,
                true,
                params.gff
            )
            ch_quast_results = QUAST.out.results
            ch_quast_tsv     = QUAST.out.tsv
            ch_versions      = ch_versions.mix(QUAST.out.versions)
        }

        if (!params.skip_pangolin) {
            PANGOLIN (
                ch_consensus
            )
            ch_pangolin_report = PANGOLIN.out.report
            ch_versions        = ch_versions.mix(PANGOLIN.out.versions.first())
        }

        if (!params.skip_nextclade) {
            NEXTCLADE_RUN (
                ch_consensus,
                nextclade_db
            )
            ch_nextclade_report = NEXTCLADE_RUN.out.csv
            ch_versions         = ch_versions.mix(NEXTCLADE_RUN.out.versions.first())
        }
    }

    //
    // Annotate variants
    //
    ch_snpeff_vcf   = Channel.empty()
    ch_snpeff_tbi   = Channel.empty()
    ch_snpeff_stats = Channel.empty()
    ch_snpeff_csv   = Channel.empty()
    ch_snpeff_txt   = Channel.empty()
    ch_snpeff_html  = Channel.empty()
    ch_snpsift_txt  = Channel.empty()
    if (params.gff && !params.skip_snpeff) {
        SNPEFF_SNPSIFT (
            VCF_BGZIP_TABIX_STATS.out.vcf,
            snpeff_db,
            snpeff_config,
            fasta
        )
        ch_snpeff_vcf   = SNPEFF_SNPSIFT.out.vcf
        ch_snpeff_tbi   = SNPEFF_SNPSIFT.out.tbi
        ch_snpeff_stats = SNPEFF_SNPSIFT.out.stats
        ch_snpeff_csv   = SNPEFF_SNPSIFT.out.csv
        ch_snpeff_txt   = SNPEFF_SNPSIFT.out.txt
        ch_snpeff_html  = SNPEFF_SNPSIFT.out.html
        ch_snpsift_txt  = SNPEFF_SNPSIFT.out.snpsift_txt
        ch_versions     = ch_versions.mix(SNPEFF_SNPSIFT.out.versions)
    }

    //
    // Variant screenshots with ASCIIGenome
    //
    ch_asciigenome_pdf = Channel.empty()
    if (!params.skip_asciigenome) {
        bam
            .join(VCF_BGZIP_TABIX_STATS.out.vcf, by: [0])
            .join(VCF_BGZIP_TABIX_STATS.out.stats, by: [0])
            .map { meta, bam, vcf, stats ->
                if (WorkflowCommons.getNumVariantsFromBCFToolsStats(stats) > 0) {
                    return [ meta, bam, vcf ]
                }
            }
            .set { ch_asciigenome }

        ASCIIGENOME (
            ch_asciigenome,
            fasta,
            sizes,
            gff,
            bed,
            params.asciigenome_window_size,
            params.asciigenome_read_depth
        )
        ch_asciigenome_pdf = ASCIIGENOME.out.pdf
        ch_versions        = ch_versions.mix(ASCIIGENOME.out.versions.first())
    }

    emit:
    tsv              = IVAR_VARIANTS.out.tsv           // channel: [ val(meta), [ tsv ] ]

    vcf_orig         = IVAR_VARIANTS_TO_VCF.out.vcf    // channel: [ val(meta), [ vcf ] ]
    log_out          = IVAR_VARIANTS_TO_VCF.out.log    // channel: [ val(meta), [ log ] ]
    multiqc_tsv      = IVAR_VARIANTS_TO_VCF.out.tsv    // channel: [ val(meta), [ tsv ] ]

    vcf              = VCF_BGZIP_TABIX_STATS.out.vcf   // channel: [ val(meta), [ vcf ] ]
    tbi              = VCF_BGZIP_TABIX_STATS.out.tbi   // channel: [ val(meta), [ tbi ] ]
    stats            = VCF_BGZIP_TABIX_STATS.out.stats // channel: [ val(meta), [ txt ] ]

    consensus        = ch_consensus                    // channel: [ val(meta), [ fasta ] ]
    consensus_qual   = ch_consensus_qual               // channel: [ val(meta), [ fasta ] ]
    bases_tsv        = ch_bases_tsv                    // channel: [ val(meta), [ tsv ] ]
    bases_pdf        = ch_bases_pdf                    // channel: [ val(meta), [ pdf ] ]

    quast_results    = ch_quast_results                // channel: [ val(meta), [ results ] ]
    quast_tsv        = ch_quast_tsv                    // channel: [ val(meta), [ tsv ] ]

    snpeff_vcf       = ch_snpeff_vcf                   // channel: [ val(meta), [ vcf.gz ] ]
    snpeff_tbi       = ch_snpeff_tbi                   // channel: [ val(meta), [ tbi ] ]
    snpeff_stats     = ch_snpeff_stats                 // channel: [ val(meta), [ txt ] ]
    snpeff_csv       = ch_snpeff_csv                   // channel: [ val(meta), [ csv ] ]
    snpeff_txt       = ch_snpeff_txt                   // channel: [ val(meta), [ txt ] ]
    snpeff_html      = ch_snpeff_html                  // channel: [ val(meta), [ html ] ]
    snpsift_txt      = ch_snpsift_txt                  // channel: [ val(meta), [ txt ] ]

    pangolin_report  = ch_pangolin_report              // channel: [ val(meta), [ csv ] ]

    nextclade_report = ch_nextclade_report             // channel: [ val(meta), [ csv ] ]

    asciigenome_pdf  = ch_asciigenome_pdf              // channel: [ val(meta), [ pdf ] ]

    versions         = ch_versions                     // channel: [ versions.yml ]
}
