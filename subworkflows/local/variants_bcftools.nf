//
// Variant calling and downstream processing for BCFTools
//

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
params.pangolin_options            = [:]
params.nextclade_options           = [:]
params.asciigenome_options         = [:]

include { BCFTOOLS_MPILEUP } from '../../modules/nf-core/software/bcftools/mpileup/main' addParams( options: params.bcftools_mpileup_options )
include { QUAST            } from '../../modules/nf-core/software/quast/main'            addParams( options: params.quast_options            )
include { PANGOLIN         } from '../../modules/nf-core/software/pangolin/main'         addParams( options: params.pangolin_options         )
include { NEXTCLADE        } from '../../modules/nf-core/software/nextclade/main'        addParams( options: params.nextclade_options        )
include { ASCIIGENOME      } from '../../modules/local/asciigenome'                      addParams( options: params.asciigenome_options      )
include { MAKE_CONSENSUS   } from './make_consensus'                                     addParams( genomecov_options: params.consensus_genomecov_options, merge_options: params.consensus_merge_options, mask_options: params.consensus_mask_options, maskfasta_options: params.consensus_maskfasta_options, bcftools_options: params.consensus_bcftools_options, plot_bases_options: params.consensus_plot_options )
include { SNPEFF_SNPSIFT   } from './snpeff_snpsift'                                     addParams( snpeff_options: params.snpeff_options, snpsift_options: params.snpsift_options, bgzip_options: params.snpeff_bgzip_options, tabix_options: params.snpeff_tabix_options, stats_options:  params.snpeff_stats_options )

workflow VARIANTS_BCFTOOLS {
    take:
    bam           // channel: [ val(meta), [ bam ] ]
    fasta         // channel: /path/to/genome.fasta
    gff           // channel: /path/to/genome.gff
    bed           // channel: /path/to/primers.bed
    snpeff_db     // channel: /path/to/snpeff_db/
    snpeff_config // channel: /path/to/snpeff.config

    main:

    //
    // Call variants
    //
    BCFTOOLS_MPILEUP ( bam, fasta )

    //
    // Create genome consensus using variants in VCF, run QUAST and pangolin
    //
    ch_consensus         = Channel.empty()
    ch_bases_tsv         = Channel.empty()
    ch_bases_pdf         = Channel.empty()
    ch_bedtools_version  = Channel.empty()
    ch_quast_results     = Channel.empty()
    ch_quast_tsv         = Channel.empty()
    ch_quast_version     = Channel.empty()
    ch_pangolin_report   = Channel.empty()
    ch_pangolin_version  = Channel.empty()
    ch_nextclade_report  = Channel.empty()
    ch_nextclade_version = Channel.empty()
    if (!params.skip_consensus) {
        MAKE_CONSENSUS ( bam.join(BCFTOOLS_MPILEUP.out.vcf, by: [0]).join(BCFTOOLS_MPILEUP.out.tbi, by: [0]), fasta )
        ch_consensus        = MAKE_CONSENSUS.out.fasta
        ch_bases_tsv        = MAKE_CONSENSUS.out.tsv
        ch_bases_pdf        = MAKE_CONSENSUS.out.pdf
        ch_bedtools_version = MAKE_CONSENSUS.out.bedtools_version

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

        if (!params.skip_nextclade) {
            NEXTCLADE ( ch_consensus, 'csv' )
            ch_nextclade_report  = NEXTCLADE.out.csv
            ch_nextclade_version = NEXTCLADE.out.version
        }
    }

    //
    // Annotate variants
    //
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
        SNPEFF_SNPSIFT ( BCFTOOLS_MPILEUP.out.vcf, snpeff_db, snpeff_config, fasta )
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

    //
    // Variant screenshots with ASCIIGenome
    //
    ch_asciigenome_pdf     = Channel.empty()
    ch_asciigenome_version = Channel.empty()
    if (!params.skip_asciigenome) {
        bam
            .join(BCFTOOLS_MPILEUP.out.vcf, by: [0])
            .join(BCFTOOLS_MPILEUP.out.stats, by: [0])
            .map { meta, bam, vcf, stats ->
                if (WorkflowCommons.getNumVariantsFromBCFToolsStats(stats) > 0) {
                    return [ meta, bam, vcf ]
                }
            }
            .set { ch_asciigenome }

        ASCIIGENOME (
            ch_asciigenome,
            fasta,
            gff,
            bed,
            params.asciigenome_window_size,
            params.asciigenome_read_depth
        )
        ch_asciigenome_pdf     = ASCIIGENOME.out.pdf
        ch_asciigenome_version = ASCIIGENOME.out.version
    }

    emit:
    vcf                 = BCFTOOLS_MPILEUP.out.vcf     // channel: [ val(meta), [ vcf ] ]
    tbi                 = BCFTOOLS_MPILEUP.out.tbi     // channel: [ val(meta), [ tbi ] ]
    stats               = BCFTOOLS_MPILEUP.out.stats   // channel: [ val(meta), [ txt ] ]
    bcftools_version    = BCFTOOLS_MPILEUP.out.version //    path: *.version.txt

    consensus           = ch_consensus                 // channel: [ val(meta), [ fasta ] ]
    bases_tsv           = ch_bases_tsv                 // channel: [ val(meta), [ tsv ] ]
    bases_pdf           = ch_bases_pdf                 // channel: [ val(meta), [ pdf ] ]
    bedtools_version    = ch_bedtools_version          //    path: *.version.txt

    quast_results       = ch_quast_results             // channel: [ val(meta), [ results ] ]
    quast_tsv           = ch_quast_tsv                 // channel: [ val(meta), [ tsv ] ]
    quast_version       = ch_quast_version             //    path: *.version.txt

    snpeff_vcf          = ch_snpeff_vcf                // channel: [ val(meta), [ vcf.gz ] ]
    snpeff_tbi          = ch_snpeff_tbi                // channel: [ val(meta), [ tbi ] ]
    snpeff_stats        = ch_snpeff_stats              // channel: [ val(meta), [ txt ] ]
    snpeff_csv          = ch_snpeff_csv                // channel: [ val(meta), [ csv ] ]
    snpeff_txt          = ch_snpeff_txt                // channel: [ val(meta), [ txt ] ]
    snpeff_html         = ch_snpeff_html               // channel: [ val(meta), [ html ] ]
    snpsift_txt         = ch_snpsift_txt               // channel: [ val(meta), [ txt ] ]
    snpeff_version      = ch_snpeff_version            //    path: *.version.txt
    snpsift_version     = ch_snpsift_version           //    path: *.version.txt

    pangolin_report     = ch_pangolin_report           // channel: [ val(meta), [ csv ] ]
    pangolin_version    = ch_pangolin_version          //    path: *.version.txt

    nextclade_report    = ch_nextclade_report          // channel: [ val(meta), [ csv ] ]
    nextclade_version   = ch_nextclade_version         //    path: *.version.txt

    asciigenome_pdf     = ch_asciigenome_pdf           // channel: [ val(meta), [ pdf ] ]
    asciigenome_version = ch_asciigenome_version       //    path: *.version.txt
}
