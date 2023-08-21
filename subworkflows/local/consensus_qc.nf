//
// Consensus calling QC
//

include { QUAST             } from '../../modules/nf-core/quast/main'
include { PANGOLIN          } from '../../modules/nf-core/pangolin/main'
include { NEXTCLADE_RUN     } from '../../modules/nf-core/nextclade/run/main'
include { PLOT_BASE_DENSITY } from '../../modules/local/plot_base_density'

workflow CONSENSUS_QC {
    take:
    consensus    // channel: [ val(meta), [ consensus ] ]
    fasta        // channel: /path/to/genome.fasta
    gff          // channel: /path/to/genome.gff
    nextclade_db // channel: /path/to/nextclade_db/

    main:

    ch_versions = Channel.empty()

    //
    // Consensus QC report across samples with QUAST
    //
    ch_quast_results = Channel.empty()
    ch_quast_tsv     = Channel.empty()
    if (!params.skip_variants_quast) {
        QUAST (
            consensus.collect{ it[1] },
            fasta,
            gff,
            true,
            params.gff
        )
        ch_quast_results = QUAST.out.results
        ch_quast_tsv     = QUAST.out.tsv
        ch_versions      = ch_versions.mix(QUAST.out.versions)
    }

    //
    // Lineage analysis with Pangolin
    //
    ch_pangolin_report = Channel.empty()
    if (!params.skip_pangolin) {
        PANGOLIN (
            consensus
        )
        ch_pangolin_report = PANGOLIN.out.report
        ch_versions        = ch_versions.mix(PANGOLIN.out.versions.first())
    }

    //
    // Lineage analysis with Nextclade
    //
    ch_nextclade_report = Channel.empty()
    if (!params.skip_nextclade) {
        NEXTCLADE_RUN (
            consensus,
            nextclade_db
        )
        ch_nextclade_report = NEXTCLADE_RUN.out.csv
        ch_versions         = ch_versions.mix(NEXTCLADE_RUN.out.versions.first())
    }

    //
    // Plot consensus base density
    //
    ch_bases_tsv = Channel.empty()
    ch_bases_pdf = Channel.empty()
    if (!params.skip_consensus_plots) {
        PLOT_BASE_DENSITY (
            consensus
        )
        ch_bases_tsv = PLOT_BASE_DENSITY.out.tsv
        ch_bases_pdf = PLOT_BASE_DENSITY.out.pdf
        ch_versions  = ch_versions.mix(PLOT_BASE_DENSITY.out.versions.first())
    }

    emit:
    quast_results    = ch_quast_results    // channel: [ val(meta), [ results ] ]
    quast_tsv        = ch_quast_tsv        // channel: [ val(meta), [ tsv ] ]

    pangolin_report  = ch_pangolin_report  // channel: [ val(meta), [ csv ] ]

    nextclade_report = ch_nextclade_report // channel: [ val(meta), [ csv ] ]

    bases_tsv        = ch_bases_tsv        // channel: [ val(meta), [ tsv ] ]
    bases_pdf        = ch_bases_pdf        // channel: [ val(meta), [ pdf ] ]

    versions         = ch_versions         // channel: [ versions.yml ]
}
