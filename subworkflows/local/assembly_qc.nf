//
// Downstream analysis for assembly scaffolds
//

include { FILTER_BLASTN } from '../../modules/local/filter_blastn'
include { ABACAS        } from '../../modules/nf-core/abacas/main'
include { BLAST_BLASTN  } from '../../modules/nf-core/blast/blastn/main'
include { PLASMIDID     } from '../../modules/nf-core/plasmidid/main'
include { QUAST         } from '../../modules/nf-core/quast/main'

workflow ASSEMBLY_QC {
    take:
    scaffolds    // channel: [ val(meta), [ scaffolds ] ]
    fasta        // channel: /path/to/genome.fasta
    gff          // channel: /path/to/genome.gff
    blast_db     // channel: /path/to/blast_db/
    blast_header // channel: /path/to/blast_header.txt

    main:

    ch_versions = Channel.empty()

    //
    // Run blastn on assembly scaffolds
    //
    ch_blast_txt        = Channel.empty()
    ch_blast_filter_txt = Channel.empty()
    if (!params.skip_blast) {
        BLAST_BLASTN (
            scaffolds,
            blast_db
        )
        ch_blast_txt = BLAST_BLASTN.out.txt
        ch_versions  = ch_versions.mix(BLAST_BLASTN.out.versions.first())

        FILTER_BLASTN (
            BLAST_BLASTN.out.txt,
            blast_header
        )
        ch_blast_filter_txt = FILTER_BLASTN.out.txt
        ch_versions         = ch_versions.mix(FILTER_BLASTN.out.versions.first())
    }

    //
    // Assembly QC across all samples with QUAST
    //
    ch_quast_results = Channel.empty()
    ch_quast_tsv     = Channel.empty()
    if (!params.skip_assembly_quast) {
        QUAST (
            scaffolds.collect{ it[1] },
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
    // Contiguate assembly with ABACAS
    //
    ch_abacas_results = Channel.empty()
    if (!params.skip_abacas) {
        ABACAS (
            scaffolds,
            fasta
        )
        ch_abacas_results = ABACAS.out.results
        ch_versions       = ch_versions.mix(ABACAS.out.versions.first())
    }

    //
    // Assembly report with PlasmidID
    //
    ch_plasmidid_html     = Channel.empty()
    ch_plasmidid_tab      = Channel.empty()
    ch_plasmidid_images   = Channel.empty()
    ch_plasmidid_logs     = Channel.empty()
    ch_plasmidid_data     = Channel.empty()
    ch_plasmidid_database = Channel.empty()
    ch_plasmidid_fasta    = Channel.empty()
    ch_plasmidid_kmer     = Channel.empty()
    if (!params.skip_plasmidid) {
        PLASMIDID (
            scaffolds,
            fasta
        )
        ch_plasmidid_html     = PLASMIDID.out.html
        ch_plasmidid_tab      = PLASMIDID.out.tab
        ch_plasmidid_images   = PLASMIDID.out.images
        ch_plasmidid_logs     = PLASMIDID.out.logs
        ch_plasmidid_data     = PLASMIDID.out.data
        ch_plasmidid_database = PLASMIDID.out.database
        ch_plasmidid_fasta    = PLASMIDID.out.fasta_files
        ch_plasmidid_kmer     = PLASMIDID.out.kmer
        ch_versions           = ch_versions.mix(PLASMIDID.out.versions.first())
    }

    emit:
    blast_txt          = ch_blast_txt          // channel: [ val(meta), [ txt ] ]
    blast_filter_txt   = ch_blast_filter_txt   // channel: [ val(meta), [ txt ] ]

    quast_results      = ch_quast_results      // channel: [ val(meta), [ results ] ]
    quast_tsv          = ch_quast_tsv          // channel: [ val(meta), [ tsv ] ]

    abacas_results     = ch_abacas_results     // channel: [ val(meta), [ results ] ]

    plasmidid_html     = ch_plasmidid_html     // channel: [ val(meta), [ html ] ]
    plasmidid_tab      = ch_plasmidid_tab      // channel: [ val(meta), [ tab ] ]
    plasmidid_images   = ch_plasmidid_images   // channel: [ val(meta), [ images/ ] ]
    plasmidid_logs     = ch_plasmidid_logs     // channel: [ val(meta), [ logs/ ] ]
    plasmidid_data     = ch_plasmidid_data     // channel: [ val(meta), [ data/ ] ]
    plasmidid_database = ch_plasmidid_database // channel: [ val(meta), [ database/ ] ]
    plasmidid_fasta    = ch_plasmidid_fasta    // channel: [ val(meta), [ fasta_files/ ] ]
    plasmidid_kmer     = ch_plasmidid_kmer     // channel: [ val(meta), [ kmer/ ] ]

    versions           = ch_versions           // channel: [ versions.yml ]
}
