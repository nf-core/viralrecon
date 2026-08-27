//
// Downstream analysis for assembly scaffolds
//

include { FILTER_BLASTN } from '../../../modules/local/filter_blastn'
include { ABACAS        } from '../../../modules/nf-core/abacas/main'
include { BLAST_BLASTN  } from '../../../modules/nf-core/blast/blastn/main'
include { PLASMIDID     } from '../../../modules/nf-core/plasmidid/main'
include { QUAST         } from '../../../modules/nf-core/quast/main'

workflow ASSEMBLY_QC {
    take:
    scaffolds             // channel: [ val(meta), [ scaffolds ] ]
    fasta                 // channel: /path/to/genome.fasta
    gff                   // channel: /path/to/genome.gff
    blast_db              // channel: /path/to/blast_db/
    blast_header          // channel: /path/to/blast_header.txt
    blast_filtered_header // channel: /path/to/blast_filtered_header.txt
    ch_taxidlist          // channel: /path/to/taxidlist.txt

    main:

    ch_versions = channel.empty()

    //
    // Run blastn on assembly scaffolds
    //
    ch_blast_txt        = channel.empty()
    ch_blast_filter_txt = channel.empty()
    if (!params.skip_blast) {
        BLAST_BLASTN (
            scaffolds,
            blast_db,
            ch_taxidlist,
            [],
            []
        )

        FILTER_BLASTN (
            BLAST_BLASTN.out.txt,
            blast_header,
            blast_filtered_header
        )
        ch_blast_txt = FILTER_BLASTN.out.blast
        ch_blast_filter_txt = FILTER_BLASTN.out.txt
    }

    //
    // Assembly QC across all samples with QUAST
    //
    ch_quast_results = channel.empty()
    ch_quast_tsv     = channel.empty()
    if (!params.skip_assembly_quast) {
        scaffolds
            .collect{ _meta, scaffolds_file -> scaffolds_file }
            .map { scaffolds_collect -> tuple([id: "quast"], scaffolds_collect) }
            .set { ch_to_quast }

        QUAST (
            ch_to_quast,
            fasta.map { fasta_files ->
                [ [:], fasta_files ] },
            gff
        )
        ch_quast_results = QUAST.out.results
        ch_quast_tsv     = QUAST.out.tsv
    }

    //
    // Contiguate assembly with ABACAS
    //
    ch_abacas_results = channel.empty()
    if (!params.skip_abacas) {
        ABACAS (
            scaffolds,
            fasta.map { fasta_files ->
                [ [:], fasta_files ] }
        )
        ch_abacas_results = ABACAS.out.results
    }

    //
    // Assembly report with PlasmidID
    //
    ch_plasmidid_html     = channel.empty()
    ch_plasmidid_tab      = channel.empty()
    ch_plasmidid_images   = channel.empty()
    ch_plasmidid_logs     = channel.empty()
    ch_plasmidid_data     = channel.empty()
    ch_plasmidid_database = channel.empty()
    ch_plasmidid_fasta    = channel.empty()
    ch_plasmidid_kmer     = channel.empty()
    ch_versions           = channel.empty()
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
        ch_versions           = PLASMIDID.out.versions
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

    versions           = ch_versions           // channel: versions.yml
}
