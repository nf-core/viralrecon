//
// Assembly and downstream processing for SPAdes scaffolds
//

include { SPADES                     } from '../../modules/nf-core/spades/main'
include { BANDAGE_IMAGE              } from '../../modules/nf-core/bandage/image/main'
include { GUNZIP as GUNZIP_SCAFFOLDS } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFA       } from '../../modules/nf-core/gunzip/main'

include { ASSEMBLY_QC   } from './assembly_qc'

workflow ASSEMBLY_SPADES {
    take:
    reads        // channel: [ val(meta), [ reads ] ]
    mode         // string : spades assembly mode e.g. 'rnaviral'
    hmm          // channel: /path/to/spades.hmm
    fasta        // channel: /path/to/genome.fasta
    gff          // channel: /path/to/genome.gff
    blast_db     // channel: /path/to/blast_db/
    blast_header // channel: /path/to/blast_header.txt

    main:

    ch_versions = Channel.empty()

    //
    // Filter for paired-end samples if running metaSPAdes / metaviralSPAdes / metaplasmidSPAdes
    //
    ch_reads = reads
    if (mode.contains('meta') || mode.contains('bio')) {
        reads
            .filter { meta, illumina, pacbio, nanopore -> !meta.single_end }
            .set { ch_reads }
    }

    //
    // Assemble reads with SPAdes
    //
    SPADES (
        ch_reads,
        [],
        hmm
    )
    ch_versions = ch_versions.mix(SPADES.out.versions.first())

    //
    // Unzip scaffolds file
    //
    GUNZIP_SCAFFOLDS (
        SPADES.out.scaffolds
    )
    ch_versions = ch_versions.mix(GUNZIP_SCAFFOLDS.out.versions.first())

    //
    // Unzip gfa file
    //
    GUNZIP_GFA (
        SPADES.out.gfa
    )

    //
    // Filter for empty scaffold files
    //
    GUNZIP_SCAFFOLDS
        .out
        .gunzip
        .filter { meta, scaffold -> scaffold.size() > 0 }
        .set { ch_scaffolds }

    GUNZIP_GFA
        .out
        .gunzip
        .filter { meta, gfa -> gfa.size() > 0 }
        .set { ch_gfa }

    //
    // Generate assembly visualisation with Bandage
    //
    ch_bandage_png = Channel.empty()
    ch_bandage_svg = Channel.empty()
    if (!params.skip_bandage) {
        BANDAGE_IMAGE (
            ch_gfa
        )
        ch_bandage_png = BANDAGE_IMAGE.out.png
        ch_bandage_svg = BANDAGE_IMAGE.out.svg
        ch_versions    = ch_versions.mix(BANDAGE_IMAGE.out.versions.first())
    }

    //
    // Downstream assembly steps
    //
    ASSEMBLY_QC (
        ch_scaffolds,
        fasta,
        gff,
        blast_db,
        blast_header
    )
    ch_versions = ch_versions.mix(ASSEMBLY_QC.out.versions)

    emit:
    scaffolds          = SPADES.out.scaffolds               // channel: [ val(meta), [ scaffolds ] ]
    contigs            = SPADES.out.contigs                 // channel: [ val(meta), [ contigs ] ]
    transcripts        = SPADES.out.transcripts             // channel: [ val(meta), [ transcripts ] ]
    gene_clusters      = SPADES.out.gene_clusters           // channel: [ val(meta), [ gene_clusters ] ]
    gfa                = SPADES.out.gfa                     // channel: [ val(meta), [ gfa ] ]
    log_out            = SPADES.out.log                     // channel: [ val(meta), [ log ] ]

    bandage_png        = ch_bandage_png                     // channel: [ val(meta), [ png ] ]
    bandage_svg        = ch_bandage_svg                     // channel: [ val(meta), [ svg ] ]

    blast_txt          = ASSEMBLY_QC.out.blast_txt          // channel: [ val(meta), [ txt ] ]
    blast_filter_txt   = ASSEMBLY_QC.out.blast_filter_txt   // channel: [ val(meta), [ txt ] ]

    quast_results      = ASSEMBLY_QC.out.quast_results      // channel: [ val(meta), [ results ] ]
    quast_tsv          = ASSEMBLY_QC.out.quast_tsv          // channel: [ val(meta), [ tsv ] ]

    abacas_results     = ASSEMBLY_QC.out.abacas_results     // channel: [ val(meta), [ results ] ]

    plasmidid_html     = ASSEMBLY_QC.out.plasmidid_html     // channel: [ val(meta), [ html ] ]
    plasmidid_tab      = ASSEMBLY_QC.out.plasmidid_tab      // channel: [ val(meta), [ tab ] ]
    plasmidid_images   = ASSEMBLY_QC.out.plasmidid_images   // channel: [ val(meta), [ images/ ] ]
    plasmidid_logs     = ASSEMBLY_QC.out.plasmidid_logs     // channel: [ val(meta), [ logs/ ] ]
    plasmidid_data     = ASSEMBLY_QC.out.plasmidid_data     // channel: [ val(meta), [ data/ ] ]
    plasmidid_database = ASSEMBLY_QC.out.plasmidid_database // channel: [ val(meta), [ database/ ] ]
    plasmidid_fasta    = ASSEMBLY_QC.out.plasmidid_fasta    // channel: [ val(meta), [ fasta_files/ ] ]
    plasmidid_kmer     = ASSEMBLY_QC.out.plasmidid_kmer     // channel: [ val(meta), [ kmer/ ] ]

    versions           = ch_versions                        // channel: [ versions.yml ]
}
