//
// Assembly and downstream processing for SPAdes scaffolds
//

include { SPADES                     } from '../../../modules/nf-core/spades/main'
include { BANDAGE_IMAGE              } from '../../../modules/nf-core/bandage/image/main'
include { GUNZIP as GUNZIP_SCAFFOLDS } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFA       } from '../../../modules/nf-core/gunzip/main'

include { ASSEMBLY_QC   } from '../assembly_qc'
include { BLAST_REPORT  } from '../../../modules/local/blast_report'

workflow ASSEMBLY_SPADES {
    take:
    reads                 // channel: [ val(meta), [ reads ] ]
    mode                  // string : spades assembly mode e.g. 'rnaviral'
    hmm                   // channel: /path/to/spades.hmm
    fasta                 // channel: /path/to/genome.fasta
    gff                   // channel: /path/to/genome.gff
    blast_db              // channel: /path/to/blast_db/
    blast_header          // channel: /path/to/blast_header.txt
    blast_filtered_header // channel: /path/to/blast_filtered_header.txt
    ch_taxidlist          // channel: /path/to/taxidlist.txt

    main:

    ch_versions = channel.empty()

    //
    // Filter for paired-end samples if running metaSPAdes / metaviralSPAdes / metaplasmidSPAdes
    //
    ch_reads = reads
    if (mode.contains('meta') || mode.contains('bio')) {
        reads
            .filter { meta, _illumina, _pacbio, _nanopore -> !meta.single_end }
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

    SPADES.out.scaffolds
        .mix(SPADES.out.contigs)
        .groupTuple(by: 0)
        .map { meta, files ->
            // Choose scaffold if it exists and is not empty, otherwise contig
            def scaffold = files.find { it.name.contains('scaffold') }
            def contig = files.find { it.name.contains('contig') }

            def assembly = scaffold ? scaffold : contig

            if (!assembly) {
                error "No assembly found for sample ${meta}"
            }

            [meta, file(assembly)]
        }
        .set { ch_assembly }
    //
    // Unzip scaffolds file
    //
    GUNZIP_SCAFFOLDS (
        ch_assembly
    )

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
        .filter { _meta, scaffold -> scaffold.size() > 0 }
        .set { ch_scaffolds }

    GUNZIP_GFA
        .out
        .gunzip
        .filter { _meta, gfa -> gfa.size() > 0 }
        .set { ch_gfa }

    //
    // Generate assembly visualisation with Bandage
    //
    ch_bandage_png = channel.empty()
    ch_bandage_svg = channel.empty()
    if (!params.skip_bandage) {
        BANDAGE_IMAGE (
            ch_gfa
        )
        ch_bandage_png = BANDAGE_IMAGE.out.png
        ch_bandage_svg = BANDAGE_IMAGE.out.svg
    }

    //
    // Downstream assembly steps
    //
    ASSEMBLY_QC (
        ch_scaffolds,
        fasta,
        gff,
        blast_db,
        blast_header,
        blast_filtered_header,
        ch_taxidlist
    )
    ch_versions = ch_versions.mix(ASSEMBLY_QC.out.versions)

    ch_blast_report = channel.empty()
    ch_reversed_fasta = channel.empty()
    ch_genotype = channel.empty()

    if (!params.skip_blast && (params.genome == 'NC_002058.3' || params.perform_ev_typing)) {
        ch_blast_report_input = ASSEMBLY_QC.out.blast_txt.join(ch_scaffolds, by: [0])
            .filter{ _meta, blast, _assembly_fasta -> blast.countLines() > 1 }
        BLAST_REPORT (
            ch_blast_report_input
        )
        ch_blast_report = BLAST_REPORT.out.blast_report
        ch_reversed_fasta = BLAST_REPORT.out.reversed_contigs
        ch_genotype = BLAST_REPORT.out.genotype
    }

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
    blast_report       = ch_blast_report                    // channel: [ val(meta), [ html ] ]
    reversed_fasta     = ch_reversed_fasta                  // channel: [ val(meta), [ fasta ] ]
    genotype           = ch_genotype                        // channel: [ val(meta), [ csv ] ]

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

    versions           = ch_versions                       // channel: versions.yml
}
