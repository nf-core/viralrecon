/*
 * Assembly and downstream processing for Unicycler scaffolds
 */

params.unicycler_options = [:]
params.bandage_options   = [:]
params.blastn_options    = [:]
params.abacas_options    = [:]
params.plasmidid_options = [:]
params.quast_options     = [:]

include { UNICYCLER     } from '../../modules/local/unicycler'                     addParams( options: params.unicycler_options ) 
include { BANDAGE_IMAGE } from '../../modules/nf-core/software/bandage/image/main' addParams( options: params.bandage_options   ) 
include { ASSEMBLY_QC   } from './assembly_qc'                                     addParams( blastn_options: params.blastn_options, abacas_options: params.abacas_options, plasmidid_options: params.plasmidid_options, quast_options: params.quast_options )

workflow ASSEMBLY_UNICYCLER {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    fasta         // channel: /path/to/genome.fasta
    gff           // channel: /path/to/genome.gff
    blast_db      // channel: /path/to/blast_db/
    
    main:
    /*
     * Assemble reads with Unicycler
     */
    UNICYCLER ( reads )

    /*
     * Filter for empty scaffold files
     */
    UNICYCLER
        .out
        .scaffolds
        .filter { meta, scaffold -> scaffold.size() > 0 }
        .set { ch_scaffolds }
    
    UNICYCLER
        .out
        .gfa
        .filter { meta, gfa -> gfa.size() > 0 }
        .set { ch_gfa }

    /*
     * Generate assembly visualisation with Bandage
     */
    ch_bandage_png     = Channel.empty()
    ch_bandage_svg     = Channel.empty()
    ch_bandage_version = Channel.empty()
    if (!params.skip_bandage) {
        BANDAGE_IMAGE ( ch_gfa )
        ch_bandage_version = BANDAGE_IMAGE.out.version
        ch_bandage_png     = BANDAGE_IMAGE.out.png
        ch_bandage_svg     = BANDAGE_IMAGE.out.svg
    }

    /*
     * Downstream assembly steps
     */
    ASSEMBLY_QC ( 
        ch_scaffolds,
        fasta,
        gff,
        blast_db
    )

    emit:
    scaffolds          = UNICYCLER.out.scaffolds            // channel: [ val(meta), [ scaffolds ] ]
    gfa                = UNICYCLER.out.gfa                  // channel: [ val(meta), [ gfa ] ]
    log_out            = UNICYCLER.out.log                  // channel: [ val(meta), [ log ] ]
    unicycler_version  = UNICYCLER.out.version              //    path: *.version.txt

    bandage_png        = ch_bandage_png                     // channel: [ val(meta), [ png ] ]
    bandage_svg        = ch_bandage_svg                     // channel: [ val(meta), [ svg ] ]
    bandage_version    = ch_bandage_version                 //    path: *.version.txt

    blast_txt          = ASSEMBLY_QC.out.blast_txt          // channel: [ val(meta), [ txt ] ]
    blast_version      = ASSEMBLY_QC.out.blast_version      //    path: *.version.txt

    quast_results      = ASSEMBLY_QC.out.quast_results      // channel: [ val(meta), [ results ] ]
    quast_tsv          = ASSEMBLY_QC.out.quast_tsv          // channel: [ val(meta), [ tsv ] ]
    quast_version      = ASSEMBLY_QC.out.quast_version      //    path: *.version.txt
    
    abacas_results     = ASSEMBLY_QC.out.abacas_results     // channel: [ val(meta), [ results ] ]
    abacas_version     = ASSEMBLY_QC.out.abacas_version     //    path: *.version.txt

    plasmidid_html     = ASSEMBLY_QC.out.plasmidid_html     // channel: [ val(meta), [ html ] ]
    plasmidid_tab      = ASSEMBLY_QC.out.plasmidid_tab      // channel: [ val(meta), [ tab ] ]
    plasmidid_images   = ASSEMBLY_QC.out.plasmidid_images   // channel: [ val(meta), [ images/ ] ]
    plasmidid_logs     = ASSEMBLY_QC.out.plasmidid_logs     // channel: [ val(meta), [ logs/ ] ]
    plasmidid_data     = ASSEMBLY_QC.out.plasmidid_data     // channel: [ val(meta), [ data/ ] ]
    plasmidid_database = ASSEMBLY_QC.out.plasmidid_database // channel: [ val(meta), [ database/ ] ]
    plasmidid_fasta    = ASSEMBLY_QC.out.plasmidid_fasta    // channel: [ val(meta), [ fasta_files/ ] ]
    plasmidid_kmer     = ASSEMBLY_QC.out.plasmidid_kmer     // channel: [ val(meta), [ kmer/ ] ]
    plasmidid_version  = ASSEMBLY_QC.out.plasmidid_version  //    path: *.version.txt

}