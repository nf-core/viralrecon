/*
 * Assembly and downstream processing for minia scaffolds
 */

params.minia_options     = [:]
params.blastn_options    = [:]
params.abacas_options    = [:]
params.plasmidid_options = [:]
params.quast_options     = [:]

include { MINIA       } from '../../modules/local/minia' addParams( options: params.minia_options ) 
include { ASSEMBLY_QC } from './assembly_qc'             addParams( blastn_options: params.blastn_options, abacas_options: params.abacas_options, plasmidid_options: params.plasmidid_options, quast_options: params.quast_options )

workflow ASSEMBLY_MINIA {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    fasta         // channel: /path/to/genome.fasta
    gff           // channel: /path/to/genome.gff
    blast_db      // channel: /path/to/blast_db/
    
    main:
    /*
     * Assemble reads with minia
     */
    MINIA ( reads )

    /*
     * Filter for empty contig files
     */
    MINIA
        .out
        .contigs
        .filter { meta, contig -> contig.size() > 0 }
        .set { ch_contigs }

    /*
     * Downstream assembly steps
     */
    ASSEMBLY_QC ( 
        ch_contigs,
        fasta,
        gff,
        blast_db
    )

    emit:
    contigs           = MINIA.out.contigs                 // channel: [ val(meta), [ contigs ] ]
    unitigs           = MINIA.out.unitigs                 // channel: [ val(meta), [ unitigs ] ]
    h5                = MINIA.out.h5                      // channel: [ val(meta), [ h5 ] ]
    minia_version     = MINIA.out.version                 //    path: *.version.txt

    blast_txt         = ASSEMBLY_QC.out.blast_txt         // channel: [ val(meta), [ txt ] ]
    blast_version     = ASSEMBLY_QC.out.blast_version     //    path: *.version.txt

    quast_results     = ASSEMBLY_QC.out.quast_results     // channel: [ val(meta), [ results ] ]
    quast_tsv         = ASSEMBLY_QC.out.quast_tsv         // channel: [ val(meta), [ tsv ] ]
    quast_version     = ASSEMBLY_QC.out.quast_version     //    path: *.version.txt
    
    abacas_results    = ASSEMBLY_QC.out.abacas_results    // channel: [ val(meta), [ results ] ]
    abacas_version    = ASSEMBLY_QC.out.abacas_version    //    path: *.version.txt

    plasmidid_results = ASSEMBLY_QC.out.plasmidid_results // channel: [ val(meta), [ results ] ]
    plasmidid_version = ASSEMBLY_QC.out.plasmidid_version //    path: *.version.txt

}