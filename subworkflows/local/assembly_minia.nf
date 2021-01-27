/*
 * Assembly and downstream processing for minia scaffolds
 */

params.minia_options        = [:]
params.blastn_options       = [:]
params.abacas_options       = [:]
params.plasmidid_options    = [:]
params.quast_options        = [:]
params.snpeff_options       = [:]
params.snpeff_bgzip_options = [:]
params.snpeff_tabix_options = [:]
params.snpeff_stats_options = [:]
params.snpsift_options      = [:]

include { MINIA   } from '../../modules/local/minia' addParams( options: params.minia_options ) 
include { ASSEMBLY_QC } from './assembly_qc'         addParams( blastn_options: params.blastn_options, abacas_options: params.abacas_options, plasmidid_options: params.plasmidid_options, quast_options: params.quast_options, snpeff_options: params.snpeff_options, snpeff_bgzip_options: params.snpeff_bgzip_options, snpeff_tabix_options: params.snpeff_tabix_options, snpeff_stats_options: params.snpeff_stats_options, snpsift_options: params.snpsift_options )

workflow ASSEMBLY_MINIA {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    fasta         // channel: /path/to/genome.fasta
    gff           // channel: /path/to/genome.gff
    blast_db      // channel: /path/to/blast_db/
    snpeff_db     // channel: /path/to/snpeff_db/
    snpeff_config // channel: /path/to/snpeff.config
    
    main:
    /*
     * Assemble reads with minia
     */
    MINIA ( reads )

    // FILTER CHANNELS HERE WHERE FASTA HAS NO CONTIGS
    //input: tuple val(sample), val(single_end), path(scaffold) from ch_unicycler_plasmidid.filter { it.size() > 0 }

    /*
     * Downstream assembly steps
     */
    ASSEMBLY_QC ( 
        MINIA.out.contigs,
        fasta,
        gff,
        blast_db,
        snpeff_db,
        snpeff_config
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