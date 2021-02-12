/*
 * Downstream analysis for assembly scaffolds using variant graphs
 */

params.snpeff_options       = [:]
params.snpeff_bgzip_options = [:]
params.snpeff_tabix_options = [:]
params.snpeff_stats_options = [:]
params.snpsift_options      = [:]

include { SNPEFF_SNPSIFT } from './snpeff_snpsift'                 addParams( snpeff_options: params.snpeff_options, snpsift_options: params.snpsift_options, bgzip_options: params.snpeff_bgzip_options, tabix_options: params.snpeff_tabix_options, stats_options: params.snpeff_stats_options )

workflow ASSEMBLY_VG {
    take:
    scaffolds     // channel: [ val(meta), [ scaffolds ] ]
    fasta         // channel: /path/to/genome.fasta
    snpeff_db     // channel: /path/to/snpeff_db/
    snpeff_config // channel: /path/to/snpeff_db/
    
    main:

    // if (!params.skip_assembly_snpeff) {
    //     SNPEFF_SNPSIFT ( vcf, snpeff_db, snpeff_config, fasta )
    // }

    // emit:
    // blast_txt         = ch_blast_txt         // channel: [ val(meta), [ txt ] ]
    // blast_version     = ch_blast_version     //    path: *.version.txt

}

