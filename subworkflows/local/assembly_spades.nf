/*
 * Assembly and downstream processing for SPAdes scaffolds
 */

params.spades_options       = [:]
params.bandage_options      = [:]
params.blastn_options       = [:]
params.abacas_options       = [:]
params.plasmidid_options    = [:]
params.quast_options        = [:]
params.snpeff_options       = [:]
params.snpeff_bgzip_options = [:]
params.snpeff_tabix_options = [:]
params.snpeff_stats_options = [:]
params.snpsift_options      = [:]

include { SPADES      } from '../../modules/local/spades'  addParams( options: params.spades_options  ) 
include { BANDAGE     } from '../../modules/local/bandage' addParams( options: params.bandage_options ) 
include { ASSEMBLY_QC } from './assembly_qc'               addParams( blastn_options: params.blastn_options, abacas_options: params.abacas_options, plasmidid_options: params.plasmidid_options, quast_options: params.quast_options, snpeff_options: params.snpeff_options, snpeff_bgzip_options: params.snpeff_bgzip_options, snpeff_tabix_options: params.snpeff_tabix_options, snpeff_stats_options: params.snpeff_stats_options, snpsift_options: params.snpsift_options )

workflow ASSEMBLY_SPADES {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    hmm           // channel: /path/to/spades.hmm
    coronaspades  // boolean: run coronaspades.py
    fasta         // channel: /path/to/genome.fasta
    gff           // channel: /path/to/genome.gff
    blast_db      // channel: /path/to/blast_db/
    snpeff_db     // channel: /path/to/snpeff_db/
    snpeff_config // channel: /path/to/snpeff.config
    
    main:
    /*
     * Filter for paired-end samples if running metaSPAdes / metaviralSPAdes / metaplasmidSPAdes
     */
    ch_reads = reads
    if (params.spades_options.args.contains('--meta')) {
        reads
            .filter { meta, fastq -> !meta.single_end }
            .set { ch_reads }
    }
    
    /*
     * Assemble reads with SPAdes
     */
    SPADES ( ch_reads, hmm, coronaspades )

    // FILTER CHANNELS HERE WHERE FASTA HAS NO CONTIGS
    //input: tuple val(sample), val(single_end), path(scaffold) from ch_unicycler_plasmidid.filter { it.size() > 0 }

    /*
     * Generate assembly visualisation with Bandage
     */
    ch_bandage_png     = Channel.empty()
    ch_bandage_svg     = Channel.empty()
    ch_bandage_version = Channel.empty()
    if (!params.skip_bandage) {
        BANDAGE ( SPADES.out.graph )
        ch_bandage_version = BANDAGE.out.version
        ch_bandage_png     = BANDAGE.out.png
        ch_bandage_svg     = BANDAGE.out.svg
    }

    /*
     * Downstream assembly steps
     */
    ASSEMBLY_QC ( 
        SPADES.out.scaffolds,
        fasta,
        gff,
        blast_db,
        snpeff_db,
        snpeff_config
    )

    emit:
    scaffolds         = SPADES.out.scaffolds              // channel: [ val(meta), [ scaffolds ] ]
    graph             = SPADES.out.graph                  // channel: [ val(meta), [ graph ] ]
    log_out           = SPADES.out.log                    // channel: [ val(meta), [ log ] ]
    spades_version    = SPADES.out.version                //    path: *.version.txt

    bandage_png       = ch_bandage_png                    // channel: [ val(meta), [ png ] ]
    bandage_svg       = ch_bandage_svg                    // channel: [ val(meta), [ svg ] ]
    bandage_version   = ch_bandage_version                //    path: *.version.txt

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