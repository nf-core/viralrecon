

// /*
//  * PREPROCESSING: Build SnpEff database for viral genome
//  */
// process MAKE_SNPEFF_DB {
//     tag "${index_base}.fa"
//     label 'process_low'
//     if (params.save_reference) {
//         publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
//     }

//     when:
//     (!params.skip_variants || !params.skip_assembly) && params.gff && !params.skip_snpeff

//     input:
//     path ("SnpEffDB/genomes/${index_base}.fa") from ch_fasta
//     path ("SnpEffDB/${index_base}/genes.gff") from ch_gff

//     output:
//     tuple path("SnpEffDB"), path("*.config") into ch_snpeff_db_varscan2,
//                                                   ch_snpeff_db_ivar,
//                                                   ch_snpeff_db_bcftools,
//                                                   ch_snpeff_db_spades,
//                                                   ch_snpeff_db_metaspades,
//                                                   ch_snpeff_db_unicycler,
//                                                   ch_snpeff_db_minia

//     script:
//     """
//     echo "${index_base}.genome : ${index_base}" > snpeff.config
//     snpEff build -config snpeff.config -dataDir ./SnpEffDB -gff3 -v ${index_base}
//     """
// }
