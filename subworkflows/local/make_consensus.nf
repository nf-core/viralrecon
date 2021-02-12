/*
 * Run various tools to generate a masked genome consensus sequence
 */

params.genomecov_options  = [:]
params.merge_options      = [:]
params.mask_options       = [:]
params.maskfasta_options  = [:]
params.bcftools_options   = [:]
params.plot_bases_options = [:]

include { BEDTOOLS_GENOMECOV } from '../../modules/nf-core/software/bedtools/genomecov/main' addParams( options: params.genomecov_options  )
include { BEDTOOLS_MERGE     } from '../../modules/nf-core/software/bedtools/merge/main'     addParams( options: params.merge_options      )
include { BEDTOOLS_MASKFASTA } from '../../modules/nf-core/software/bedtools/maskfasta/main' addParams( options: params.maskfasta_options  )
include { BCFTOOLS_CONSENSUS } from '../../modules/nf-core/software/bcftools/consensus/main' addParams( options: params.bcftools_options   )
include { MAKE_BED_MASK      } from '../../modules/local/make_bed_mask'                      addParams( options: params.mask_options       )
include { PLOT_BASE_DENSITY  } from '../../modules/local/plot_base_density'                  addParams( options: params.plot_bases_options )

workflow MAKE_CONSENSUS {
    take:
    bam_vcf // channel: [ val(meta), [ bam ], [ vcf ], [ tbi ] ]
    fasta
    
    main:
    BEDTOOLS_GENOMECOV ( bam_vcf.map { meta, bam, vcf, tbi -> [ meta, bam ] } )
    
    BEDTOOLS_MERGE ( BEDTOOLS_GENOMECOV.out.bed )

    MAKE_BED_MASK ( bam_vcf.map { meta, bam, vcf, tbi -> [ meta, vcf ] }.join( BEDTOOLS_MERGE.out.bed, by: [0] ) )
    
    BEDTOOLS_MASKFASTA ( MAKE_BED_MASK.out.bed, fasta )

    BCFTOOLS_CONSENSUS ( bam_vcf.map { meta, bam, vcf, tbi -> [ meta, vcf, tbi ] }.join( BEDTOOLS_MASKFASTA.out.fasta, by: [0] ) )

    PLOT_BASE_DENSITY ( BCFTOOLS_CONSENSUS.out.fasta )

    emit:
    fasta            = BCFTOOLS_CONSENSUS.out.fasta   // channel: [ val(meta), [ fasta ] ]
    tsv              = PLOT_BASE_DENSITY.out.tsv      // channel: [ val(meta), [ tsv ] ]
    pdf              = PLOT_BASE_DENSITY.out.pdf      // channel: [ val(meta), [ pdf ] ]
    bedtools_version = BEDTOOLS_MERGE.out.version     //    path: *.version.txt 
    bcftools_version = BCFTOOLS_CONSENSUS.out.version //    path: *.version.txt

}

