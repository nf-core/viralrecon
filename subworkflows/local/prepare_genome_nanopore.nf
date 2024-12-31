//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA      } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF        } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRIMER_BED } from '../../modules/nf-core/gunzip/main'
include { UNTAR                       } from '../../modules/nf-core/untar/main'
include { CUSTOM_GETCHROMSIZES        } from '../../modules/nf-core/custom/getchromsizes/main'
include { NEXTCLADE_DATASETGET        } from '../../modules/nf-core/nextclade/datasetget/main'
include { COLLAPSE_PRIMERS            } from '../../modules/local/collapse_primers'
include { SNPEFF_BUILD                } from '../../modules/local/snpeff_build'

workflow PREPARE_GENOME {

    take:
    fasta
    gff
    primer_bed
    bowtie2_index
    nextclade_dataset
    nextclade_dataset_name
    nextclade_dataset_tag

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        GUNZIP_FASTA (
            [ [:], fasta ]
        )
        ch_fasta    = GUNZIP_FASTA.out.gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta))
    }

    //
    // Uncompress GFF annotation file
    //
    ch_gff = Channel.empty()
    if (gff) {
        if (gff.endsWith('.gz')) {
            GUNZIP_GFF (
                [ [:], gff ]
            )
            ch_gff      = GUNZIP_GFF.out.gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = Channel.value(file(gff))
        }
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES (
        ch_fasta.map { [ [:], it ] }
    )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Uncompress primer BED file
    //
    ch_primer_bed = Channel.empty()
    if (primer_bed) {
        if (primer_bed.endsWith('.gz')) {
            GUNZIP_PRIMER_BED (
                [ [:], primer_bed ]
            )
            ch_primer_bed = GUNZIP_PRIMER_BED.out.gunzip.map { it[1] }
            ch_versions   = ch_versions.mix(GUNZIP_PRIMER_BED.out.versions)
        } else {
            ch_primer_bed = Channel.value(file(primer_bed))
        }
    }

    //
    // Generate collapsed BED file
    //
    ch_primer_collapsed_bed = Channel.empty()
    if (!params.skip_mosdepth) {
        COLLAPSE_PRIMERS (
            ch_primer_bed,
            params.primer_left_suffix,
            params.primer_right_suffix
        )
        ch_primer_collapsed_bed = COLLAPSE_PRIMERS.out.bed
        ch_versions             = ch_versions.mix(COLLAPSE_PRIMERS.out.versions)
    }

    //
    // Prepare Nextclade dataset
    //
    ch_nextclade_db = Channel.empty()
    if (!params.skip_consensus && !params.skip_nextclade) {
        if (nextclade_dataset) {
            if (nextclade_dataset.endsWith('.tar.gz')) {
                UNTAR (
                    [ [:], nextclade_dataset ]
                )
                ch_nextclade_db = UNTAR.out.untar.map { it[1] }
                ch_versions     = ch_versions.mix(UNTAR.out.versions)
            } else {
                ch_nextclade_db = Channel.value(file(nextclade_dataset))
            }
        } else if (nextclade_dataset_name) {
            NEXTCLADE_DATASETGET (
                nextclade_dataset_name,
                nextclade_dataset_tag
            )
            ch_nextclade_db = NEXTCLADE_DATASETGET.out.dataset
            ch_versions     = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)
        }
    }

    //
    // Make snpEff database
    //
    ch_snpeff_db     = Channel.empty()
    ch_snpeff_config = Channel.empty()
    if (!params.skip_snpeff) {
        SNPEFF_BUILD (
            ch_fasta,
            ch_gff
        )
        ch_snpeff_db     = SNPEFF_BUILD.out.db
        ch_snpeff_config = SNPEFF_BUILD.out.config
        ch_versions      = ch_versions.mix(SNPEFF_BUILD.out.versions)
    }

    emit:
    fasta                = ch_fasta                // path: genome.fasta
    gff                  = ch_gff                  // path: genome.gff
    fai                  = ch_fai                  // path: genome.fai
    chrom_sizes          = ch_chrom_sizes          // path: genome.sizes
    primer_bed           = ch_primer_bed           // path: primer.bed
    primer_collapsed_bed = ch_primer_collapsed_bed // path: primer.collapsed.bed
    nextclade_db         = ch_nextclade_db         // path: nextclade_db
    snpeff_db            = ch_snpeff_db            // path: snpeff_db
    snpeff_config        = ch_snpeff_config        // path: snpeff.config

    versions             = ch_versions             // channel: [ versions.yml ]
}
