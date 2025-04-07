//
// Run snpEff, bgzip, tabix, stats and SnpSift commands
//

include { SNPEFF_SNPEFF } from '../../modules/nf-core/snpeff/snpeff/main'
include { SNPSIFT_EXTRACTFIELDS } from '../../modules/local/snpsift_extractfields'

include { VCF_BGZIP_TABIX_STATS } from './vcf_bgzip_tabix_stats'

workflow SNPEFF_SNPSIFT {
    take:
    vcf    // channel: [ val(meta), [ vcf ] ]
    db     // path   : snpEff database
    config // path   : snpEff config
    fasta_path  // path   : genome.fasta

    main:

    ch_versions = Channel.empty()

    // Obtain genome ID from FASTA name
    genome_id = fasta_path.map { input ->
        def file = input instanceof List ? input.flatten()[0] : input
        def filename = file.getName()
        return filename.replaceAll(/\.f(ast|na)?(\.gz)?$/, '')
    }

    genome_ids = vcf.map { genome_id.value }
    snpeff_cache_per_sample = vcf.map { [ [ id: genome_id.value ], db.collect().value[0] ] }

    SNPEFF_SNPEFF(
        vcf,
        genome_ids,
        snpeff_cache_per_sample,
        config
    )
    ch_versions = ch_versions.mix(SNPEFF_SNPEFF.out.versions)

    VCF_BGZIP_TABIX_STATS (
        SNPEFF_SNPEFF.out.vcf,
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ]
    )
    ch_versions = ch_versions.mix(VCF_BGZIP_TABIX_STATS.out.versions)

    SNPSIFT_EXTRACTFIELDS (
        VCF_BGZIP_TABIX_STATS.out.vcf
    )
    ch_versions = ch_versions.mix(SNPSIFT_EXTRACTFIELDS.out.versions.first())

    emit:
    csv         = SNPEFF_SNPEFF.out.report           // channel: [ val(meta), [ csv ] ]
    txt         = SNPEFF_SNPEFF.out.genes_txt           // channel: [ val(meta), [ txt ] ]
    html        = SNPEFF_SNPEFF.out.summary_html         // channel: [ val(meta), [ html ] ]

    vcf         = VCF_BGZIP_TABIX_STATS.out.vcf   // channel: [ val(meta), [ vcf.gz ] ]
    tbi         = VCF_BGZIP_TABIX_STATS.out.tbi   // channel: [ val(meta), [ tbi ] ]
    csi         = VCF_BGZIP_TABIX_STATS.out.csi   // channel: [ val(meta), [ csi ] ]
    stats       = VCF_BGZIP_TABIX_STATS.out.stats // channel: [ val(meta), [ txt ] ]

    snpsift_txt = SNPSIFT_EXTRACTFIELDS.out.txt   // channel: [ val(meta), [ txt ] ]

    versions    = ch_versions                     // channel: [ versions.yml ]
}
