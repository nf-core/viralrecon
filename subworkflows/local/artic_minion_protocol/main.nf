//
// Run snpEff, bgzip, tabix, stats and SnpSift commands
//

include { ARTIC_MINION   } from '../../../modules/nf-core/artic/minion/main'
include { VCFLIB_VCFUNIQ } from '../../../modules/nf-core/vcflib/vcfuniq/main'
include { BCFTOOLS_INDEX } from '../../../modules/nf-core/bcftools/index/main'

workflow ARTIC_MINION_PROTOCOL {
    take:
    reads      // channel: [ val(meta), [ fastq ] ]
    model      // channel: [ val(meta), model_dir, model ]
    fasta_bed  // channel: [ val(meta), fasta, bed ]

    main:

    ARTIC_MINION (
        reads,
        model,
        fasta_bed,
        []
    )

    //
    // MODULE: Remove duplicate variants
    //

    VCFLIB_VCFUNIQ (
        ARTIC_MINION.out.vcf.join(ARTIC_MINION.out.tbi, by: [0])
    )

    //
    // MODULE: Index VCF file
    //
    BCFTOOLS_INDEX (
        VCFLIB_VCFUNIQ.out.vcf
    )

    emit:
    bam       = ARTIC_MINION.out.bam_primertrimmed
    bai       = ARTIC_MINION.out.bai_primertrimmed

    vcf       = VCFLIB_VCFUNIQ.out.vcf
    tbi       = BCFTOOLS_INDEX.out.index

    consensus = ARTIC_MINION.out.fasta

    artic_minion_report = ARTIC_MINION.out.json

}
