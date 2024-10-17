/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog       } from 'plugin/nf-schema'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_viralrecon_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    protocols         : ['metagenomic', 'amplicon'],
    variant_callers   : ['ivar', 'bcftools'],
    consensus_callers : ['ivar', 'bcftools'],
    assemblers        : ['spades', 'unicycler', 'minia'],
    spades_modes      : ['rnaviral', 'corona', 'metaviral', 'meta', 'metaplasmid', 'plasmid', 'isolate', 'rna', 'bio']
]

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.fasta, params.gff, params.bowtie2_index,
    params.kraken2_db, params.primer_bed, params.primer_fasta,
    params.blast_db, params.spades_hmm, params.multiqc_config,
    params.freyja_barcodes, params.freyja_lineages, params.additional_annotation
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.input)                 { ch_input          = file(params.input)                 } else { exit 1, 'Input samplesheet file not specified!' }
if (params.spades_hmm)            { ch_spades_hmm     = file(params.spades_hmm)            } else { ch_spades_hmm = []                              }
if (params.additional_annotation) { ch_additional_gtf = file(params.additional_annotation) } else { additional_annotation = []                      }

def assemblers = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []

def variant_caller = params.variant_caller
if (!variant_caller) { variant_caller = params.protocol == 'amplicon' ? 'ivar' : 'bcftools' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config_illumina.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? file(params.multiqc_config) : []

// Header files
ch_blast_outfmt6_header          = file("$projectDir/assets/headers/blast_outfmt6_header.txt", checkIfExists: true)
ch_blast_filtered_outfmt6_header = file("$projectDir/assets/headers/blast_filtered_outfmt6_header.txt", checkIfExists: true)
ch_ivar_variants_header_mqc      = file("$projectDir/assets/headers/ivar_variants_header_mqc.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { MULTIQC                                                 } from '../modules/local/multiqc_illumina'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_GENOME   } from '../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_AMPLICON } from '../modules/local/plot_mosdepth_regions'
include { PREPARE_PRIMER_FASTA                                    } from '../modules/local/prepare_primer_fasta'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_GENOME          } from '../subworkflows/local/prepare_genome_illumina'
include { VARIANTS_IVAR           } from '../subworkflows/local/variants_ivar'
include { VARIANTS_BCFTOOLS       } from '../subworkflows/local/variants_bcftools'
include { CONSENSUS_IVAR          } from '../subworkflows/local/consensus_ivar'
include { CONSENSUS_BCFTOOLS      } from '../subworkflows/local/consensus_bcftools'
include { VARIANTS_LONG_TABLE     } from '../subworkflows/local/variants_long_table'
include { ADDITIONAL_ANNOTATION   } from '../subworkflows/local/additional_annotation'
include { ASSEMBLY_SPADES         } from '../subworkflows/local/assembly_spades'
include { ASSEMBLY_UNICYCLER      } from '../subworkflows/local/assembly_unicycler'
include { ASSEMBLY_MINIA          } from '../subworkflows/local/assembly_minia'
include { BAM_TRIM_PRIMERS_IVAR   } from '../subworkflows/local/bam_trim_primers_ivar'
include { FASTQ_TRIM_FASTP_FASTQC } from '../subworkflows/local/fastq_trim_fastp_fastqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                     } from '../modules/nf-core/cat/fastq/main'
include { CUTADAPT                      } from '../modules/nf-core/cutadapt/main'
include { FASTQC                        } from '../modules/nf-core/fastqc/main'
include { KRAKEN2_KRAKEN2               } from '../modules/nf-core/kraken2/kraken2/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from '../modules/nf-core/picard/collectmultiplemetrics/main'
include { MOSDEPTH as MOSDEPTH_GENOME   } from '../modules/nf-core/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_AMPLICON } from '../modules/nf-core/mosdepth/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { FASTQ_ALIGN_BOWTIE2           } from '../subworkflows/nf-core/fastq_align_bowtie2/main'
include { BAM_MARKDUPLICATES_PICARD     } from '../subworkflows/nf-core/bam_markduplicates_picard/main'
include { BAM_VARIANT_DEMIX_BOOT_FREYJA } from '../subworkflows/nf-core/bam_variant_demix_boot_freyja/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def pass_mapped_reads = [:]
def fail_mapped_reads = [:]

workflow ILLUMINA {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_genome_fasta
    ch_genome_gff
    ch_primer_bed
    ch_bowtie2_index
    ch_nextclade_dataset
    ch_nextclade_dataset_name
    ch_nextclade_dataset_reference
    ch_nextclade_dataset_tag

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    multiqc_report   = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME (
        ch_genome_fasta,
        ch_genome_gff,
        ch_primer_bed,
        ch_bowtie2_index,
        ch_nextclade_dataset,
        ch_nextclade_dataset_name,
        ch_nextclade_dataset_reference,
        ch_nextclade_dataset_tag
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // Check genome fasta only contains a single contig
    PREPARE_GENOME
        .out
        .fasta
        .map { WorkflowIllumina.isMultiFasta(it, log) }

    if (params.protocol == 'amplicon' && !params.skip_variants) {
        // Check primer BED file only contains suffixes provided --primer_left_suffix / --primer_right_suffix
        PREPARE_GENOME
            .out
            .primer_bed
            .map { WorkflowCommons.checkPrimerSuffixes(it, params.primer_left_suffix, params.primer_right_suffix, log) }

        // Check whether the contigs in the primer BED file are present in the reference genome
        PREPARE_GENOME
            .out
            .primer_bed
            .map { [ WorkflowCommons.getColFromFile(it, col=0, uniqify=true, sep='\t') ] }
            .set { ch_bed_contigs }

        PREPARE_GENOME
            .out
            .fai
            .map { [ WorkflowCommons.getColFromFile(it, col=0, uniqify=true, sep='\t') ] }
            .concat(ch_bed_contigs)
            .collect()
            .map { fai, bed -> WorkflowCommons.checkContigsInBED(fai, bed, log) }

        // Check whether the primer BED file supplied to the pipeline is from the SWIFT/SNAP protocol
        if (!params.ivar_trim_offset) {
            PREPARE_GENOME
                .out
                .primer_bed
                .map { WorkflowIllumina.checkIfSwiftProtocol(it, 'covid19genome', log) }
        }
    }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_samplesheet
    )
    .reads
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQ_TRIM_FASTP_FASTQC (
        ch_cat_fastq,
        [],
        false,
        params.save_trimmed_fail,
        false
    )
    ch_variants_fastq = FASTQ_TRIM_FASTP_FASTQC.out.reads
    ch_versions = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions)

    //
    // Filter empty FastQ files after adapter trimming
    //
    ch_fail_reads_multiqc = Channel.empty()
    if (!params.skip_fastp) {
        ch_variants_fastq
            .join(FASTQ_TRIM_FASTP_FASTQC.out.trim_json)
            .map {
                meta, reads, json ->
                    pass = WorkflowIllumina.getFastpReadsAfterFiltering(json) > 0
                    [ meta, reads, json, pass ]
            }
            .set { ch_pass_fail_reads }

        ch_pass_fail_reads
            .map { meta, reads, json, pass -> if (pass) [ meta, reads ] }
            .set { ch_variants_fastq }

        ch_pass_fail_reads
            .map {
                meta, reads, json, pass ->
                if (!pass) {
                    fail_mapped_reads[meta.id] = 0
                    num_reads = WorkflowIllumina.getFastpReadsBeforeFiltering(json)
                    return [ "$meta.id\t$num_reads" ]
                }
            }
            .collect()
            .map {
                tsv_data ->
                    def header = ['Sample', 'Reads before trimming']
                    WorkflowCommons.multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_reads_multiqc }
    }

    //
    // MODULE: Run Kraken2 for removal of host reads
    //
    ch_assembly_fastq  = ch_variants_fastq
    ch_kraken2_multiqc = Channel.empty()
    if (!params.skip_kraken2) {
        KRAKEN2_KRAKEN2 (
            ch_variants_fastq,
            PREPARE_GENOME.out.kraken2_db,
            params.kraken2_variants_host_filter || params.kraken2_assembly_host_filter,
            params.kraken2_variants_host_filter || params.kraken2_assembly_host_filter
        )
        ch_kraken2_multiqc = KRAKEN2_KRAKEN2.out.report
        ch_versions        = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())

        if (params.kraken2_variants_host_filter) {
            ch_variants_fastq = KRAKEN2_KRAKEN2.out.unclassified_reads_fastq
        }

        if (params.kraken2_assembly_host_filter) {
            ch_assembly_fastq = KRAKEN2_KRAKEN2.out.unclassified_reads_fastq
        }
    }

    //
    // SUBWORKFLOW: Alignment with Bowtie2
    //
    ch_bam                      = Channel.empty()
    ch_bai                      = Channel.empty()
    ch_bowtie2_multiqc          = Channel.empty()
    ch_bowtie2_flagstat_multiqc = Channel.empty()
    if (!params.skip_variants) {
        FASTQ_ALIGN_BOWTIE2 (
            ch_variants_fastq,
            PREPARE_GENOME.out.bowtie2_index,
            params.save_unaligned,
            false,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] }
        )
        ch_bam                      = FASTQ_ALIGN_BOWTIE2.out.bam
        ch_bai                      = FASTQ_ALIGN_BOWTIE2.out.bai
        ch_bowtie2_multiqc          = FASTQ_ALIGN_BOWTIE2.out.log_out
        ch_bowtie2_flagstat_multiqc = FASTQ_ALIGN_BOWTIE2.out.flagstat
        ch_versions                 = ch_versions.mix(FASTQ_ALIGN_BOWTIE2.out.versions)
    }

    //
    // Filter channels to get samples that passed Bowtie2 minimum mapped reads threshold
    //
    ch_fail_mapping_multiqc = Channel.empty()
    if (!params.skip_variants) {
        ch_bowtie2_flagstat_multiqc
            .map { meta, flagstat -> [ meta ] + WorkflowIllumina.getFlagstatMappedReads(flagstat, params) }
            .set { ch_mapped_reads }

        ch_bam
            .join(ch_mapped_reads, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_bam }

        ch_bai
            .join(ch_mapped_reads, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_bai }

        ch_mapped_reads
            .branch { meta, mapped, pass ->
                pass: pass
                    pass_mapped_reads[meta.id] = mapped
                    return [ "$meta.id\t$mapped" ]
                fail: !pass
                    fail_mapped_reads[meta.id] = mapped
                    return [ "$meta.id\t$mapped" ]
            }
            .set { ch_pass_fail_mapped }

        ch_pass_fail_mapped
            .fail
            .collect()
            .map {
                tsv_data ->
                    def header = ['Sample', 'Mapped reads']
                    WorkflowCommons.multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_mapping_multiqc }
    }

    //
    // SUBWORKFLOW: Trim primer sequences from reads with iVar
    //
    ch_ivar_trim_flagstat_multiqc = Channel.empty()
    if (!params.skip_variants && !params.skip_ivar_trim && params.protocol == 'amplicon') {
        BAM_TRIM_PRIMERS_IVAR (
            ch_bam.join(ch_bai, by: [0]),
            PREPARE_GENOME.out.primer_bed,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] }
        )
        ch_bam                        = BAM_TRIM_PRIMERS_IVAR.out.bam
        ch_bai                        = BAM_TRIM_PRIMERS_IVAR.out.bai
        ch_ivar_trim_flagstat_multiqc = BAM_TRIM_PRIMERS_IVAR.out.flagstat
        ch_versions                   = ch_versions.mix(BAM_TRIM_PRIMERS_IVAR.out.versions)
    }

    //
    // SUBWORKFLOW: Mark duplicate reads
    //
    ch_markduplicates_flagstat_multiqc = Channel.empty()
    if (!params.skip_variants && !params.skip_markduplicates) {
        BAM_MARKDUPLICATES_PICARD (
            ch_bam,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] },
            PREPARE_GENOME.out.fai
        )
        ch_bam                             = BAM_MARKDUPLICATES_PICARD.out.bam
        ch_bai                             = BAM_MARKDUPLICATES_PICARD.out.bai
        ch_markduplicates_flagstat_multiqc = BAM_MARKDUPLICATES_PICARD.out.flagstat
        ch_versions                        = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)
    }

    //
    // MODULE: Picard metrics
    //
    if (!params.skip_variants && !params.skip_picard_metrics) {
        PICARD_COLLECTMULTIPLEMETRICS (
            ch_bam.join(ch_bai, by: [0]),
            PREPARE_GENOME.out.fasta.map { [ [:], it ] },
            [ [:], [] ]
        )
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
    }

    //
    // MODULE: Genome-wide and amplicon-specific coverage QC plots
    //
    ch_mosdepth_multiqc         = Channel.empty()
    ch_amplicon_heatmap_multiqc = Channel.empty()
    if (!params.skip_variants && !params.skip_mosdepth) {
        MOSDEPTH_GENOME (
            ch_bam
                .join(ch_bai, by: [0])
                .map { meta, bam, bai -> [ meta, bam, bai, [] ] },
            [ [:], [] ],
        )
        ch_mosdepth_multiqc = MOSDEPTH_GENOME.out.global_txt
        ch_versions         = ch_versions.mix(MOSDEPTH_GENOME.out.versions.first())

        PLOT_MOSDEPTH_REGIONS_GENOME (
            MOSDEPTH_GENOME.out.regions_bed.collect { it[1] }
        )
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_GENOME.out.versions)

        if (params.protocol == 'amplicon') {
            MOSDEPTH_AMPLICON (
                ch_bam
                    .join(ch_bai, by: [0])
                    .combine(PREPARE_GENOME.out.primer_collapsed_bed),
                [ [:], [] ],
            )
            ch_versions = ch_versions.mix(MOSDEPTH_AMPLICON.out.versions.first())

            PLOT_MOSDEPTH_REGIONS_AMPLICON (
                MOSDEPTH_AMPLICON.out.regions_bed.collect { it[1] }
            )
            ch_amplicon_heatmap_multiqc = PLOT_MOSDEPTH_REGIONS_AMPLICON.out.heatmap_tsv
            ch_versions                 = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_AMPLICON.out.versions)
        }
    }

    //
    // SUBWORKFLOW: Call variants with IVar
    //
    ch_vcf                    = Channel.empty()
    ch_tbi                    = Channel.empty()
    ch_ivar_counts_multiqc    = Channel.empty()
    ch_bcftools_stats_multiqc = Channel.empty()
    ch_snpsift_txt            = Channel.empty()
    ch_snpeff_multiqc         = Channel.empty()
    if (!params.skip_variants && variant_caller == 'ivar') {
        VARIANTS_IVAR (
            ch_bam,
            PREPARE_GENOME.out.fasta,
            (params.protocol == 'amplicon' || !params.skip_asciigenome || !params.skip_markduplicates) ? PREPARE_GENOME.out.fai : [],
            (params.protocol == 'amplicon' || !params.skip_asciigenome || !params.skip_markduplicates) ? PREPARE_GENOME.out.chrom_sizes : [],
            ch_genome_gff ? PREPARE_GENOME.out.gff : [],
            (params.protocol == 'amplicon' && ch_primer_bed) ? PREPARE_GENOME.out.primer_bed : [],
            PREPARE_GENOME.out.snpeff_db,
            PREPARE_GENOME.out.snpeff_config,
            ch_ivar_variants_header_mqc
        )
        ch_vcf                    = VARIANTS_IVAR.out.vcf
        ch_tbi                    = VARIANTS_IVAR.out.tbi
        ch_ivar_counts_multiqc    = VARIANTS_IVAR.out.multiqc_tsv
        ch_bcftools_stats_multiqc = VARIANTS_IVAR.out.stats
        ch_snpeff_multiqc         = VARIANTS_IVAR.out.snpeff_csv
        ch_snpsift_txt            = VARIANTS_IVAR.out.snpsift_txt
        ch_versions               = ch_versions.mix(VARIANTS_IVAR.out.versions)
    }

    //
    // SUBWORKFLOW: Call variants with BCFTools
    //
    if (!params.skip_variants && variant_caller == 'bcftools') {
        VARIANTS_BCFTOOLS (
            ch_bam,
            PREPARE_GENOME.out.fasta,
            (params.protocol == 'amplicon' || !params.skip_asciigenome || !params.skip_markduplicates) ? PREPARE_GENOME.out.chrom_sizes : [],
            ch_genome_gff ? PREPARE_GENOME.out.gff : [],
            (params.protocol == 'amplicon' && ch_primer_bed) ? PREPARE_GENOME.out.primer_bed : [],
            PREPARE_GENOME.out.snpeff_db,
            PREPARE_GENOME.out.snpeff_config
        )
        ch_vcf                    = VARIANTS_BCFTOOLS.out.vcf
        ch_tbi                    = VARIANTS_BCFTOOLS.out.tbi
        ch_bcftools_stats_multiqc = VARIANTS_BCFTOOLS.out.stats
        ch_snpeff_multiqc         = VARIANTS_BCFTOOLS.out.snpeff_csv
        ch_snpsift_txt            = VARIANTS_BCFTOOLS.out.snpsift_txt
        ch_versions               = ch_versions.mix(VARIANTS_BCFTOOLS.out.versions)
    }

    //
    // SUBWORKFLOW: Determine variants with Freyja
    //
    ch_freyja_multiqc = Channel.empty()
    if (!params.skip_variants && !params.skip_freyja) {
        BAM_VARIANT_DEMIX_BOOT_FREYJA(
            ch_bam,
            PREPARE_GENOME.out.fasta,
            params.skip_freyja_boot,
            params.freyja_repeats,
            params.freyja_db_name,
            params.freyja_barcodes,
            params.freyja_lineages,
        )
        ch_versions       = ch_versions.mix(BAM_VARIANT_DEMIX_BOOT_FREYJA.out.versions)
        ch_freyja_multiqc = BAM_VARIANT_DEMIX_BOOT_FREYJA.out.demix
    }

    //
    // SUBWORKFLOW: Call consensus with iVar and downstream QC
    //
    ch_quast_multiqc    = Channel.empty()
    ch_pangolin_multiqc = Channel.empty()
    ch_nextclade_report = Channel.empty()
    if (!params.skip_consensus && params.consensus_caller == 'ivar') {
        CONSENSUS_IVAR (
            ch_bam,
            PREPARE_GENOME.out.fasta,
            ch_genome_gff ? PREPARE_GENOME.out.gff.map { [ [:], it ] } : [ [:], [] ],
            PREPARE_GENOME.out.nextclade_db
        )

        ch_quast_multiqc    = CONSENSUS_IVAR.out.quast_tsv
        ch_pangolin_multiqc = CONSENSUS_IVAR.out.pangolin_report
        ch_nextclade_report = CONSENSUS_IVAR.out.nextclade_report
        ch_versions         = ch_versions.mix(CONSENSUS_IVAR.out.versions)
    }

    //
    // SUBWORKFLOW: Call consensus with BCFTools
    //
    if (!params.skip_consensus && params.consensus_caller == 'bcftools' && variant_caller) {
        CONSENSUS_BCFTOOLS (
            ch_bam,
            ch_vcf,
            ch_tbi,
            PREPARE_GENOME.out.fasta,
            ch_genome_gff ? PREPARE_GENOME.out.gff.map { [ [:], it ] } : [ [:], [] ],
            PREPARE_GENOME.out.nextclade_db
        )

        ch_quast_multiqc    = CONSENSUS_BCFTOOLS.out.quast_tsv
        ch_pangolin_multiqc = CONSENSUS_BCFTOOLS.out.pangolin_report
        ch_nextclade_report = CONSENSUS_BCFTOOLS.out.nextclade_report
        ch_versions         = ch_versions.mix(CONSENSUS_BCFTOOLS.out.versions)
    }

    //
    // MODULE: Get Nextclade clade information for MultiQC report
    //
    ch_nextclade_multiqc = Channel.empty()
    if (!params.skip_nextclade) {
        ch_nextclade_report
            .map { meta, csv ->
                def clade = WorkflowCommons.getNextcladeFieldMapFromCsv(csv)['clade']
                return [ "$meta.id\t$clade" ]
            }
            .collect()
            .map {
                tsv_data ->
                    def header = ['Sample', 'clade']
                    WorkflowCommons.multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_nextclade_multiqc }
    }

    //
    // SUBWORKFLOW: Create variants long table report
    //
    if (!params.skip_variants && !params.skip_variants_long_table && ch_genome_gff && !params.skip_snpeff) {
        VARIANTS_LONG_TABLE (
            ch_vcf,
            ch_tbi,
            ch_snpsift_txt,
            ch_pangolin_multiqc
        )
        ch_versions = ch_versions.mix(VARIANTS_LONG_TABLE.out.versions)
    }

    //
    // SUBWORKFLOW: Create variants long table report for additional annotation file
    //
    if (params.additional_annotation) {
        ADDITIONAL_ANNOTATION (
            ch_vcf,
            ch_tbi,
            PREPARE_GENOME.out.fasta,
            ch_additional_gtf,
            ch_pangolin_multiqc

        )
        ch_versions = ch_versions.mix(ADDITIONAL_ANNOTATION.out.versions)
    }

    //
    // MODULE: Primer trimming with Cutadapt
    //
    ch_cutadapt_multiqc = Channel.empty()
    if (params.protocol == 'amplicon' && !params.skip_assembly && !params.skip_cutadapt) {
        if (!params.skip_noninternal_primers){
            PREPARE_PRIMER_FASTA(
                PREPARE_GENOME.out.primer_fasta.collect { it[1] }
                )
            ch_assembly_fastq
                .map { info, reads ->
                    def meta = info +
                        [primers: PREPARE_PRIMER_FASTA.out.adapters.value]
                    return [meta, reads] }
                .set{ ch_assembly_fastq }
        } else {
            ch_assembly_fastq
                .map { info, reads ->
                    def meta = info +
                        [primers: PREPARE_GENOME.out.primer_fasta.value[1]]
                    return [meta, reads] }
                .set{ ch_assembly_fastq }
        }

        CUTADAPT (
            ch_assembly_fastq
        )
        ch_assembly_fastq   = CUTADAPT.out.reads
        ch_cutadapt_multiqc = CUTADAPT.out.log
        ch_versions         = ch_versions.mix(CUTADAPT.out.versions.first())

        if (!params.skip_fastqc) {
            FASTQC (
                CUTADAPT.out.reads
            )
            ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        }
    }

    //
    // SUBWORKFLOW: Run SPAdes assembly and downstream analysis
    //
    ch_spades_quast_multiqc = Channel.empty()
    if (!params.skip_assembly && 'spades' in assemblers) {
        ASSEMBLY_SPADES (
            ch_assembly_fastq.map { meta, fastq -> [ meta, fastq, [], [] ] },
            params.spades_mode,
            ch_spades_hmm,
            PREPARE_GENOME.out.fasta,
            ch_genome_gff ? PREPARE_GENOME.out.gff.map { [ [:], it ] } : [ [:], [] ],
            PREPARE_GENOME.out.blast_db,
            ch_blast_outfmt6_header,
            ch_blast_filtered_outfmt6_header
        )
        ch_spades_quast_multiqc = ASSEMBLY_SPADES.out.quast_tsv
        ch_versions             = ch_versions.mix(ASSEMBLY_SPADES.out.versions)
    }

    //
    // SUBWORKFLOW: Run Unicycler assembly and downstream analysis
    //
    ch_unicycler_quast_multiqc = Channel.empty()
    if (!params.skip_assembly && 'unicycler' in assemblers) {
        ASSEMBLY_UNICYCLER (
            ch_assembly_fastq.map { meta, fastq -> [ meta, fastq, [] ] },
            PREPARE_GENOME.out.fasta,
            ch_genome_gff ? PREPARE_GENOME.out.gff.map { [ [:], it ] } : [ [:], [] ],
            PREPARE_GENOME.out.blast_db,
            ch_blast_outfmt6_header,
            ch_blast_filtered_outfmt6_header
        )
        ch_unicycler_quast_multiqc = ASSEMBLY_UNICYCLER.out.quast_tsv
        ch_versions                = ch_versions.mix(ASSEMBLY_UNICYCLER.out.versions)
    }

    //
    // SUBWORKFLOW: Run minia assembly and downstream analysis
    //
    ch_minia_quast_multiqc = Channel.empty()
    if (!params.skip_assembly && 'minia' in assemblers) {
        ASSEMBLY_MINIA (
            ch_assembly_fastq,
            PREPARE_GENOME.out.fasta,
            ch_genome_gff ? PREPARE_GENOME.out.gff.map { [ [:], it ] } : [ [:], [] ],
            PREPARE_GENOME.out.blast_db,
            ch_blast_outfmt6_header,
            ch_blast_filtered_outfmt6_header
        )
        ch_minia_quast_multiqc = ASSEMBLY_MINIA.out.quast_tsv
        ch_versions            = ch_versions.mix(ASSEMBLY_MINIA.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        summary_params                        = paramsSummaryMap(
            workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
            file(params.multiqc_methods_description, checkIfExists: true) :
            file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = Channel.value(
            methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_logo                       = params.multiqc_logo ?
            Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
            Channel.empty()

        ch_multiqc_files                      = ch_multiqc_files.mix(
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files                      = ch_multiqc_files.mix(
            ch_methods_description.collectFile(
                name: 'methods_description_mqc.yaml',
                sort: false))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config,
            ch_multiqc_custom_config,
            ch_multiqc_logo.toList(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_fail_reads_multiqc.collectFile(name: 'fail_mapped_reads_mqc.tsv').ifEmpty([]),
            ch_fail_mapping_multiqc.collectFile(name: 'fail_mapped_samples_mqc.tsv').ifEmpty([]),
            ch_amplicon_heatmap_multiqc.ifEmpty([]),
            FASTQ_TRIM_FASTP_FASTQC.out.fastqc_raw_zip.collect{it[1]}.ifEmpty([]),
            FASTQ_TRIM_FASTP_FASTQC.out.trim_json.collect{it[1]}.ifEmpty([]),
            ch_kraken2_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_flagstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_multiqc.collect{it[1]}.ifEmpty([]),
            ch_ivar_trim_flagstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_flagstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_mosdepth_multiqc.collect{it[1]}.ifEmpty([]),
            ch_ivar_counts_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bcftools_stats_multiqc.collect{it[1]}.ifEmpty([]),
            ch_snpeff_multiqc.collect{it[1]}.ifEmpty([]),
            ch_quast_multiqc.collect{it[1]}.ifEmpty([]),
            ch_pangolin_multiqc.collect{it[1]}.ifEmpty([]),
            ch_nextclade_multiqc.collectFile(name: 'nextclade_clade_mqc.tsv').ifEmpty([]),
            ch_cutadapt_multiqc.collect{it[1]}.ifEmpty([]),
            ch_spades_quast_multiqc.collect{it[1]}.ifEmpty([]),
            ch_unicycler_quast_multiqc.collect{it[1]}.ifEmpty([]),
            ch_minia_quast_multiqc.collect{it[1]}.ifEmpty([]),
            ch_freyja_multiqc.collect{it[1]}.ifEmpty([]),
        )

        multiqc_report = MULTIQC.out.report.toList()
    }

    emit:
    multiqc_report                  // channel: /path/to/multiqc_report.html
    versions         = ch_versions  // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
