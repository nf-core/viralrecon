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
    artic_minion_caller  : ['nanopolish', 'medaka'],
    artic_minion_aligner : ['minimap2', 'bwa']
]

def checkPathParamList = [
    params.input, params.fastq_dir, params.fast5_dir,
    params.sequencing_summary, params.gff,
    params.freyja_barcodes, params.freyja_lineages, params.additional_annotation
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.fast5_dir)               { ch_fast5_dir          = file(params.fast5_dir)               } else { ch_fast5_dir          = [] }
if (params.sequencing_summary)      { ch_sequencing_summary = file(params.sequencing_summary)      } else { ch_sequencing_summary = [] }
if (params.additional_annotation)   { ch_additional_gtf     = file(params.additional_annotation)   } else { additional_annotation = [] }

// Need to stage medaka model properly depending on whether it is a string or a file
ch_medaka_model = Channel.empty()
if (params.artic_minion_caller == 'medaka') {
    if (file(params.artic_minion_medaka_model).exists()) {
        ch_medaka_model = Channel.fromPath(params.artic_minion_medaka_model)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config_nanopore.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? file(params.multiqc_config) : []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { ASCIIGENOME } from '../modules/local/asciigenome'
include { MULTIQC     } from '../modules/local/multiqc_nanopore'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_GENOME   } from '../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_AMPLICON } from '../modules/local/plot_mosdepth_regions'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_GENOME                } from '../subworkflows/local/prepare_genome_nanopore'
include { SNPEFF_SNPSIFT                } from '../subworkflows/local/snpeff_snpsift'
include { ADDITIONAL_ANNOTATION         } from '../subworkflows/local/additional_annotation'
include { VARIANTS_LONG_TABLE           } from '../subworkflows/local/variants_long_table'
include { FILTER_BAM_SAMTOOLS           } from '../subworkflows/local/filter_bam_samtools'
include { BAM_VARIANT_DEMIX_BOOT_FREYJA } from '../subworkflows/nf-core/bam_variant_demix_boot_freyja/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { PYCOQC                        } from '../modules/nf-core/pycoqc/main'
include { NANOPLOT                      } from '../modules/nf-core/nanoplot/main'
include { ARTIC_GUPPYPLEX               } from '../modules/nf-core/artic/guppyplex/main'
include { ARTIC_MINION                  } from '../modules/nf-core/artic/minion/main'
include { VCFLIB_VCFUNIQ                } from '../modules/nf-core/vcflib/vcfuniq/main'
include { TABIX_TABIX                   } from '../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_STATS                } from '../modules/nf-core/bcftools/stats/main'
include { QUAST                         } from '../modules/nf-core/quast/main'
include { PANGOLIN                      } from '../modules/nf-core/pangolin/main'
include { NEXTCLADE_RUN                 } from '../modules/nf-core/nextclade/run/main'
include { MOSDEPTH as MOSDEPTH_GENOME   } from '../modules/nf-core/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_AMPLICON } from '../modules/nf-core/mosdepth/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def pass_barcode_reads = [:]
def fail_barcode_reads = [:]

workflow NANOPORE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_genome_fasta
    ch_genome_gff
    ch_primer_bed
    ch_artic_scheme
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
    // MODULE: PycoQC on sequencing summary file
    //
    ch_pycoqc_multiqc = Channel.empty()
    if (params.sequencing_summary && !params.skip_pycoqc) {
        PYCOQC (
            Channel.of(ch_sequencing_summary).map { [ [:], it ] }
        )
        ch_pycoqc_multiqc = PYCOQC.out.json
        ch_versions       = ch_versions.mix(PYCOQC.out.versions)
    }

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

    barcode_dirs       = file("${params.fastq_dir}/barcode*", type: 'dir' , maxdepth: 1)
    single_barcode_dir = file("${params.fastq_dir}/*.fastq" , type: 'file', maxdepth: 1)
    ch_custom_no_sample_name_multiqc = Channel.empty()
    ch_custom_no_barcodes_multiqc    = Channel.empty()
    if (barcode_dirs) {
        Channel
            .fromPath( barcode_dirs )
            .filter( ~/.*barcode[0-9]{1,4}$/ )
            .map { dir ->
                def count = 0
                for (x in dir.listFiles()) {
                    if (x.isFile() && x.toString().contains('.fastq')) {
                        count += x.countFastq()
                    }
                }
                return [ dir.baseName , dir, count ]
            }
            .set { ch_fastq_dirs }

        //
        // SUBWORKFLOW: Read in samplesheet containing sample to barcode mappings
        //
        if (params.input) {
            ch_samplesheet
            .join(ch_fastq_dirs, remainder: true)
            .set { ch_fastq_dirs }

            //
            // MODULE: Create custom content file for MultiQC to report barcodes were allocated reads >= params.min_barcode_reads but no sample name in samplesheet
            //
            ch_fastq_dirs
                .filter { it[1] == null }
                .filter { it[-1] >= params.min_barcode_reads }
                .map { it -> [ "${it[0]}\t${it[-1]}" ] }
                .collect()
                .map {
                    tsv_data ->
                        def header = ['Barcode', 'Read count']
                        WorkflowCommons.multiqcTsvFromList(tsv_data, header)
                }
                .set { ch_custom_no_sample_name_multiqc }

            //
            // MODULE: Create custom content file for MultiQC to report samples that were in samplesheet but have no barcodes
            //
            ch_fastq_dirs
                .filter { it[-1] == null }
                .map { it -> [ "${it[1]}\t${it[0]}" ] }
                .collect()
                .map {
                    tsv_data ->
                        def header = ['Sample', 'Missing barcode']
                        WorkflowCommons.multiqcTsvFromList(tsv_data, header)
                }
                .set { ch_custom_no_barcodes_multiqc }

            ch_fastq_dirs
                .filter { (it[1] != null)  }
                .filter { (it[-1] != null) }
                .set { ch_fastq_dirs }

        } else {
            ch_fastq_dirs
                .map { barcode, dir, count -> [ barcode, barcode, dir, count ] }
                .set { ch_fastq_dirs }
        }
    } else if (single_barcode_dir) {
        Channel
            .fromPath("${params.fastq_dir}", type: 'dir', maxDepth: 1)
            .map { it -> [ 'SAMPLE_1', 'single_barcode', it, 10000000 ] }
            .set{ ch_fastq_dirs }
    } else {
        log.error "Please specify a valid folder containing ONT basecalled, barcoded fastq files generated by guppy_barcoder or guppy_basecaller e.g. '--fastq_dir ./20191023_1522_MC-110615_0_FAO93606_12bf9b4f/fastq_pass/"
        System.exit(1)
    }

    //
    // MODULE: Create custom content file for MultiQC to report samples with reads < params.min_barcode_reads
    //
    ch_fastq_dirs
        .branch { barcode, sample, dir, count  ->
            pass: count > params.min_barcode_reads
                pass_barcode_reads[sample] = count
                return [ "$sample\t$count" ]
            fail: count < params.min_barcode_reads
                fail_barcode_reads[sample] = count
                return [ "$sample\t$count" ]
        }
        .set { ch_pass_fail_barcode_count }

    ch_pass_fail_barcode_count
        .fail
        .collect()
        .map {
            tsv_data ->
                def header = ['Sample', 'Barcode count']
                WorkflowCommons.multiqcTsvFromList(tsv_data, header)
        }
        .set { ch_custom_fail_barcodes_count_multiqc }

    // Re-arrange channels to have meta map of information for sample
    ch_fastq_dirs
        .filter { it[-1] > params.min_barcode_reads }
        .map { barcode, sample, dir, count -> [ [ id: sample, barcode:barcode ], dir ] }
        .set { ch_fastq_dirs }

    //
    // MODULE: Run Artic Guppyplex
    //
    ARTIC_GUPPYPLEX (
        ch_fastq_dirs
    )
    ch_versions = ch_versions.mix(ARTIC_GUPPYPLEX.out.versions.first())

    //
    // MODULE: Create custom content file for MultiQC to report samples with reads < params.min_guppyplex_reads
    //
    ARTIC_GUPPYPLEX
        .out
        .fastq
        .branch { meta, fastq  ->
            def count = fastq.countFastq()
            pass: count > params.min_guppyplex_reads
                return [ "$meta.id\t$count" ]
            fail: count < params.min_guppyplex_reads
                return [ "$meta.id\t$count" ]
        }
        .set { ch_pass_fail_guppyplex_count }

    ch_pass_fail_guppyplex_count
        .fail
        .collect()
        .map {
            tsv_data ->
                def header = ['Sample', 'Read count']
                WorkflowCommons.multiqcTsvFromList(tsv_data, header)
        }
        .set { ch_custom_fail_guppyplex_count_multiqc }

    //
    // MODULE: Nanoplot QC for FastQ files
    //
    if (!params.skip_nanoplot) {
        NANOPLOT (
            ARTIC_GUPPYPLEX.out.fastq
        )
        ch_versions = ch_versions.mix(NANOPLOT.out.versions.first())
    }

    //
    // MODULE: Run Artic minion
    //
    ARTIC_MINION (
        ARTIC_GUPPYPLEX.out.fastq.filter { it[-1].countFastq() > params.min_guppyplex_reads },
        ch_fast5_dir,
        ch_sequencing_summary,
        PREPARE_GENOME.out.fasta.collect(),
        PREPARE_GENOME.out.primer_bed.collect(),
        ch_medaka_model.collect().ifEmpty([]),
        params.artic_minion_medaka_model ?: '',
        ch_artic_scheme,
        params.primer_set_version
    )
    ch_versions = ch_versions.mix(ARTIC_MINION.out.versions.first())

    //
    // MODULE: Remove duplicate variants
    //
    VCFLIB_VCFUNIQ (
        ARTIC_MINION.out.vcf.join(ARTIC_MINION.out.tbi, by: [0]),
    )
    ch_versions = ch_versions.mix(VCFLIB_VCFUNIQ.out.versions.first())

    //
    // MODULE: Index VCF file
    //
    TABIX_TABIX (
        VCFLIB_VCFUNIQ.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    //
    // MODULE: VCF stats with bcftools stats
    //
    BCFTOOLS_STATS (
        VCFLIB_VCFUNIQ.out.vcf.join(TABIX_TABIX.out.tbi, by: [0]),
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ]
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions.first())

    //
    // SUBWORKFLOW: Filter unmapped reads from BAM
    //
    FILTER_BAM_SAMTOOLS (
        ARTIC_MINION.out.bam.join(ARTIC_MINION.out.bai, by: [0]),
        [ [:], [] ]
    )
    ch_versions = ch_versions.mix(FILTER_BAM_SAMTOOLS.out.versions)

    //
    // MODULE: Genome-wide and amplicon-specific coverage QC plots
    //
    ch_mosdepth_multiqc         = Channel.empty()
    ch_amplicon_heatmap_multiqc = Channel.empty()
    if (!params.skip_mosdepth) {

        MOSDEPTH_GENOME (
            ARTIC_MINION.out.bam_primertrimmed
                .join(ARTIC_MINION.out.bai_primertrimmed, by: [0])
                .map { meta, bam, bai -> [ meta, bam, bai, [] ] },
            [ [:], [] ]
        )
        ch_mosdepth_multiqc  = MOSDEPTH_GENOME.out.global_txt
        ch_versions          = ch_versions.mix(MOSDEPTH_GENOME.out.versions.first())

        PLOT_MOSDEPTH_REGIONS_GENOME (
            MOSDEPTH_GENOME.out.regions_bed.collect { it[1] }
        )
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_GENOME.out.versions)

        MOSDEPTH_AMPLICON (
            ARTIC_MINION.out.bam_primertrimmed.join(ARTIC_MINION.out.bai_primertrimmed, by: [0]).join(PREPARE_GENOME.out.primer_collapsed_bed),
            [ [:], [] ]
        )
        ch_versions = ch_versions.mix(MOSDEPTH_AMPLICON.out.versions.first())

        PLOT_MOSDEPTH_REGIONS_AMPLICON (
            MOSDEPTH_AMPLICON.out.regions_bed.collect { it[1] }
        )
        ch_amplicon_heatmap_multiqc = PLOT_MOSDEPTH_REGIONS_AMPLICON.out.heatmap_tsv
        ch_versions                 = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_AMPLICON.out.versions)
    }

    //
    // MODULE: Lineage analysis with Pangolin
    //
    ch_pangolin_multiqc = Channel.empty()
    if (!params.skip_pangolin) {
        PANGOLIN (
            ARTIC_MINION.out.fasta
        )
        ch_pangolin_multiqc = PANGOLIN.out.report
        ch_versions         = ch_versions.mix(PANGOLIN.out.versions.first())
    }

    //
    // MODULE: Clade assignment, mutation calling, and sequence quality checks with Nextclade
    //
    ch_nextclade_multiqc = Channel.empty()
    if (!params.skip_nextclade) {
        NEXTCLADE_RUN (
            ARTIC_MINION.out.fasta,
            PREPARE_GENOME.out.nextclade_db.collect()
        )
        ch_versions = ch_versions.mix(NEXTCLADE_RUN.out.versions.first())

        //
        // MODULE: Get Nextclade clade information for MultiQC report
        //
        NEXTCLADE_RUN
            .out
            .csv
            .map {
                meta, csv ->
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
    // SUBWORKFLOW: Determine variants with Freyja
    //
    ch_freyja_multiqc = Channel.empty()
    if (!params.skip_freyja) {
        BAM_VARIANT_DEMIX_BOOT_FREYJA(
            ARTIC_MINION.out.bam_primertrimmed,
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
    // MODULE: Consensus QC across all samples with QUAST
    //
    ch_quast_multiqc = Channel.empty()
    if (!params.skip_variants_quast) {
        ARTIC_MINION.out.fasta
            .collect{ it[1] }
            .map { consensus_collect -> tuple([id: "quast"], consensus_collect) }
            .set { ch_to_quast }
        QUAST (
            ch_to_quast,
            PREPARE_GENOME.out.fasta.collect().map { [ [:], it ] },
            ch_genome_gff ? PREPARE_GENOME.out.gff.map { [ [:], it ] } : [ [:], [] ],
        )
        ch_quast_multiqc = QUAST.out.tsv
        ch_versions      = ch_versions.mix(QUAST.out.versions)
    }

    //
    // SUBWORKFLOW: Annotate variants with snpEff
    //
    ch_snpeff_multiqc = Channel.empty()
    ch_snpsift_txt    = Channel.empty()
    if (ch_genome_gff && !params.skip_snpeff) {
        SNPEFF_SNPSIFT (
            VCFLIB_VCFUNIQ.out.vcf,
            PREPARE_GENOME.out.snpeff_db.collect(),
            PREPARE_GENOME.out.snpeff_config.collect(),
            PREPARE_GENOME.out.fasta.collect()
        )
        ch_snpeff_multiqc = SNPEFF_SNPSIFT.out.csv
        ch_snpsift_txt    = SNPEFF_SNPSIFT.out.snpsift_txt
        ch_versions       = ch_versions.mix(SNPEFF_SNPSIFT.out.versions)
    }

    //
    // MODULE: Variant screenshots with ASCIIGenome
    //
    if (!params.skip_asciigenome) {
        ARTIC_MINION
            .out
            .bam_primertrimmed
            .join(VCFLIB_VCFUNIQ.out.vcf, by: [0])
            .join(BCFTOOLS_STATS.out.stats, by: [0])
            .map { meta, bam, vcf, stats ->
                if (WorkflowCommons.getNumVariantsFromBCFToolsStats(stats) > 0) {
                    return [ meta, bam, vcf ]
                }
            }
            .set { ch_asciigenome }

        ASCIIGENOME (
            ch_asciigenome,
            PREPARE_GENOME.out.fasta.collect(),
            PREPARE_GENOME.out.chrom_sizes.collect(),
            ch_genome_gff ? PREPARE_GENOME.out.gff : [],
            PREPARE_GENOME.out.primer_bed.collect(),
            params.asciigenome_window_size,
            params.asciigenome_read_depth
        )
        ch_versions = ch_versions.mix(ASCIIGENOME.out.versions.first())
    }

    //
    // SUBWORKFLOW: Create variants long table report
    //
    if (!params.skip_variants_long_table && ch_genome_gff && !params.skip_snpeff) {
        VARIANTS_LONG_TABLE (
            VCFLIB_VCFUNIQ.out.vcf,
            TABIX_TABIX.out.tbi,
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
            VCFLIB_VCFUNIQ.out.vcf,
            TABIX_TABIX.out.tbi,
            PREPARE_GENOME.out.fasta,
            ch_additional_gtf,
            ch_pangolin_multiqc

        )
        ch_versions = ch_versions.mix(ADDITIONAL_ANNOTATION.out.versions)
    }

    //
    // MODULE: Pipeline reporting
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
            ch_custom_no_sample_name_multiqc.collectFile(name: 'fail_barcodes_no_sample_mqc.tsv').ifEmpty([]),
            ch_custom_no_barcodes_multiqc.collectFile(name: 'fail_no_barcode_samples_mqc.tsv').ifEmpty([]),
            ch_custom_fail_barcodes_count_multiqc.collectFile(name: 'fail_barcode_count_samples_mqc.tsv').ifEmpty([]),
            ch_custom_fail_guppyplex_count_multiqc.collectFile(name: 'fail_guppyplex_count_samples_mqc.tsv').ifEmpty([]),
            ch_amplicon_heatmap_multiqc.ifEmpty([]),
            ch_pycoqc_multiqc.collect{it[1]}.ifEmpty([]),
            ARTIC_MINION.out.json.collect{it[1]}.ifEmpty([]),
            FILTER_BAM_SAMTOOLS.out.flagstat.collect{it[1]}.ifEmpty([]),
            BCFTOOLS_STATS.out.stats.collect{it[1]}.ifEmpty([]),
            ch_mosdepth_multiqc.collect{it[1]}.ifEmpty([]),
            ch_quast_multiqc.collect{it[1]}.ifEmpty([]),
            ch_snpeff_multiqc.collect{it[1]}.ifEmpty([]),
            ch_pangolin_multiqc.collect{it[1]}.ifEmpty([]),
            ch_nextclade_multiqc.collectFile(name: 'nextclade_clade_mqc.tsv').ifEmpty([]),
            ch_freyja_multiqc.collect{it[1]}.ifEmpty([]),
        )

        multiqc_report = MULTIQC.out.report.toList()
    }

    emit:
    multiqc_report                // channel: /path/to/multiqc_report.html
    versions       = ch_versions  // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
