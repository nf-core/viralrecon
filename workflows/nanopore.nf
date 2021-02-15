////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

def valid_params = [
    artic_minion_caller   : ['nanopolish', 'medaka'],
    artic_minion_aligner  : ['minimap2', 'bwa']
]

// Validate input parameters
Workflow.validate_nanopore_params(params, log, valid_params)

def checkPathParamList = [
    params.input, params.fastq_dir, params.fast5_dir, 
    params.sequencing_summary, params.gff
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

if (params.input)              { ch_input              = file(params.input)              }
if (params.gff)                { ch_gff                = file(params.gff)                } else { ch_gff                = ch_dummy_file }
if (params.fast5_dir)          { ch_fast5_dir          = file(params.fast5_dir)          } else { ch_fast5_dir          = ch_dummy_file }
if (params.sequencing_summary) { ch_sequencing_summary = file(params.sequencing_summary) } else { ch_sequencing_summary = ch_dummy_file }

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config_nanopore.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def artic_minion_options   = modules['artic_minion']
artic_minion_options.args += params.artic_minion_caller  == 'medaka' ? " --medaka" : ""
artic_minion_options.args += params.artic_minion_aligner == 'bwa'    ? " --bwa"    : " --minimap2"

def multiqc_options   = modules['artic_multiqc']
multiqc_options.args += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

include { ARTIC_GUPPYPLEX                     } from '../modules/local/artic_guppyplex'                     addParams( options: modules['artic_guppyplex']      )
include { ARTIC_MINION                        } from '../modules/local/artic_minion'                        addParams( options: artic_minion_options            )
include { MULTIQC_CUSTOM_FAIL_NO_BARCODES     } from '../modules/local/multiqc_custom_fail_no_barcodes'     addParams( options: [publish_files: false]          )
include { MULTIQC_CUSTOM_FAIL_BARCODE_COUNT   } from '../modules/local/multiqc_custom_fail_barcode_count'   addParams( options: [publish_files: false]          )
include { MULTIQC_CUSTOM_FAIL_GUPPYPLEX_COUNT } from '../modules/local/multiqc_custom_fail_guppyplex_count' addParams( options: [publish_files: false]          )
include { GET_SOFTWARE_VERSIONS               } from '../modules/local/get_software_versions'               addParams( options: [publish_files: ['csv':'']]     )
include { MULTIQC                             } from '../modules/local/multiqc_nanopore'                    addParams( options: multiqc_options                 )

/*
 * SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
 */
def publish_genome_options   = params.save_reference ? [publish_dir: 'genome'] : [publish_files: false]
def clone_scheme_options     = modules['artic_clone_scheme']
def collapse_primers_options = modules['artic_collapse_primers']
def snpeff_build_options     = modules['artic_snpeff_build']
if (!params.save_reference) {
    clone_scheme_options['publish_files']  = false
    collapse_primers_options['publish_files']  = false
    snpeff_build_options['publish_files']      = false
}

include { INPUT_CHECK    } from '../subworkflows/local/input_check'             addParams( options: [:] )
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome_nanopore' addParams( genome_options: publish_genome_options, clone_scheme_options: clone_scheme_options, collapse_primers_options: collapse_primers_options, snpeff_build_options: snpeff_build_options )
include { SNPEFF_SNPSIFT } from '../subworkflows/local/snpeff_snpsift'          addParams( snpeff_options: modules['artic_snpeff'], snpsift_options: modules['artic_snpsift'], bgzip_options: modules['artic_snpeff_bgzip'], tabix_options: modules['artic_snpeff_tabix'], stats_options: modules['artic_snpeff_stats'] )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

/*
 * MODULE: Installed directly from nf-core/modules
 */
include { BCFTOOLS_STATS } from '../modules/nf-core/software/bcftools/stats/main' addParams( options: modules['artic_bcftools_stats'] )
include { QUAST          } from '../modules/nf-core/software/quast/main'          addParams( options: modules['artic_quast']          )
include { PANGOLIN       } from '../modules/nf-core/software/pangolin/main'       addParams( options: modules['artic_pangolin']       )

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */
include { FILTER_BAM_SAMTOOLS } from '../subworkflows/nf-core/filter_bam_samtools' addParams( samtools_view_options: modules['artic_filter_bam'], samtools_index_options: modules['artic_filter_bam_stats'] )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report     = []
def pass_barcode_reads = [:]
def fail_barcode_reads = [:]

workflow NANOPORE {

    /*
     * SUBWORKFLOW: Uncompress and prepare reference genome files
     */
    PREPARE_GENOME (
        ch_dummy_file
    )

    // Check primer BED file only contains suffixes provided --primer_left_suffix / --primer_right_suffix
    Workflow.check_primer_suffixes(PREPARE_GENOME.out.primer_bed, params.primer_left_suffix, params.primer_right_suffix, log)
    
    barcode_dirs       = file("${params.fastq_dir}/barcode*", type: 'dir' , maxdepth: 1)
    single_barcode_dir = file("${params.fastq_dir}/*.fastq" , type: 'file', maxdepth: 1)
    if (barcode_dirs) {
        Channel
            .fromPath( barcode_dirs )
            .filter( ~/.*barcode[0-9]{1,4}$/ )
            .map { dir ->
                def count = 0
                for (x in dir.listFiles()) {
                    if (x.isFile()) {
                        count += x.countFastq()
                    }
                }
                return [ dir.baseName , dir, count ]
            }
            .set { ch_fastq_dirs }

        /*
         * SUBWORKFLOW: Read in samplesheet containing sample to barcode mappings
         */
        if (params.input) {
            INPUT_CHECK ( 
                ch_input,
                params.platform
            )
            .join(ch_fastq_dirs, remainder: true)
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
    
    /*
     * MODULE: Create custom content file for MultiQC to report samples that were in samplesheet but have no barcodes
     */
    ch_fastq_dirs
        .filter { it[-1] == null }
        .map { it -> [ "${it[1]}\t${it[0]}" ] }
        .set { ch_samples_no_barcode }

    MULTIQC_CUSTOM_FAIL_NO_BARCODES ( 
        ch_samples_no_barcode.collect()
    )
    
    /*
     * MODULE: Create custom content file for MultiQC to report samples with reads < params.min_barcode_reads
     */
    ch_fastq_dirs
        .filter { it[-1] != null }
        .branch { barcode, sample, dir, count  ->
            pass: count > params.min_barcode_reads
                pass_barcode_reads[sample] = count
                return [ "$sample\t$count" ]
            fail: count < params.min_barcode_reads
                fail_barcode_reads[sample] = count
                return [ "$sample\t$count" ]
        }
        .set { ch_pass_fail_barcode_count }

    MULTIQC_CUSTOM_FAIL_BARCODE_COUNT ( 
        ch_pass_fail_barcode_count.fail.collect()
    )

    // Re-arrange channels to have meta map of information for sample
    ch_fastq_dirs
        .filter { it[-1] > params.min_barcode_reads }
        .map { barcode, sample, dir, count -> [ [ id: sample, barcode:barcode ], dir ] }
        .set { ch_fastq_dirs }

    /*
     * MODULE: Run Artic Guppyplex
     */
    ARTIC_GUPPYPLEX (
        ch_fastq_dirs
    )
    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(ARTIC_GUPPYPLEX.out.version.first().ifEmpty(null))

    /*
     * MODULE: Create custom content file for MultiQC to report samples with reads < params.min_guppyplex_reads
     */
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

    MULTIQC_CUSTOM_FAIL_GUPPYPLEX_COUNT ( 
        ch_pass_fail_guppyplex_count.fail.collect()
    )
    
    /*
     * MODULE: Run Artic minion
     */
    ARTIC_MINION (
        ARTIC_GUPPYPLEX.out.fastq.filter { it[-1].countFastq() > params.min_guppyplex_reads },
        ch_fast5_dir,
        ch_sequencing_summary,
        PREPARE_GENOME.out.scheme_dir,
        params.artic_scheme,
        params.artic_scheme_version
    )

    /*
     * SUBWORKFLOW: Filter unmapped reads from BAM and trim reads with iVar
     */
    FILTER_BAM_SAMTOOLS ( 
        ARTIC_MINION.out.bam 
    )
    ch_software_versions = ch_software_versions.mix(FILTER_BAM_SAMTOOLS.out.samtools_version.first().ifEmpty(null))

    /*
     * MODULE: VCF stats with bcftools stats
     */
    BCFTOOLS_STATS ( 
        ARTIC_MINION.out.vcf    
    )
    ch_software_versions = ch_software_versions.mix(BCFTOOLS_STATS.out.version.ifEmpty(null))

    /*
     * MODULE: Consensus QC across all samples with QUAST
     */
    ch_quast_multiqc = Channel.empty()
    if (!params.skip_variants_quast) {
        QUAST ( 
            ARTIC_MINION.out.fasta.collect{ it[1] },
            PREPARE_GENOME.out.fasta, 
            PREPARE_GENOME.out.gff, 
            true, 
            params.gff
        )
        ch_quast_multiqc = QUAST.out.tsv
        ch_software_versions = ch_software_versions.mix(QUAST.out.version.ifEmpty(null))
    }

    /*
     * MODULE: Lineage analysis with Pangolin
     */
    if (!params.skip_variants_pangolin) {
        PANGOLIN ( 
            ARTIC_MINION.out.fasta
        )
        ch_software_versions = ch_software_versions.mix(PANGOLIN.out.version.ifEmpty(null))
    }
    
    /*
     * SUBWORKFLOW: Annotate variants with snpEff
     */
    ch_snpeff_multiqc = Channel.empty()
    if (params.gff && !params.skip_variants_snpeff) {
        SNPEFF_SNPSIFT ( 
            ARTIC_MINION.out.vcf, 
            PREPARE_GENOME.out.snpeff_db, 
            PREPARE_GENOME.out.snpeff_config, 
            PREPARE_GENOME.out.fasta
        )
        ch_snpeff_multiqc = SNPEFF_SNPSIFT.out.csv
        ch_software_versions = ch_software_versions.mix(SNPEFF_SNPSIFT.out.snpeff_version.ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(SNPEFF_SNPSIFT.out.snpsift_version.ifEmpty(null))
    }

    /*
     * MODULE: Pipeline reporting
     */
    // Get unique list of files containing version information
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions
    )

    /*
     * MODULE: MultiQC
     */
    if (!params.skip_multiqc) {
        workflow_summary    = Schema.params_summary_multiqc(workflow, params.summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            GET_SOFTWARE_VERSIONS.out.yaml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            MULTIQC_CUSTOM_FAIL_NO_BARCODES.out.ifEmpty([]),
            MULTIQC_CUSTOM_FAIL_BARCODE_COUNT.out.ifEmpty([]),
            MULTIQC_CUSTOM_FAIL_GUPPYPLEX_COUNT.out.ifEmpty([]),
            FILTER_BAM_SAMTOOLS.out.stats.collect{it[1]}.ifEmpty([]),
            FILTER_BAM_SAMTOOLS.out.flagstat.collect{it[1]}.ifEmpty([]),
            FILTER_BAM_SAMTOOLS.out.idxstats.collect{it[1]}.ifEmpty([]),
            BCFTOOLS_STATS.out.stats.collect{it[1]}.ifEmpty([]),
            ch_quast_multiqc.collect().ifEmpty([]),
            ch_snpeff_multiqc.collect{it[1]}.ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////