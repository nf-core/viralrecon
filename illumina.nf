////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.fasta, params.gff, 
    params.bowtie2_index, params.amplicon_fasta, params.amplicon_bed,
    params.multiqc_config
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet file not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Genome fasta file not specified!' }
if (params.gff)   { ch_gff   = file(params.gff)   }

// Check optional parameters
def protocolList = ['metagenomic', 'amplicon']
if (!protocolList.contains(params.protocol)) {
    exit 1, "Invalid protocol option: ${params.protocol}. Valid options: ${protocolList.join(', ')}"
}

// Variant calling parameter validation
def callerList = [ 'varscan2', 'ivar', 'bcftools', 'none']
def callers = params.callers ? params.callers.split(',').collect{ it.trim().toLowerCase() } : []
if ((callerList + callers).unique().size() != callerList.size()) {
    exit 1, "Invalid variant calller option: ${params.callers}. Valid options: ${callerList.join(', ')}"
}

if (params.protocol == 'amplicon' && !params.skip_variants && !params.amplicon_bed) {
    exit 1, "To perform variant calling in 'amplicon' mode please provide a valid amplicon BED file!"
}
if (params.amplicon_bed) { ch_amplicon_bed = file(params.amplicon_bed) }

// Assembly parameter validation
def assemblerList = [ 'spades', 'metaspades', 'unicycler', 'minia', 'none' ]
def assemblers = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []
if ((assemblerList + assemblers).unique().size() != assemblerList.size()) {
    exit 1, "Invalid assembler option: ${params.assemblers}. Valid options: ${assemblerList.join(', ')}"
}

if (params.protocol == 'amplicon' && !params.skip_assembly && !params.amplicon_fasta) {
    exit 1, "To perform de novo assembly in 'amplicon' mode please provide a valid amplicon fasta file!"
}
if (params.amplicon_fasta) { ch_amplicon_fasta = file(params.amplicon_fasta) }

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

// Header files
ch_blast_outfmt6_header     = file("$projectDir/assets/headers/blast_outfmt6_header.txt", checkIfExists: true)
ch_ivar_variants_header_mqc = file("$projectDir/assets/headers/ivar_variants_header_mqc.txt", checkIfExists: true)

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def publish_genome_options = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options  = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]

// def cat_fastq_options          = modules['cat_fastq']
// if (!params.save_merged_fastq) { cat_fastq_options['publish_files'] = false }

// def multiqc_options         = modules['multiqc']
// multiqc_options.args       += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''
// if (params.skip_alignment)  { multiqc_options['publish_dir'] = '' }

// include { CAT_FASTQ                          } from './modules/local/process/cat_fastq'                   addParams( options: cat_fastq_options                                           ) 
// include { MULTIQC                            } from './modules/local/process/multiqc'                     addParams( options: multiqc_options                                             )
include { GET_SOFTWARE_VERSIONS              } from './modules/local/process/get_software_versions'       addParams( options: [publish_files : ['csv':'']]                                )

// /*
//  * SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//  */
def bowtie2_build_options    = modules['bowtie2_build']
// if (!params.save_reference) { bowtie2_build_options['publish_files'] = false }

// def bowtie2_align_options         = modules['bowtie2_align']
// if (params.save_align_intermeds) { bowtie2_align_options.publish_files.put('bam','') }
// if (params.save_unaligned)       { bowtie2_align_options.publish_files.put('fastq.gz','unmapped') }

// def samtools_sort_options = modules['samtools_sort']
// if (['star_salmon','bowtie2'].contains(params.aligner)) {
//     if (params.save_align_intermeds || (!params.with_umi && params.skip_markduplicates)) {
//         samtools_sort_options.publish_files.put('bam','')
//         samtools_sort_options.publish_files.put('bai','')
//     }
// } else {
//     if (params.save_align_intermeds || params.skip_markduplicates) {
//         samtools_sort_options.publish_files.put('bam','')
//         samtools_sort_options.publish_files.put('bai','')
//     }
// }
        
include { INPUT_CHECK     } from './modules/local/subworkflow/input_check'    addParams( options: [:] )
include { PREPARE_GENOME  } from './modules/local/subworkflow/prepare_genome' addParams( genome_options: publish_genome_options, index_options: publish_index_options, bowtie2_index_options: bowtie2_build_options )
// include { ALIGN_BOWTIE2    } from './modules/local/subworkflow/align_bowtie2'   addParams( align_options: bowtie2_align_options, samtools_options: samtools_sort_options )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

/*
 * MODULE: Installed directly from nf-core/modules
 */


/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */
// def trimgalore_options    = modules['trimgalore']
// trimgalore_options.args  += params.trim_nextseq > 0 ? " --nextseq ${params.trim_nextseq}" : ''
// if (params.save_trimmed)  { trimgalore_options.publish_files.put('fq.gz','') }

// include { FASTQC_FASTP           } from './modules/nf-core/subworkflow/fastqc_fastp' addParams( fastqc_options: modules['fastqc'], umitools_options: umitools_extract_options, trimgalore_options: trimgalore_options               )
// include { MARK_DUPLICATES_PICARD } from './modules/nf-core/subworkflow/mark_duplicates_picard'     addParams( markduplicates_options: modules['picard_markduplicates'], samtools_options: modules['picard_markduplicates_samtools'] )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report      = []
// def pass_percent_mapped = [:]
// def fail_percent_mapped = [:]

workflow ILLUMINA {

    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Uncompress and prepare reference genome files
     */
    PREPARE_GENOME (
        ch_dummy_file
    )
    
    
    // /*
    //  * SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //  */
    // INPUT_CHECK ( 
    //     ch_input
    // )
    // .map {
    //     meta, fastq ->
    //         meta.id = meta.id.split('_')[0..-2].join('_')
    //         [ meta, fastq ] }
    // .groupTuple(by: [0])
    // .map { it ->  [ it[0], it[1].flatten() ] }
    // .set { ch_cat_fastq }

    // /*
    //  * MODULE: Concatenate FastQ files from same sample if required
    //  */
    // CAT_FASTQ ( 
    //     ch_cat_fastq
    // )

    // /*
    //  * SUBWORKFLOW: Read QC, extract UMI and trim adapters
    //  */
    // FASTQC_UMITOOLS_TRIMGALORE (
    //     CAT_FASTQ.out.reads,
    //     params.skip_fastqc || params.skip_qc,
    //     params.with_umi,
    //     params.skip_trimming
    // )
    // ch_software_versions = ch_software_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_version.first().ifEmpty(null))
    // ch_software_versions = ch_software_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.umitools_version.first().ifEmpty(null))
    // ch_software_versions = ch_software_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.trimgalore_version.first().ifEmpty(null))

    // /*
    //  * MODULE: Remove ribosomal RNA reads
    //  */
    // ch_trimmed_reads     = FASTQC_UMITOOLS_TRIMGALORE.out.reads
    // ch_sortmerna_multiqc = Channel.empty()
    // if (params.remove_ribo_rna) {
    //     ch_sortmerna_fasta = Channel.from(ch_ribo_db.readLines()).map { row -> file(row) }.collect()

    //     SORTMERNA ( 
    //         ch_trimmed_reads, 
    //         ch_sortmerna_fasta
    //     )
    //     .reads
    //     .set { ch_trimmed_reads }

    //     ch_sortmerna_multiqc = SORTMERNA.out.log
    //     ch_software_versions = ch_software_versions.mix(SORTMERNA.out.version.first().ifEmpty(null))
    // }

    // /*
    //  * SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    //  */
    
    // /*
    //  * SUBWORKFLOW: Alignment with BOWTIE2
    //  */
    // ch_genome_bam        = Channel.empty()
    // ch_genome_bai        = Channel.empty()
    // ch_samtools_stats    = Channel.empty()
    // ch_samtools_flagstat = Channel.empty()
    // ch_samtools_idxstats = Channel.empty()
    // ch_bowtie2_multiqc    = Channel.empty()
    // if (!params.skip_alignment && params.aligner == 'bowtie2') {        
    //     ALIGN_BOWTIE2 (
    //         ch_trimmed_reads,
    //         PREPARE_GENOME.out.bowtie2_index,
    //         PREPARE_GENOME.out.splicesites
    //     )
    //     ch_genome_bam        = ALIGN_BOWTIE2.out.bam
    //     ch_genome_bai        = ALIGN_BOWTIE2.out.bai
    //     ch_samtools_stats    = ALIGN_BOWTIE2.out.stats
    //     ch_samtools_flagstat = ALIGN_BOWTIE2.out.flagstat
    //     ch_samtools_idxstats = ALIGN_BOWTIE2.out.idxstats
    //     ch_bowtie2_multiqc    = ALIGN_BOWTIE2.out.summary
    //     ch_software_versions = ch_software_versions.mix(ALIGN_BOWTIE2.out.bowtie2_version.first().ifEmpty(null))
    //     ch_software_versions = ch_software_versions.mix(ALIGN_BOWTIE2.out.samtools_version.first().ifEmpty(null))
    // }

    // /*
    //  * Filter channels to get samples that passed STAR minimum mapping percentage
    //  */
    // ch_fail_mapping_multiqc = Channel.empty()
    // if (!params.skip_alignment && params.aligner.contains('star')) {
    //     ch_star_multiqc
    //         .map { meta, align_log -> [ meta ] + Checks.get_star_percent_mapped(workflow, params, log, align_log) }
    //         .set { ch_percent_mapped }

    //     ch_genome_bam
    //         .join(ch_percent_mapped, by: [0])
    //         .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
    //         .set { ch_genome_bam }

    //     ch_genome_bai
    //         .join(ch_percent_mapped, by: [0])
    //         .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
    //         .set { ch_genome_bai }

    //     ch_percent_mapped
    //         .branch { meta, mapped, pass ->
    //             pass: pass
    //                 pass_percent_mapped[meta.id] = mapped
    //                 return [ "$meta.id\t$mapped" ]
    //             fail: !pass
    //                 fail_percent_mapped[meta.id] = mapped
    //                 return [ "$meta.id\t$mapped" ]
    //         }
    //         .set { ch_pass_fail_mapped }

    //     MULTIQC_CUSTOM_FAIL_MAPPED ( 
    //         ch_pass_fail_mapped.fail.collect()
    //     )
    //     .set { ch_fail_mapping_multiqc }
    // }

    // /*
    //  * SUBWORKFLOW: Mark duplicate reads
    //  */
    // ch_markduplicates_multiqc = Channel.empty()
    // if (!params.skip_alignment && !params.skip_markduplicates) {
    //     MARK_DUPLICATES_PICARD (
    //         ch_genome_bam
    //     )
    //     ch_genome_bam             = MARK_DUPLICATES_PICARD.out.bam
    //     ch_genome_bai             = MARK_DUPLICATES_PICARD.out.bai
    //     ch_samtools_stats         = MARK_DUPLICATES_PICARD.out.stats
    //     ch_samtools_flagstat      = MARK_DUPLICATES_PICARD.out.flagstat
    //     ch_samtools_idxstats      = MARK_DUPLICATES_PICARD.out.idxstats
    //     ch_markduplicates_multiqc = MARK_DUPLICATES_PICARD.out.metrics
    //     ch_software_versions      = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.picard_version.first().ifEmpty(null))
    // }

    /*
     * MODULE: Pipeline reporting
     */
    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions.map { it }.collect()
    )

    // /*
    //  * MultiQC
    //  */
    // if (!params.skip_multiqc) {
    //     workflow_summary    = Schema.params_summary_multiqc(workflow, params.summary_params)
    //     ch_workflow_summary = Channel.value(workflow_summary)

    //     MULTIQC (
    //         ch_multiqc_config,
    //         ch_multiqc_custom_config.collect().ifEmpty([]),
    //         GET_SOFTWARE_VERSIONS.out.yaml.collect(),
    //         ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
    //         ch_fail_mapping_multiqc.ifEmpty([]),
    //         ch_fail_strand_multiqc.ifEmpty([]),
    //         FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
    //         FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]),
    //         FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),
    //         ch_sortmerna_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_star_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_bowtie2_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_rsem_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_salmon_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_samtools_stats.collect{it[1]}.ifEmpty([]),
    //         ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
    //         ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
    //         ch_markduplicates_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_featurecounts_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_aligner_pca_multiqc.collect().ifEmpty([]),
    //         ch_aligner_clustering_multiqc.collect().ifEmpty([]),
    //         ch_pseudoaligner_pca_multiqc.collect().ifEmpty([]),
    //         ch_pseudoaligner_clustering_multiqc.collect().ifEmpty([]),
    //         ch_preseq_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_qualimap_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_dupradar_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_bamstat_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_inferexperiment_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_innerdistance_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_junctionannotation_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_junctionsaturation_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_readdistribution_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_readduplication_multiqc.collect{it[1]}.ifEmpty([])
    //     )
    //     multiqc_report = MULTIQC.out.report.toList()
    // }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

// workflow.onComplete {
//     Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report, fail_percent_mapped)
//     Completion.summary(workflow, params, log, fail_percent_mapped, pass_percent_mapped)
// }

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////

// // Print warning if viral genome fasta has more than one sequence
// def count = 0
// ch_fasta.withReader { reader ->
//     while (line = reader.readLine()) {
//         if (line.contains('>')) {
//             count++
//             if (count > 1) {
//                 log.info "[nf-core/viralrecon] WARNING: This pipeline does not support multi-fasta genome files. Please amend the '--fasta' parameter."
//                 break
//             }
//         }
//     }
// }

// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                 UNZIP/UNTAR REFERENCE FILES                         -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////

// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                     PARSE DESIGN FILE                               -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////

// /*
//  * PREPROCESSING: Reformat samplesheet and check validity
//  */
// process CHECK_SAMPLESHEET {
//     tag "$samplesheet"
//     publishDir "${params.outdir}/", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".tsv")) "preprocess/sra/$filename"
//                       else "pipeline_info/$filename"
//                 }

//     input:
//     path samplesheet from ch_input

//     output:
//     path "samplesheet.valid.csv" into ch_samplesheet_reformat
//     path "sra_run_info.tsv" optional true

//     script:  // These scripts are bundled with the pipeline, in nf-core/viralrecon/bin/
//     run_sra = !params.skip_sra && !isOffline()
//     """
//     awk -F, '{if(\$1 != "" && \$2 != "") {print \$0}}' $samplesheet > nonsra_id.csv
//     check_samplesheet.py nonsra_id.csv nonsra.samplesheet.csv

//     awk -F, '{if(\$1 != "" && \$2 == "" && \$3 == "") {print \$1}}' $samplesheet > sra_id.list
//     if $run_sra && [ -s sra_id.list ]
//     then
//         fetch_sra_runinfo.py sra_id.list sra_run_info.tsv --platform ILLUMINA --library_layout SINGLE,PAIRED
//         sra_runinfo_to_samplesheet.py sra_run_info.tsv sra.samplesheet.csv
//     fi

//     if [ -f nonsra.samplesheet.csv ]
//     then
//         head -n 1 nonsra.samplesheet.csv > samplesheet.valid.csv
//     else
//         head -n 1 sra.samplesheet.csv > samplesheet.valid.csv
//     fi
//     tail -n +2 -q *sra.samplesheet.csv >> samplesheet.valid.csv
//     """
// }

// // Function to get list of [ sample, single_end?, is_sra?, is_ftp?, [ fastq_1, fastq_2 ], [ md5_1, md5_2] ]
// def validate_input(LinkedHashMap sample) {
//     def sample_id = sample.sample_id
//     def single_end = sample.single_end.toBoolean()
//     def is_sra = sample.is_sra.toBoolean()
//     def is_ftp = sample.is_ftp.toBoolean()
//     def fastq_1 = sample.fastq_1
//     def fastq_2 = sample.fastq_2
//     def md5_1 = sample.md5_1
//     def md5_2 = sample.md5_2

//     def array = []
//     if (!is_sra) {
//         if (single_end) {
//             array = [ sample_id, single_end, is_sra, is_ftp, [ file(fastq_1, checkIfExists: true) ] ]
//         } else {
//             array = [ sample_id, single_end, is_sra, is_ftp, [ file(fastq_1, checkIfExists: true), file(fastq_2, checkIfExists: true) ] ]
//         }
//     } else {
//         array = [ sample_id, single_end, is_sra, is_ftp, [ fastq_1, fastq_2 ], [ md5_1, md5_2 ] ]
//     }

//     return array
// }

// /*
//  * Create channels for input fastq files
//  */
// ch_samplesheet_reformat
//     .splitCsv(header:true, sep:',')
//     .map { validate_input(it) }
//     .into { ch_reads_all
//             ch_reads_sra }

// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                     MERGE RESEQUENCED FASTQ                         -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////

// /*
//  * STEP 2: Merge FastQ files with the same sample identifier
//  */
// process CAT_FASTQ {
//     tag "$sample"

//     input:
//     tuple val(sample), val(single_end), path(reads) from ch_reads_all

//     output:
//     tuple val(sample), val(single_end), path("*.merged.fastq.gz") into ch_cat_fastqc,
//                                                                        ch_cat_fastp

//     script:
//     readList = reads.collect{it.toString()}
//     if (!single_end) {
//         if (readList.size > 2) {
//             def read1 = []
//             def read2 = []
//             readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
//             """
//             cat ${read1.sort().join(' ')} > ${sample}_1.merged.fastq.gz
//             cat ${read2.sort().join(' ')} > ${sample}_2.merged.fastq.gz
//             """
//         } else {
//             """
//             ln -s ${reads[0]} ${sample}_1.merged.fastq.gz
//             ln -s ${reads[1]} ${sample}_2.merged.fastq.gz
//             """
//         }
//     } else {
//         if (readList.size > 1) {
//             """
//             cat ${readList.sort().join(' ')} > ${sample}.merged.fastq.gz
//             """
//         } else {
//             """
//             ln -s $reads ${sample}.merged.fastq.gz
//             """
//         }
//     }
// }

// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                        FASTQ QC                                     -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////

// /*
//  * STEP 3: FastQC on input reads after merging libraries from the same sample
//  */
// process FASTQC {
//     tag "$sample"
//     label 'process_medium'
//     publishDir "${params.outdir}/preprocess/fastqc", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       filename.endsWith(".zip") ? "zips/$filename" : filename
//                 }

//     when:
//     !params.skip_fastqc

//     input:
//     tuple val(sample), val(single_end), path(reads) from ch_cat_fastqc

//     output:
//     path "*.{zip,html}" into ch_fastqc_raw_reports_mqc

//     script:
//     """
//     fastqc --quiet --threads $task.cpus *.fastq.gz
//     """
// }

// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                        ADAPTER TRIMMING                             -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////

// /*
//  * STEP 4: Fastp adapter trimming and quality filtering
//  */
// if (!params.skip_adapter_trimming) {
//     process FASTP {
//         tag "$sample"
//         label 'process_medium'
//         publishDir "${params.outdir}/preprocess/fastp", mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                           if (filename.endsWith(".json")) filename
//                           else if (filename.endsWith(".fastp.html")) filename
//                           else if (filename.endsWith("_fastqc.html")) "fastqc/$filename"
//                           else if (filename.endsWith(".zip")) "fastqc/zips/$filename"
//                           else if (filename.endsWith(".log")) "log/$filename"
//                           else params.save_trimmed ? filename : null
//                     }

//         when:
//         !params.skip_variants || !params.skip_assembly

//         input:
//         tuple val(sample), val(single_end), path(reads) from ch_cat_fastp

//         output:
//         tuple val(sample), val(single_end), path("*.trim.fastq.gz") into ch_fastp_bowtie2,
//                                                                          ch_fastp_cutadapt,
//                                                                          ch_fastp_kraken2
//         path "*.json" into ch_fastp_mqc
//         path "*_fastqc.{zip,html}" into ch_fastp_fastqc_mqc
//         path "*.{log,fastp.html}"
//         path "*.fail.fastq.gz"

//         script:
//         // Added soft-links to original fastqs for consistent naming in MultiQC
//         autodetect = single_end ? "" : "--detect_adapter_for_pe"
//         """
//         IN_READS='--in1 ${sample}.fastq.gz'
//         OUT_READS='--out1 ${sample}.trim.fastq.gz --failed_out ${sample}.fail.fastq.gz'
//         if $single_end; then
//             [ ! -f  ${sample}.fastq.gz ] && ln -s $reads ${sample}.fastq.gz
//         else
//             [ ! -f  ${sample}_1.fastq.gz ] && ln -s ${reads[0]} ${sample}_1.fastq.gz
//             [ ! -f  ${sample}_2.fastq.gz ] && ln -s ${reads[1]} ${sample}_2.fastq.gz
//             IN_READS='--in1 ${sample}_1.fastq.gz --in2 ${sample}_2.fastq.gz'
//             OUT_READS='--out1 ${sample}_1.trim.fastq.gz --out2 ${sample}_2.trim.fastq.gz --unpaired1 ${sample}_1.fail.fastq.gz --unpaired2 ${sample}_2.fail.fastq.gz'
//         fi

//         fastp \\
//             \$IN_READS \\
//             \$OUT_READS \\
//             $autodetect \\
//             --cut_front \\
//             --cut_tail \\
//             --cut_mean_quality $params.cut_mean_quality \\
//             --qualified_quality_phred $params.qualified_quality_phred \\
//             --unqualified_percent_limit $params.unqualified_percent_limit \\
//             --length_required $params.min_trim_length \\
//             --trim_poly_x \\
//             --thread $task.cpus \\
//             --json ${sample}.fastp.json \\
//             --html ${sample}.fastp.html \\
//             2> ${sample}.fastp.log

//         fastqc --quiet --threads $task.cpus *.trim.fastq.gz
//         """
//     }
// } else {
//     ch_cat_fastp
//         .into { ch_fastp_bowtie2
//                 ch_fastp_cutadapt
//                 ch_fastp_kraken2 }
//     ch_fastp_mqc = Channel.empty()
//     ch_fastp_fastqc_mqc = Channel.empty()
// }

// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                  VARIANT CALLING PROCESSES                          -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////

// /*
//  * PREPROCESSING: Build Bowtie2 index for viral genome
//  */
// process BOWTIE2_INDEX {
//     tag "$fasta"
//     label 'process_medium'
//     if (params.save_reference) {
//         publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
//     }

//     when:
//     !params.skip_variants

//     input:
//     path fasta from ch_fasta

//     output:
//     path "Bowtie2Index" into ch_index

//     script:
//     """
//     bowtie2-build \\
//         --seed 1 \\
//         --threads $task.cpus \\
//         $fasta \\
//         $index_base
//     mkdir Bowtie2Index && mv ${index_base}* Bowtie2Index
//     """
// }

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

// /*
//  * STEP 5.1: Map read(s) with Bowtie 2
//  */
// process BOWTIE2 {
//     tag "$sample"
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/bam", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".log")) "log/$filename"
//                       else params.save_align_intermeds ? filename : null
//                 }

//     when:
//     !params.skip_variants

//     input:
//     tuple val(sample), val(single_end), path(reads) from ch_fastp_bowtie2
//     path index from ch_index

//     output:
//     tuple val(sample), val(single_end), path("*.bam") into ch_bowtie2_bam
//     path "*.log" into ch_bowtie2_mqc

//     script:
//     input_reads = single_end ? "-U $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
//     filter = params.filter_unmapped ? "-F4" : ""
//     """
//     bowtie2 \\
//         --threads $task.cpus \\
//         --local \\
//         --very-sensitive-local \\
//         -x ${index}/${index_base} \\
//         $input_reads \\
//         2> ${sample}.bowtie2.log \\
//         | samtools view -@ $task.cpus -b -h -O BAM -o ${sample}.bam $filter -
//     """
// }

// /*
//  * STEP 5.2: Convert BAM to coordinate sorted BAM
//  */
// process SORT_BAM {
//     tag "$sample"
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/bam", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
//                       else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
//                       else if (filename.endsWith(".stats")) "samtools_stats/$filename"
//                       else (params.protocol != 'amplicon' && params.skip_markduplicates) || params.save_align_intermeds ? filename : null
//                 }

//     when:
//     !params.skip_variants

//     input:
//     tuple val(sample), val(single_end), path(bam) from ch_bowtie2_bam

//     output:
//     tuple val(sample), val(single_end), path("*.sorted.{bam,bam.bai}"), path("*.flagstat") into ch_sort_bam
//     path "*.{flagstat,idxstats,stats}" into ch_sort_bam_flagstat_mqc

//     script:
//     """
//     samtools sort -@ $task.cpus -o ${sample}.sorted.bam -T $sample $bam
//     samtools index ${sample}.sorted.bam
//     samtools flagstat ${sample}.sorted.bam > ${sample}.sorted.bam.flagstat
//     samtools idxstats ${sample}.sorted.bam > ${sample}.sorted.bam.idxstats
//     samtools stats ${sample}.sorted.bam > ${sample}.sorted.bam.stats
//     """
// }

// // Get total number of mapped reads from flagstat file
// def get_mapped_from_flagstat(flagstat) {
//     def mapped = 0
//     flagstat.eachLine { line ->
//         if (line.contains(' mapped (')) {
//             mapped = line.tokenize().first().toInteger()
//         }
//     }
//     return mapped
// }

// // Function that checks the number of mapped reads from flagstat output
// // and returns true if > params.min_mapped_reads and otherwise false
// pass_mapped_reads = [:]
// fail_mapped_reads = [:]
// def check_mapped(sample,flagstat,min_mapped_reads=500) {
//     mapped = get_mapped_from_flagstat(flagstat)
//     c_reset = params.monochrome_logs ? '' : "\033[0m";
//     c_green = params.monochrome_logs ? '' : "\033[0;32m";
//     c_red = params.monochrome_logs ? '' : "\033[0;31m";
//     if (mapped < min_mapped_reads.toInteger()) {
//         log.info ">${c_red}>>>> $sample FAILED MAPPED READ THRESHOLD: ${mapped} < ${params.min_mapped_reads}. IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! <<<<${c_reset}<"
//         fail_mapped_reads[sample] = mapped
//         return false
//     } else {
//         //log.info "-${c_green}           Passed mapped read threshold > bowtie2 ($sample)   >> ${mapped} <<${c_reset}"
//         pass_mapped_reads[sample] = mapped
//         return true
//     }
// }

// // Remove samples that failed mapped read threshold
// ch_sort_bam
//     .filter { sample, single_end, bam, flagstat -> check_mapped(sample,flagstat,params.min_mapped_reads) }
//     .map { it[0..2] }
//     .set { ch_sort_bam }

// /*
//  * STEP 5.3: Trim amplicon sequences with iVar
//  */
// if (params.protocol != 'amplicon') {
//     ch_sort_bam
//         .set { ch_ivar_trim_bam }
//     ch_ivar_trim_flagstat_mqc = Channel.empty()
//     ch_ivar_trim_log_mqc = Channel.empty()
// } else {
//     process IVAR_TRIM {
//         tag "$sample"
//         label 'process_medium'
//         publishDir "${params.outdir}/variants/bam", mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                           if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
//                           else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
//                           else if (filename.endsWith(".stats")) "samtools_stats/$filename"
//                           else if (filename.endsWith(".log")) "log/$filename"
//                           else params.skip_markduplicates || params.save_align_intermeds ? filename : null
//                     }

//         when:
//         !params.skip_variants

//         input:
//         tuple val(sample), val(single_end), path(bam) from ch_sort_bam
//         path bed from ch_amplicon_bed

//         output:
//         tuple val(sample), val(single_end), path("*.sorted.{bam,bam.bai}") into ch_ivar_trim_bam
//         path "*.{flagstat,idxstats,stats}" into ch_ivar_trim_flagstat_mqc
//         path "*.log" into ch_ivar_trim_log_mqc

//         script:
//         exclude_reads = params.ivar_trim_noprimer ? "" : "-e"
//         prefix = "${sample}.trim"
//         """
//         samtools view -b -F 4 ${bam[0]} > ${sample}.mapped.bam
//         samtools index ${sample}.mapped.bam

//         ivar trim \\
//             -i ${sample}.mapped.bam \\
//             -b $bed \\
//             -m $params.ivar_trim_min_len \\
//             -q $params.ivar_trim_min_qual \\
//             -s $params.ivar_trim_window_width \\
//             $exclude_reads \\
//             -p $prefix > ${prefix}.ivar.log

//         samtools sort -@ $task.cpus -o ${prefix}.sorted.bam -T $prefix ${prefix}.bam
//         samtools index ${prefix}.sorted.bam
//         samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
//         samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
//         samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
//         """
//     }
// }

// /*
//  * STEP 5.4: Picard MarkDuplicates
//  */
// if (params.skip_markduplicates) {
//     ch_ivar_trim_bam
//         .into { ch_markdup_bam_metrics
//                 ch_markdup_bam_mosdepth_genome
//                 ch_markdup_bam_mosdepth_amplicon
//                 ch_markdup_bam_mpileup
//                 ch_markdup_bam_varscan2_consensus
//                 ch_markdup_bam_bcftools
//                 ch_markdup_bam_bcftools_consensus }
//     ch_markdup_bam_flagstat_mqc = Channel.empty()
//     ch_markdup_bam_metrics_mqc = Channel.empty()
// } else {
//     process PICARD_MARKDUPLICATES {
//         tag "$sample"
//         label 'process_medium'
//         publishDir "${params.outdir}/variants/bam", mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                           if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
//                           else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
//                           else if (filename.endsWith(".stats")) "samtools_stats/$filename"
//                           else if (filename.endsWith(".metrics.txt")) "picard_metrics/$filename"
//                           else filename
//                     }

//         when:
//         !params.skip_variants

//         input:
//         tuple val(sample), val(single_end), path(bam) from ch_ivar_trim_bam
//         path fasta from ch_fasta

//         output:
//         tuple val(sample), val(single_end), path("*.sorted.{bam,bam.bai}") into ch_markdup_bam_metrics,
//                                                                                 ch_markdup_bam_mosdepth_genome,
//                                                                                 ch_markdup_bam_mosdepth_amplicon,
//                                                                                 ch_markdup_bam_mpileup,
//                                                                                 ch_markdup_bam_varscan2_consensus,
//                                                                                 ch_markdup_bam_bcftools,
//                                                                                 ch_markdup_bam_bcftools_consensus
//         path "*.{flagstat,idxstats,stats}" into ch_markdup_bam_flagstat_mqc
//         path "*.txt" into ch_markdup_bam_metrics_mqc

//         script:
//         def avail_mem = 3
//         if (!task.memory) {
//             log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
//         } else {
//             avail_mem = task.memory.toGiga()
//         }
//         prefix = params.protocol == 'amplicon' ? "${sample}.trim.mkD" : "${sample}.mkD"
//         keep_dup = params.filter_dups ? "true" : "false"
//         """
//         picard -Xmx${avail_mem}g MarkDuplicates \\
//             INPUT=${bam[0]} \\
//             OUTPUT=${prefix}.sorted.bam \\
//             ASSUME_SORTED=true \\
//             REMOVE_DUPLICATES=$keep_dup \\
//             METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt \\
//             VALIDATION_STRINGENCY=LENIENT \\
//             TMP_DIR=tmp
//         samtools index ${prefix}.sorted.bam
//         samtools idxstats ${prefix}.sorted.bam > ${prefix}.sorted.bam.idxstats
//         samtools flagstat ${prefix}.sorted.bam > ${prefix}.sorted.bam.flagstat
//         samtools stats ${prefix}.sorted.bam > ${prefix}.sorted.bam.stats
//         """
//     }
// }

// /*
//  * STEP 5.5: Picard CollectMultipleMetrics and CollectWgsMetrics
//  */
// process PICARD_METRICS {
//     tag "$sample"
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/bam/picard_metrics", mode: params.publish_dir_mode

//     when:
//     !params.skip_variants && !params.skip_picard_metrics

//     input:
//     tuple val(sample), val(single_end), path(bam) from ch_markdup_bam_metrics
//     path fasta from ch_fasta

//     output:
//     path "*metrics" into ch_picard_metrics_mqc
//     path "*.pdf"

//     script:
//     def avail_mem = 3
//     if (!task.memory) {
//         log.info "[Picard CollectMultipleMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
//     } else {
//         avail_mem = task.memory.toGiga()
//     }
//     suffix = params.skip_markduplicates ? "" : ".mkD"
//     prefix = params.protocol == 'amplicon' ? "${sample}.trim${suffix}" : "${sample}${suffix}"
//     """
//     picard -Xmx${avail_mem}g CollectMultipleMetrics \\
//         INPUT=${bam[0]} \\
//         OUTPUT=${prefix}.CollectMultipleMetrics \\
//         REFERENCE_SEQUENCE=$fasta \\
//         VALIDATION_STRINGENCY=LENIENT \\
//         TMP_DIR=tmp

//     picard -Xmx${avail_mem}g CollectWgsMetrics \\
//         COVERAGE_CAP=1000000 \\
//         INPUT=${bam[0]} \\
//         OUTPUT=${prefix}.CollectWgsMetrics.coverage_metrics \\
//         REFERENCE_SEQUENCE=$fasta \\
//         VALIDATION_STRINGENCY=LENIENT \\
//         TMP_DIR=tmp
//     """
// }

// /*
//  * STEP 5.6.1: mosdepth genome-wide coverage
//  */
// process MOSDEPTH_GENOME {
//     tag "$sample"
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/bam/mosdepth/genome", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".pdf")) "plots/$filename"
//                       else if (filename.endsWith(".tsv")) "plots/$filename"
//                       else filename
//                 }

//     when:
//     !params.skip_variants && !params.skip_mosdepth

//     input:
//     tuple val(sample), val(single_end), path(bam) from ch_markdup_bam_mosdepth_genome

//     output:
//     path "*.global.dist.txt" into ch_mosdepth_genome_mqc
//     path "*.{txt,gz,csi,tsv,pdf}"

//     script:
//     suffix = params.skip_markduplicates ? "" : ".mkD"
//     prefix = params.protocol == 'amplicon' ? "${sample}.trim${suffix}.genome" : "${sample}${suffix}.genome"
//     plot_suffix = params.protocol == 'amplicon' ? ".trim${suffix}.genome" : "${suffix}.genome"
//     """
//     mosdepth \\
//         --by 200 \\
//         --fast-mode \\
//         $prefix \\
//         ${bam[0]}

//     plot_mosdepth_regions.r \\
//         --input_files ${prefix}.regions.bed.gz \\
//         --input_suffix ${plot_suffix}.regions.bed.gz \\
//         --output_dir ./ \\
//         --output_suffix ${plot_suffix}.regions
//     """
// }

// /*
//  * STEP 5.6.2: mosdepth amplicon coverage and plots
//  */
// if (params.protocol == 'amplicon') {
//     process MOSDEPTH_AMPLICON {
//         tag "$sample"
//         label 'process_medium'
//         publishDir "${params.outdir}/variants/bam/mosdepth/amplicon", mode: params.publish_dir_mode

//         when:
//         !params.skip_variants && !params.skip_mosdepth

//         input:
//         tuple val(sample), val(single_end), path(bam) from ch_markdup_bam_mosdepth_amplicon
//         path bed from ch_amplicon_bed

//         output:
//         path "*.regions.bed.gz" into ch_mosdepth_amplicon_region_bed
//         path "*.{txt,gz,csi}"

//         script:
//         suffix = params.skip_markduplicates ? "" : ".mkD"
//         prefix = "${sample}.trim${suffix}.amplicon"
//         """
//         collapse_amplicon_bed.py \\
//             --left_primer_suffix $params.amplicon_left_suffix \\
//             --right_primer_suffix $params.amplicon_right_suffix \\
//             $bed \\
//             amplicon.collapsed.bed

//         mosdepth \\
//             --by amplicon.collapsed.bed \\
//             --fast-mode \\
//             --use-median \\
//             --thresholds 0,1,10,50,100,500 \\
//             ${prefix} \\
//             ${bam[0]}
//         """
//     }

//     process MOSDEPTH_AMPLICON_PLOT {
//         label 'process_medium'
//         publishDir "${params.outdir}/variants/bam/mosdepth/amplicon/plots", mode: params.publish_dir_mode

//         when:
//         !params.skip_variants && !params.skip_mosdepth

//         input:
//         path bed from ch_mosdepth_amplicon_region_bed.collect()

//         output:
//         path "*.{tsv,pdf}"

//         script:
//         suffix = params.skip_markduplicates ? "" : ".mkD"
//         suffix = ".trim${suffix}.amplicon"
//         """
//         plot_mosdepth_regions.r \\
//             --input_files ${bed.join(',')} \\
//             --input_suffix ${suffix}.regions.bed.gz \\
//             --output_dir ./ \\
//             --output_suffix ${suffix}.regions
//         """
//     }
// }

// ////////////////////////////////////////////////////
// /* --              VARSCAN2                    -- */
// ////////////////////////////////////////////////////

// /*
//  * STEP 5.7: Create mpileup file for all variant callers
//  */
// process SAMTOOLS_MPILEUP {
//     tag "$sample"
//     label 'process_medium'
//     if (params.save_mpileup) {
//         publishDir "${params.outdir}/variants/bam/mpileup", mode: params.publish_dir_mode
//     }

//     when:
//     !params.skip_variants

//     input:
//     tuple val(sample), val(single_end), path(bam) from ch_markdup_bam_mpileup
//     path fasta from ch_fasta

//     output:
//     tuple val(sample), val(single_end), path("*.mpileup") into ch_mpileup_varscan2,
//                                                                ch_mpileup_ivar_variants,
//                                                                ch_mpileup_ivar_consensus,
//                                                                ch_mpileup_ivar_bcftools

//     script:
//     suffix = params.skip_markduplicates ? "" : ".mkD"
//     prefix = params.protocol == 'amplicon' ? "${sample}.trim${suffix}" : "${sample}${suffix}"
//     """
//     samtools mpileup \\
//         --count-orphans \\
//         --ignore-overlaps \\
//         --no-BAQ \\
//         --max-depth $params.mpileup_depth \\
//         --fasta-ref $fasta \\
//         --min-BQ $params.min_base_qual \\
//         --output ${prefix}.mpileup \\
//         ${bam[0]}
//     """
// }

// /*
//  * STEP 5.7.1: Variant calling with VarScan 2
//  */
// process VARSCAN2 {
//     tag "$sample"
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/varscan2", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".log")) "log/$filename"
//                       else if (filename.endsWith(".txt")) "bcftools_stats/$filename"
//                       else filename
//                 }

//     when:
//     !params.skip_variants && 'varscan2' in callers

//     input:
//     tuple val(sample), val(single_end), path(mpileup) from ch_mpileup_varscan2
//     path fasta from ch_fasta

//     output:
//     tuple val(sample), val(single_end), path("${prefix}.vcf.gz*") into ch_varscan2_highfreq_consensus,
//                                                                        ch_varscan2_highfreq_snpeff,
//                                                                        ch_varscan2_highfreq_intersect
//     tuple val(sample), val(single_end), path("${sample}.vcf.gz*") into ch_varscan2_lowfreq_snpeff
//     path "${prefix}.bcftools_stats.txt" into ch_varscan2_bcftools_highfreq_mqc
//     path "*.varscan2.log" into ch_varscan2_log_mqc
//     path "${sample}.bcftools_stats.txt"

//     script:
//     prefix = "${sample}.AF${params.max_allele_freq}"
//     strand = params.protocol != 'amplicon' && params.varscan2_strand_filter ? "--strand-filter 1" : "--strand-filter 0"
//     """
//     echo "$sample" > sample_name.list
//     varscan mpileup2cns \\
//         $mpileup \\
//         --min-coverage $params.min_coverage \\
//         --min-reads2 5 \\
//         --min-avg-qual $params.min_base_qual \\
//         --min-var-freq $params.min_allele_freq \\
//         --p-value 0.99 \\
//         --output-vcf 1 \\
//         --vcf-sample-list sample_name.list \\
//         --variants \\
//         $strand \\
//         2> ${sample}.varscan2.log \\
//         | bgzip -c > ${sample}.vcf.gz
//     tabix -p vcf -f ${sample}.vcf.gz
//     bcftools stats ${sample}.vcf.gz > ${sample}.bcftools_stats.txt
//     sed -i.bak '/LC_ALL/d' ${sample}.varscan2.log

//     bcftools filter \\
//         -i 'FORMAT/AD / (FORMAT/AD + FORMAT/RD) >= $params.max_allele_freq' \\
//         --output-type z \\
//         --output ${prefix}.vcf.gz \\
//         ${sample}.vcf.gz
//     tabix -p vcf -f ${prefix}.vcf.gz
//     bcftools stats ${prefix}.vcf.gz > ${prefix}.bcftools_stats.txt
//     """
// }

// /*
//  * STEP 5.7.1.1: Genome consensus generation with BCFtools and masked with BEDTools
//  */
// process VARSCAN2_CONSENSUS {
//     tag "$sample"
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/varscan2/consensus", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".tsv")) "base_qc/$filename"
//                       else if (filename.endsWith(".pdf")) "base_qc/$filename"
//                       else filename
//                 }

//     when:
//     !params.skip_variants && 'varscan2' in callers

//     input:
//     tuple val(sample), val(single_end), path(bam), path(vcf) from ch_markdup_bam_varscan2_consensus.join(ch_varscan2_highfreq_consensus, by: [0,1])
//     path fasta from ch_fasta

//     output:
//     tuple val(sample), val(single_end), path("*consensus.masked.fa") into ch_varscan2_consensus
//     path "*.{consensus.fa,tsv,pdf}"

//     script:
//     prefix = "${sample}.AF${params.max_allele_freq}"
//     """
//     bedtools genomecov \\
//         -bga \\
//         -ibam ${bam[0]} \\
//         -g $fasta \\
//         | awk '\$4 < $params.min_coverage' > ${prefix}.lowcov.bed

//     parse_mask_bed.py ${vcf[0]} ${prefix}.lowcov.bed ${prefix}.lowcov.fix.bed

//     bedtools merge -i ${prefix}.lowcov.fix.bed > ${prefix}.mask.bed

//     bedtools maskfasta \\
//         -fi $fasta \\
//         -bed ${prefix}.mask.bed \\
//         -fo ${index_base}.ref.masked.fa

//     cat ${index_base}.ref.masked.fa | bcftools consensus ${vcf[0]} > ${prefix}.consensus.masked.fa

//     header=\$(head -n 1 ${prefix}.consensus.masked.fa | sed 's/>//g')
//     sed -i "s/\${header}/${sample}/g" ${prefix}.consensus.masked.fa

//     plot_base_density.r --fasta_files ${prefix}.consensus.masked.fa --prefixes $prefix --output_dir ./
//     """
// }

// /*
//  * STEP 5.7.1.2: VarScan 2 variant calling annotation with SnpEff and SnpSift
//  */
// process VARSCAN2_SNPEFF {
//     tag "$sample"
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/varscan2/snpeff", mode: params.publish_dir_mode

//     when:
//     !params.skip_variants && 'varscan2' in callers && params.gff && !params.skip_snpeff

//     input:
//     tuple val(sample), val(single_end), path(highfreq_vcf), path(lowfreq_vcf) from ch_varscan2_highfreq_snpeff.join(ch_varscan2_lowfreq_snpeff, by: [0,1])
//     tuple file(db), file(config) from ch_snpeff_db_varscan2

//     output:
//     path "${prefix}.snpEff.csv" into ch_varscan2_snpeff_highfreq_mqc
//     path "${sample}.snpEff.csv"
//     path "*.vcf.gz*"
//     path "*.{txt,html}"

//     script:
//     prefix = "${sample}.AF${params.max_allele_freq}"
//     """
//     snpEff ${index_base} \\
//         -config $config \\
//         -dataDir $db \\
//         ${lowfreq_vcf[0]} \\
//         -csvStats ${sample}.snpEff.csv \\
//         | bgzip -c > ${sample}.snpEff.vcf.gz
//     tabix -p vcf -f ${sample}.snpEff.vcf.gz
//     mv snpEff_summary.html ${sample}.snpEff.summary.html

//     SnpSift extractFields -s "," \\
//         -e "." \\
//         ${sample}.snpEff.vcf.gz \\
//         CHROM POS REF ALT \\
//         "ANN[*].GENE" "ANN[*].GENEID" \\
//         "ANN[*].IMPACT" "ANN[*].EFFECT" \\
//         "ANN[*].FEATURE" "ANN[*].FEATUREID" \\
//         "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \\
//         "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \\
//         "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \\
//         "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \\
//         "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \\
//         > ${sample}.snpSift.table.txt

//     snpEff ${index_base} \\
//         -config $config \\
//         -dataDir $db \\
//         ${highfreq_vcf[0]} \\
//         -csvStats ${prefix}.snpEff.csv \\
//         | bgzip -c > ${prefix}.snpEff.vcf.gz
//     tabix -p vcf -f ${prefix}.snpEff.vcf.gz
//     mv snpEff_summary.html ${prefix}.snpEff.summary.html

//     SnpSift extractFields -s "," \\
//         -e "." \\
//         ${prefix}.snpEff.vcf.gz \\
//         CHROM POS REF ALT \\
//         "ANN[*].GENE" "ANN[*].GENEID" \\
//         "ANN[*].IMPACT" "ANN[*].EFFECT" \\
//         "ANN[*].FEATURE" "ANN[*].FEATUREID" \\
//         "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \\
//         "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \\
//         "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \\
//         "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \\
//         "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \\
//         > ${prefix}.snpSift.table.txt
//     	"""
// }

// /*
//  * STEP 5.7.1.3: VarScan 2 consensus sequence report with QUAST
//  */
// process VARSCAN2_QUAST {
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/varscan2/quast", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (!filename.endsWith(".tsv")) filename
//                 }

//     when:
//     !params.skip_variants && 'varscan2' in callers && !params.skip_variants_quast

//     input:
//     path consensus from ch_varscan2_consensus.collect{ it[2] }
//     path fasta from ch_fasta
//     path gff from ch_gff

//     output:
//     path "AF${params.max_allele_freq}"
//     path "report.tsv" into ch_varscan2_quast_mqc

//     script:
//     features = params.gff ? "--features $gff" : ""
//     """
//     quast.py \\
//         --output-dir AF${params.max_allele_freq} \\
//         -r $fasta \\
//         $features \\
//         --threads $task.cpus \\
//         ${consensus.join(' ')}
//     ln -s AF${params.max_allele_freq}/report.tsv
//     """
// }

// ////////////////////////////////////////////////////
// /* --                IVAR                      -- */
// ////////////////////////////////////////////////////

// /*
//  * STEP 5.7.2: Variant calling with iVar
//  */
// process IVAR_VARIANTS {
//     tag "$sample"
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/ivar", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".bcftools_stats.txt")) "bcftools_stats/$filename"
//                       else if (filename.endsWith(".log")) "log/$filename"
//                       else if (filename.endsWith("_mqc.tsv")) null
//                       else filename
//                 }

//     when:
//     !params.skip_variants && 'ivar' in callers

//     input:
//     tuple val(sample), val(single_end), path(mpileup) from ch_mpileup_ivar_variants
//     path header from ch_ivar_variants_header_mqc
//     path fasta from ch_fasta
//     path gff from ch_gff

//     output:
//     tuple val(sample), val(single_end), path("${prefix}.vcf.gz*") into ch_ivar_highfreq_snpeff,
//                                                                        ch_ivar_highfreq_intersect
//     tuple val(sample), val(single_end), path("${sample}.vcf.gz*") into ch_ivar_lowfreq_snpeff
//     path "${prefix}.bcftools_stats.txt" into ch_ivar_bcftools_highfreq_mqc
//     path "${sample}.variant.counts_mqc.tsv" into ch_ivar_count_mqc
//     path "${sample}.bcftools_stats.txt"
//     path "${sample}.tsv"
//     path "*.log"

//     script:
//     features = params.gff ? "-g $gff" : ""
//     prefix = "${sample}.AF${params.max_allele_freq}"
//     """
//     cat $mpileup | ivar variants -q $params.min_base_qual -t $params.min_allele_freq -m $params.min_coverage -r $fasta $features -p $sample

//     ivar_variants_to_vcf.py ${sample}.tsv ${sample}.vcf > ${sample}.variant.counts.log
//     bgzip -c ${sample}.vcf > ${sample}.vcf.gz
//     tabix -p vcf -f ${sample}.vcf.gz
//     bcftools stats ${sample}.vcf.gz > ${sample}.bcftools_stats.txt
//     cat $header ${sample}.variant.counts.log > ${sample}.variant.counts_mqc.tsv

//     ivar_variants_to_vcf.py ${sample}.tsv ${prefix}.vcf --pass_only --allele_freq_thresh $params.max_allele_freq > ${prefix}.variant.counts.log
//     bgzip -c ${prefix}.vcf > ${prefix}.vcf.gz
//     tabix -p vcf -f ${prefix}.vcf.gz
//     bcftools stats ${prefix}.vcf.gz > ${prefix}.bcftools_stats.txt
//     """
// }

// /*
//  * STEP 5.7.2.1: Generate consensus sequence with iVar
//  */
// process IVAR_CONSENSUS {
//     tag "$sample"
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/ivar/consensus", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".tsv")) "base_qc/$filename"
//                       else if (filename.endsWith(".pdf")) "base_qc/$filename"
//                       else filename
//                 }

//     when:
//     !params.skip_variants && 'ivar' in callers

//     input:
//     tuple val(sample), val(single_end), path(mpileup) from ch_mpileup_ivar_consensus
//     path fasta from ch_fasta

//     output:
//     tuple val(sample), val(single_end), path("*.fa") into ch_ivar_consensus
//     path "*.{txt,tsv,pdf}"

//     script:
//     prefix = "${sample}.AF${params.max_allele_freq}"
//     """
//     cat $mpileup | ivar consensus -q $params.min_base_qual -t $params.max_allele_freq -m $params.min_coverage -n N -p ${prefix}.consensus
//     header=\$(head -n1 ${prefix}.consensus.fa | sed 's/>//g')
//     sed -i "s/\${header}/${sample}/g" ${prefix}.consensus.fa

//     plot_base_density.r --fasta_files ${prefix}.consensus.fa --prefixes $prefix --output_dir ./
//     """
// }

// /*
//  * STEP 5.7.2.2: iVar variant calling annotation with SnpEff and SnpSift
//  */
// process IVAR_SNPEFF {
//     tag "$sample"
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/ivar/snpeff", mode: params.publish_dir_mode

//     when:
//     !params.skip_variants && 'ivar' in callers && params.gff && !params.skip_snpeff

//     input:
//     tuple val(sample), val(single_end), path(highfreq_vcf), path(lowfreq_vcf) from ch_ivar_highfreq_snpeff.join(ch_ivar_lowfreq_snpeff, by: [0,1])
//     tuple file(db), file(config) from ch_snpeff_db_ivar

//     output:
//     path "${prefix}.snpEff.csv" into ch_ivar_snpeff_highfreq_mqc
//     path "${sample}.snpEff.csv"
//     path "*.vcf.gz*"
//     path "*.{txt,html}"

//     script:
//     prefix = "${sample}.AF${params.max_allele_freq}"
//     """
//     snpEff ${index_base} \\
//         -config $config \\
//         -dataDir $db \\
//         ${lowfreq_vcf[0]} \\
//         -csvStats ${sample}.snpEff.csv \\
//         | bgzip -c > ${sample}.snpEff.vcf.gz
//     tabix -p vcf -f ${sample}.snpEff.vcf.gz
//     mv snpEff_summary.html ${sample}.snpEff.summary.html

//     SnpSift extractFields -s "," \\
//         -e "." \\
//         ${sample}.snpEff.vcf.gz \\
//         CHROM POS REF ALT \\
//         "ANN[*].GENE" "ANN[*].GENEID" \\
//         "ANN[*].IMPACT" "ANN[*].EFFECT" \\
//         "ANN[*].FEATURE" "ANN[*].FEATUREID" \\
//         "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \\
//         "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \\
//         "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \\
//         "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \\
//         "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \\
//         > ${sample}.snpSift.table.txt

//     snpEff ${index_base} \\
//         -config $config \\
//         -dataDir $db \\
//         ${highfreq_vcf[0]} \\
//         -csvStats ${prefix}.snpEff.csv \\
//         | bgzip -c > ${prefix}.snpEff.vcf.gz
//     tabix -p vcf -f ${prefix}.snpEff.vcf.gz
//     mv snpEff_summary.html ${prefix}.snpEff.summary.html

//     SnpSift extractFields -s "," \\
//         -e "." \\
//         ${prefix}.snpEff.vcf.gz \\
//         CHROM POS REF ALT \\
//         "ANN[*].GENE" "ANN[*].GENEID" \\
//         "ANN[*].IMPACT" "ANN[*].EFFECT" \\
//         "ANN[*].FEATURE" "ANN[*].FEATUREID" \\
//         "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \\
//         "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \\
//         "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \\
//         "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \\
//         "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \\
//         > ${prefix}.snpSift.table.txt
//    	"""
// }

// /*
//  * STEP 5.7.2.3: iVar consensus sequence report with QUAST
//  */
// process IVAR_QUAST {
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/ivar/quast", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (!filename.endsWith(".tsv")) filename
//                 }

//     when:
//     !params.skip_variants && 'ivar' in callers && !params.skip_variants_quast

//     input:
//     path consensus from ch_ivar_consensus.collect{ it[2] }
//     path fasta from ch_fasta
//     path gff from ch_gff

//     output:
//     path "AF${params.max_allele_freq}"
//     path "report.tsv" into ch_ivar_quast_mqc

//     script:
//     features = params.gff ? "--features $gff" : ""
//     """
//     quast.py \\
//         --output-dir AF${params.max_allele_freq} \\
//         -r $fasta \\
//         $features \\
//         --threads $task.cpus \\
//         ${consensus.join(' ')}
//     ln -s AF${params.max_allele_freq}/report.tsv
//     """
// }

// ////////////////////////////////////////////////////
// /* --              BCFTOOLS                    -- */
// ////////////////////////////////////////////////////

// /*
//  * STEP 5.7.3: Variant calling with BCFTools
//  */
// process BCFTOOLS_VARIANTS {
//     tag "$sample"
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/bcftools", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".txt")) "bcftools_stats/$filename"
//                       else filename
//                 }

//     when:
//     !params.skip_variants && 'bcftools' in callers

//     input:
//     tuple val(sample), val(single_end), path(bam) from ch_markdup_bam_bcftools
//     path fasta from ch_fasta

//     output:
//     tuple val(sample), val(single_end), path("*.vcf.gz*") into ch_bcftools_variants_consensus,
//                                                                ch_bcftools_variants_snpeff,
//                                                                ch_bcftools_variants_intersect
//     path "*.bcftools_stats.txt" into ch_bcftools_variants_mqc

//     script:
//     """
//     echo "$sample" > sample_name.list
//     bcftools mpileup \\
//         --count-orphans \\
//         --no-BAQ \\
//         --max-depth $params.mpileup_depth \\
//         --fasta-ref $fasta \\
//         --min-BQ $params.min_base_qual \\
//         --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
//         ${bam[0]} \\
//         | bcftools call --output-type v --ploidy 1 --keep-alts --keep-masked-ref --multiallelic-caller --variants-only \\
//         | bcftools reheader --samples sample_name.list \\
//         | bcftools view --output-file ${sample}.vcf.gz --output-type z --include 'INFO/DP>=$params.min_coverage'
//     tabix -p vcf -f ${sample}.vcf.gz
//     bcftools stats ${sample}.vcf.gz > ${sample}.bcftools_stats.txt
//     """
// }

// /*
//  * STEP 5.7.3.1: Genome consensus generation with BCFtools and masked with BEDTools
//  */
// process BCFTOOLS_CONSENSUS {
//     tag "$sample"
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/bcftools/consensus", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".tsv")) "base_qc/$filename"
//                       else if (filename.endsWith(".pdf")) "base_qc/$filename"
//                       else filename
//                 }

//     when:
//     !params.skip_variants && 'bcftools' in callers

//     input:
//     tuple val(sample), val(single_end), path(bam), path(vcf) from ch_markdup_bam_bcftools_consensus.join(ch_bcftools_variants_consensus, by: [0,1])
//     path fasta from ch_fasta

//     output:
//     tuple val(sample), val(single_end), path("*consensus.masked.fa") into ch_bcftools_consensus_masked
//     path "*.{consensus.fa,tsv,pdf}"

//     script:
//     """
//     bedtools genomecov \\
//         -bga \\
//         -ibam ${bam[0]} \\
//         -g $fasta \\
//         | awk '\$4 < $params.min_coverage' > ${sample}.lowcov.bed

//     parse_mask_bed.py ${vcf[0]} ${sample}.lowcov.bed ${sample}.lowcov.fix.bed

//     bedtools merge -i ${sample}.lowcov.fix.bed > ${sample}.mask.bed

//     bedtools maskfasta \\
//         -fi $fasta \\
//         -bed ${sample}.mask.bed \\
//         -fo ${index_base}.ref.masked.fa
        
//     cat ${index_base}.ref.masked.fa | bcftools consensus ${vcf[0]} > ${sample}.consensus.masked.fa

//     sed -i 's/${index_base}/${sample}/g' ${sample}.consensus.masked.fa
//     header=\$(head -n1 ${sample}.consensus.masked.fa | sed 's/>//g')
//     sed -i "s/\${header}/${sample}/g" ${sample}.consensus.masked.fa

//     plot_base_density.r --fasta_files ${sample}.consensus.masked.fa --prefixes $sample --output_dir ./
//     """
// }

// /*
//  * STEP 5.7.3.2: BCFTools variant calling annotation with SnpEff and SnpSift
//  */
// process BCFTOOLS_SNPEFF {
//     tag "$sample"
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/bcftools/snpeff", mode: params.publish_dir_mode

//     when:
//     !params.skip_variants && 'bcftools' in callers && params.gff && !params.skip_snpeff

//     input:
//     tuple val(sample), val(single_end), path(vcf) from ch_bcftools_variants_snpeff
//     tuple file(db), file(config) from ch_snpeff_db_bcftools

//     output:
//     path "*.snpEff.csv" into ch_bcftools_snpeff_mqc
//     path "*.vcf.gz*"
//     path "*.{txt,html}"

//     script:
//     """
//     snpEff ${index_base} \\
//         -config $config \\
//         -dataDir $db \\
//         ${vcf[0]} \\
//         -csvStats ${sample}.snpEff.csv \\
//         | bgzip -c > ${sample}.snpEff.vcf.gz
//     tabix -p vcf -f ${sample}.snpEff.vcf.gz
//     mv snpEff_summary.html ${sample}.snpEff.summary.html

//     SnpSift extractFields -s "," \\
//         -e "." \\
//         ${sample}.snpEff.vcf.gz \\
//         CHROM POS REF ALT \\
//         "ANN[*].GENE" "ANN[*].GENEID" \\
//         "ANN[*].IMPACT" "ANN[*].EFFECT" \\
//         "ANN[*].FEATURE" "ANN[*].FEATUREID" \\
//         "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \\
//         "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \\
//         "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \\
//         "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \\
//         "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \\
//         > ${sample}.snpSift.table.txt
//     	"""
// }

// /*
//  * STEP 5.7.3.3: BCFTools consensus sequence report with QUAST
//  */
// process BCFTOOLS_QUAST {
//     label 'process_medium'
//     publishDir "${params.outdir}/variants/bcftools", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (!filename.endsWith(".tsv")) filename
//                 }

//     when:
//     !params.skip_variants && 'bcftools' in callers && !params.skip_variants_quast

//     input:
//     path consensus from ch_bcftools_consensus_masked.collect{ it[2] }
//     path fasta from ch_fasta
//     path gff from ch_gff

//     output:
//     path "quast"
//     path "report.tsv" into ch_bcftools_quast_mqc

//     script:
//     features = params.gff ? "--features $gff" : ""
//     """
//     quast.py \\
//         --output-dir quast \\
//         -r $fasta \\
//         $features \\
//         --threads $task.cpus \\
//         ${consensus.join(' ')}
//     ln -s quast/report.tsv
//     """
// }

// ////////////////////////////////////////////////////
// /* --            INTERSECT VARIANTS            -- */
// ////////////////////////////////////////////////////

// /*
//  * STEP 5.8: Intersect variants with BCFTools
//  */
// if (!params.skip_variants && callers.size() > 2) {

//     ch_varscan2_highfreq_intersect
//         .join(ch_ivar_highfreq_intersect, by: [0,1])
//         .join(ch_bcftools_variants_intersect, by: [0,1])
//         .set { ch_varscan2_highfreq_intersect }

//     process BCFTOOLS_ISEC {
//         tag "$sample"
//         label 'process_medium'
//         label 'error_ignore'
//         publishDir "${params.outdir}/variants/intersect", mode: params.publish_dir_mode

//         input:
//         tuple val(sample), val(single_end), path('varscan2/*'), path('ivar/*'), path('bcftools/*') from ch_varscan2_highfreq_intersect

//         output:
//         path "$sample"

//         script:
//         """
//         bcftools isec  \\
//             --nfiles +2 \\
//             --output-type z \\
//             -p $sample \\
//             */*.vcf.gz
//         """
//     }
// }

// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                    DENOVO ASSEMBLY PROCESSES                        -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////

// /*
//  * PREPROCESSING: Build Blast database for viral genome
//  */
// process MAKE_BLAST_DB {
//     tag "$fasta"
//     label 'process_medium'
//     if (params.save_reference) {
//         publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
//     }

//     when:
//     !params.skip_assembly && !params.skip_blast

//     input:
//     path fasta from ch_fasta

//     output:
//     path "BlastDB" into ch_blast_db

//     script:
//     """
//     makeblastdb \\
//         -in $fasta \\
//         -parse_seqids \\
//         -dbtype nucl
//     mkdir BlastDB && mv ${fasta}* BlastDB
//     """
// }

// /*
//  * PREPROCESSING: Build Kraken2 database for host genome
//  */
// if (!isOffline()) {
//     if (!params.skip_kraken2 && !params.kraken2_db) {
//         if (!params.kraken2_db_name) { exit 1, "Please specify a valid name to build Kraken2 database for host e.g. 'human'!" }

//         process KRAKEN2_BUILD {
//             tag "$db"
//             label 'process_high'
//             if (params.save_reference) {
//                 publishDir "${params.outdir}/genome", mode: params.publish_dir_mode
//             }

//             when:
//             !params.skip_assembly

//             output:
//             path "$db" into ch_kraken2_db

//             script:
//             db = "kraken2_${params.kraken2_db_name}"
//             ftp = params.kraken2_use_ftp ? "--use-ftp" : ""
//             """
//             kraken2-build --db $db --threads $task.cpus $ftp --download-taxonomy
//             kraken2-build --db $db --threads $task.cpus $ftp --download-library $params.kraken2_db_name
//             kraken2-build --db $db --threads $task.cpus $ftp --build
//             """
//         }
//     }
// } else {
//     exit 1, "NXF_OFFLINE=true or -offline has been set so cannot download files required to build Kraken2 database!"
// }

// /*
//  * STEP 6.1: Amplicon trimming with Cutadapt
//  */
// if (params.protocol == 'amplicon' && !params.skip_assembly && !params.skip_amplicon_trimming) {
//     process CUTADAPT {
//         tag "$sample"
//         label 'process_medium'
//         publishDir "${params.outdir}/assembly/cutadapt", mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                           if (filename.endsWith(".html")) "fastqc/$filename"
//                           else if (filename.endsWith(".zip")) "fastqc/zips/$filename"
//                           else if (filename.endsWith(".log")) "log/$filename"
//                           else params.save_trimmed ? filename : null
//                     }

//         input:
//         tuple val(sample), val(single_end), path(reads) from ch_fastp_cutadapt
//         path amplicons from ch_amplicon_fasta

//         output:
//         tuple val(sample), val(single_end), path("*.ptrim.fastq.gz") into ch_cutadapt_kraken2
//         path "*.{zip,html}" into ch_cutadapt_fastqc_mqc
//         path "*.log" into ch_cutadapt_mqc

//         script:
//         adapters = single_end ? "-a file:primers.fasta" : "-a file:primers.fasta -A file:primers.fasta"
//         out_reads = single_end ? "-o ${sample}.ptrim.fastq.gz" : "-o ${sample}_1.ptrim.fastq.gz -p ${sample}_2.ptrim.fastq.gz"
//         """
//         sed -r '/^[ACTGactg]+\$/ s/\$/X/g' $amplicons > primers.fasta

//         cutadapt \\
//             --cores $task.cpus \\
//             --overlap 5 \\
//             --minimum-length 30 \\
//             --error-rate 0.1 \\
//             $adapters \\
//             $out_reads \\
//             $reads \\
//             > ${sample}.cutadapt.log

//         fastqc --quiet --threads $task.cpus *.ptrim.fastq.gz
//         """
//     }
//     ch_fastp_kraken2 = ch_cutadapt_kraken2

// } else {
//     ch_cutadapt_mqc = Channel.empty()
//     ch_cutadapt_fastqc_mqc = Channel.empty()
// }

// /*
//  * STEP 6.2: Filter reads with Kraken2
//  */
// if (!params.skip_kraken2 && !params.skip_assembly) {
//     process KRAKEN2 {
//         tag "$sample"
//         label 'process_high'
//         publishDir "${params.outdir}/assembly/kraken2", mode: params.publish_dir_mode,
//             saveAs: { filename ->
//                           if (filename.endsWith(".txt")) filename
//                           else params.save_kraken2_fastq ? filename : null
//                     }

//         input:
//         tuple val(sample), val(single_end), path(reads) from ch_fastp_kraken2
//         path db from ch_kraken2_db

//         output:
//         tuple val(sample), val(single_end), path("*.viral*") into ch_kraken2_spades,
//                                                                   ch_kraken2_metaspades,
//                                                                   ch_kraken2_unicycler,
//                                                                   ch_kraken2_minia
//         path "*.report.txt" into ch_kraken2_report_mqc
//         path "*.host*"


//         script:
//         pe = single_end ? "" : "--paired"
//         classified = single_end ? "${sample}.host.fastq" : "${sample}.host#.fastq"
//         unclassified = single_end ? "${sample}.viral.fastq" : "${sample}.viral#.fastq"
//         """
//         kraken2 \\
//             --db $db \\
//             --threads $task.cpus \\
//             --unclassified-out $unclassified \\
//             --classified-out $classified \\
//             --report ${sample}.kraken2.report.txt \\
//             --report-zero-counts \\
//             $pe \\
//             --gzip-compressed \\
//             $reads
//         pigz -p $task.cpus *.fastq
//         """
//     }
// } else {
//     ch_fastp_kraken2
//         .into { ch_kraken2_spades
//                 ch_kraken2_metaspades
//                 ch_kraken2_unicycler
//                 ch_kraken2_minia }
//     ch_kraken2_report_mqc = Channel.empty()
// }

// ////////////////////////////////////////////////////
// /* --                SPADES                    -- */
// ////////////////////////////////////////////////////

// /*
//  * STEP 6.3: De novo assembly with SPAdes
//  */
// process SPADES {
//     tag "$sample"
//     label 'process_high'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/spades", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".png")) "bandage/$filename"
//                       else if (filename.endsWith(".svg")) "bandage/$filename"
//                       else filename
//                 }

//     when:
//     !params.skip_assembly && 'spades' in assemblers

//     input:
//     tuple val(sample), val(single_end), path(reads) from ch_kraken2_spades

//     output:
//     tuple val(sample), val(single_end), path("*scaffolds.fa") into ch_spades_blast,
//                                                                    ch_spades_abacas,
//                                                                    ch_spades_plasmidid,
//                                                                    ch_spades_quast,
//                                                                    ch_spades_vg
//     path "*assembly.{gfa,png,svg}"


//     script:
//     input_reads = single_end ? "-s $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
//     """
//     spades.py \\
//         --threads $task.cpus \\
//         $input_reads \\
//         -o ./
//     mv scaffolds.fasta ${sample}.scaffolds.fa
//     mv assembly_graph_with_scaffolds.gfa ${sample}.assembly.gfa

//     if [ -s ${sample}.assembly.gfa ]
//     then
//         Bandage image ${sample}.assembly.gfa ${sample}.assembly.png --height 1000
//         Bandage image ${sample}.assembly.gfa ${sample}.assembly.svg --height 1000
//     fi
//     """
// }

// /*
//  * STEP 6.3.1: Run Blast on SPAdes de novo assembly
//  */
// process SPADES_BLAST {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/spades/blast", mode: params.publish_dir_mode

//     when:
//     !params.skip_assembly && 'spades' in assemblers && !params.skip_blast

//     input:
//     tuple val(sample), val(single_end), path(scaffold) from ch_spades_blast
//     path db from ch_blast_db
//     path header from ch_blast_outfmt6_header

//     output:
//     path "*.blast*"

//     script:
//     """
//     blastn \\
//         -num_threads $task.cpus \\
//         -db $db/$fasta_base \\
//         -query $scaffold \\
//         -outfmt \'6 stitle std slen qlen qcovs\' \\
//         -out ${sample}.blast.txt

//     awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"}{print \$0,\$5/\$15,\$5/\$14}' ${sample}.blast.txt | awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"} \$15 > 200 && \$17 > 0.7 && \$1 !~ /phage/ {print \$0}' > ${sample}.blast.filt.txt
//     cat $header ${sample}.blast.filt.txt > ${sample}.blast.filt.header.txt
//     """
// }

// /*
//  * STEP 6.3.2: Run ABACAS on SPAdes de novo assembly
//  */
// process SPADES_ABACAS {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/spades/abacas", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.indexOf("nucmer") > 0) "nucmer/$filename"
//                       else filename
//                 }

//     when:
//     !params.skip_assembly && 'spades' in assemblers && !params.skip_abacas

//     input:
//     tuple val(sample), val(single_end), path(scaffold) from ch_spades_abacas
//     path fasta from ch_fasta

//     output:
//     path "*.abacas*"

//     script:
//     """
//     abacas.pl -r $fasta -q $scaffold -m -p nucmer -o ${sample}.abacas
//     mv nucmer.delta ${sample}.abacas.nucmer.delta
//     mv nucmer.filtered.delta ${sample}.abacas.nucmer.filtered.delta
//     mv nucmer.tiling ${sample}.abacas.nucmer.tiling
//     mv unused_contigs.out ${sample}.abacas.unused.contigs.out
//     """
// }

// /*
//  * STEP 6.3.3: Run PlasmidID on SPAdes de novo assembly
//  */
// process SPADES_PLASMIDID {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/spades/plasmidid", mode: params.publish_dir_mode

//     when:
//     !params.skip_assembly && 'spades' in assemblers && !params.skip_plasmidid

//     input:
//     tuple val(sample), val(single_end), path(scaffold) from ch_spades_plasmidid.filter { it.size() > 0 }
//     path fasta from ch_fasta

//     output:
//     path "$sample"

//     script:
//     """
//     plasmidID -d $fasta -s $sample -c $scaffold --only-reconstruct -C 47 -S 47 -i 60 --no-trim -o .
//     mv NO_GROUP/$sample ./$sample
//     """
// }

// /*
//  * STEP 6.3.4: Run Quast on SPAdes de novo assembly
//  */
// process SPADES_QUAST {
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/spades", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (!filename.endsWith(".tsv")) filename
//                 }

//     when:
//     !params.skip_assembly && 'spades' in assemblers && !params.skip_assembly_quast

//     input:
//     path scaffolds from ch_spades_quast.collect{ it[2] }
//     path fasta from ch_fasta
//     path gff from ch_gff

//     output:
//     path "quast"
//     path "report.tsv" into ch_quast_spades_mqc

//     script:
//     features = params.gff ? "--features $gff" : ""
//     """
//     quast.py \\
//         --output-dir quast \\
//         -r $fasta \\
//         $features \\
//         --threads $task.cpus \\
//         ${scaffolds.join(' ')}
//     ln -s quast/report.tsv
//     """
// }

// /*
//  * STEP 6.3.5: Overlap scaffolds with Minimap2, induce and polish assembly, and call variants with seqwish and vg
//  */
// process SPADES_VG {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/spades/variants", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".txt")) "bcftools_stats/$filename"
//                       else if (filename.endsWith(".png")) "bandage/$filename"
//                       else if (filename.endsWith(".svg")) "bandage/$filename"
//                       else filename
//                 }

//     when:
//     !params.skip_assembly && 'spades' in assemblers && !params.skip_vg

//     input:
//     tuple val(sample), val(single_end), path(scaffolds) from ch_spades_vg
//     path fasta from ch_fasta

//     output:
//     tuple val(sample), val(single_end), path("${sample}.vcf.gz*") into ch_spades_vg_vcf
//     path "*.bcftools_stats.txt" into ch_spades_vg_bcftools_mqc
//     path "*.{gfa,png,svg}"

//     script:
//     """
//     minimap2 -c -t $task.cpus -x asm20 $fasta $scaffolds > ${sample}.paf

//     cat $scaffolds $fasta > ${sample}.withRef.fasta
//     seqwish --paf-alns ${sample}.paf --seqs ${sample}.withRef.fasta --gfa ${sample}.gfa --threads $task.cpus

//     vg view -Fv ${sample}.gfa --threads $task.cpus > ${sample}.vg
//     vg convert -x ${sample}.vg > ${sample}.xg

//     samtools faidx $fasta
//     vg snarls ${sample}.xg > ${sample}.snarls
//     for chrom in `cat ${fasta}.fai | cut -f1`
//     do
//         vg deconstruct -p \$chrom ${sample}.xg -r ${sample}.snarls --threads $task.cpus \\
//             | bcftools sort -O v -T ./ \\
//             | bgzip -c > ${sample}.\$chrom.vcf.gz
//     done
//     bcftools concat --output-type z --output ${sample}.vcf.gz *.vcf.gz
//     tabix -p vcf -f ${sample}.vcf.gz
//     bcftools stats ${sample}.vcf.gz > ${sample}.bcftools_stats.txt

//     if [ -s ${sample}.gfa ]
//     then
//         Bandage image ${sample}.gfa ${sample}.png --height 1000
//         Bandage image ${sample}.gfa ${sample}.svg --height 1000
//     fi
//     """
// }

// /*
//  * STEP 6.3.6: Variant annotation with SnpEff and SnpSift
//  */
// process SPADES_SNPEFF {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/spades/variants/snpeff", mode: params.publish_dir_mode

//     when:
//     !params.skip_assembly && 'spades' in assemblers && !params.skip_vg && params.gff && !params.skip_snpeff

//     input:
//     tuple val(sample), val(single_end), path(vcf) from ch_spades_vg_vcf
//     tuple file(db), file(config) from ch_snpeff_db_spades

//     output:
//     path "*.snpEff.csv" into ch_spades_snpeff_mqc
//     path "*.vcf.gz*"
//     path "*.{txt,html}"

//     script:
//     """
//     snpEff ${index_base} \\
//         -config $config \\
//         -dataDir $db \\
//         ${vcf[0]} \\
//         -csvStats ${sample}.snpEff.csv \\
//         | bgzip -c > ${sample}.snpEff.vcf.gz
//     tabix -p vcf -f ${sample}.snpEff.vcf.gz
//     mv snpEff_summary.html ${sample}.snpEff.summary.html

//     SnpSift extractFields -s "," \\
//         -e "." \\
//         ${sample}.snpEff.vcf.gz \\
//         CHROM POS REF ALT \\
//         "ANN[*].GENE" "ANN[*].GENEID" \\
//         "ANN[*].IMPACT" "ANN[*].EFFECT" \\
//         "ANN[*].FEATURE" "ANN[*].FEATUREID" \\
//         "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \\
//         "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \\
//         "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \\
//         "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \\
//         "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \\
//         > ${sample}.snpSift.table.txt
//     	"""
// }

// ////////////////////////////////////////////////////
// /* --               METASPADES                 -- */
// ////////////////////////////////////////////////////

// /*
//  * STEP 6.3: De novo assembly with MetaSPAdes
//  */
// process METASPADES {
//     tag "$sample"
//     label 'process_high'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/metaspades", mode: params.publish_dir_mode,
//     saveAs: { filename ->
//                   if (filename.endsWith(".png")) "bandage/$filename"
//                   else if (filename.endsWith(".svg")) "bandage/$filename"
//                   else filename
//             }

//     when:
//     !params.skip_assembly && 'metaspades' in assemblers && !single_end

//     input:
//     tuple val(sample), val(single_end), path(reads) from ch_kraken2_metaspades

//     output:
//     tuple val(sample), val(single_end), path("*scaffolds.fa") into ch_metaspades_blast,
//                                                                    ch_metaspades_abacas,
//                                                                    ch_metaspades_plasmidid,
//                                                                    ch_metaspades_quast,
//                                                                    ch_metaspades_vg
//     path "*assembly.{gfa,png,svg}"


//     script:
//     """
//     spades.py \\
//         --meta \\
//         --threads $task.cpus \\
//         -1 ${reads[0]} \\
//         -2 ${reads[1]} \\
//         -o ./
//     mv scaffolds.fasta ${sample}.scaffolds.fa
//     mv assembly_graph_with_scaffolds.gfa ${sample}.assembly.gfa

//     if [ -s ${sample}.assembly.gfa ]
//     then
//         Bandage image ${sample}.assembly.gfa ${sample}.assembly.png --height 1000
//         Bandage image ${sample}.assembly.gfa ${sample}.assembly.svg --height 1000
//     fi
//     """
// }

// /*
//  * STEP 6.3.1: Run Blast on MetaSPAdes de novo assembly
//  */
// process METASPADES_BLAST {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/metaspades/blast", mode: params.publish_dir_mode

//     when:
//     !params.skip_assembly && 'metaspades' in assemblers && !single_end && !params.skip_blast

//     input:
//     tuple val(sample), val(single_end), path(scaffold) from ch_metaspades_blast
//     path db from ch_blast_db
//     path header from ch_blast_outfmt6_header

//     output:
//     path "*.blast*"

//     script:
//     """
//     blastn \\
//         -num_threads $task.cpus \\
//         -db $db/$fasta_base \\
//         -query $scaffold \\
//         -outfmt \'6 stitle std slen qlen qcovs\' \\
//         -out ${sample}.blast.txt

//     awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"}{print \$0,\$5/\$15,\$5/\$14}' ${sample}.blast.txt | awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"} \$15 > 200 && \$17 > 0.7 && \$1 !~ /phage/ {print \$0}' > ${sample}.blast.filt.txt
//     cat $header ${sample}.blast.filt.txt > ${sample}.blast.filt.header.txt
//     """
// }

// /*
//  * STEP 6.3.2: Run ABACAS on MetaSPAdes de novo assembly
//  */
// process METASPADES_ABACAS {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/metaspades/abacas", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.indexOf("nucmer") > 0) "nucmer/$filename"
//                       else filename
//                 }

//     when:
//     !params.skip_assembly && 'metaspades' in assemblers && !single_end && !params.skip_abacas

//     input:
//     tuple val(sample), val(single_end), path(scaffold) from ch_metaspades_abacas
//     path fasta from ch_fasta

//     output:
//     path "*.abacas*"

//     script:
//     """
//     abacas.pl -r $fasta -q $scaffold -m -p nucmer -o ${sample}.abacas
//     mv nucmer.delta ${sample}.abacas.nucmer.delta
//     mv nucmer.filtered.delta ${sample}.abacas.nucmer.filtered.delta
//     mv nucmer.tiling ${sample}.abacas.nucmer.tiling
//     mv unused_contigs.out ${sample}.abacas.unused.contigs.out
//     """
// }

// /*
//  * STEP 6.3.3: Run PlasmidID on MetaSPAdes de novo assembly
//  */
// process METASPADES_PLASMIDID {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/metaspades/plasmidid", mode: params.publish_dir_mode

//     when:
//     !params.skip_assembly && 'metaspades' in assemblers && !single_end && !params.skip_plasmidid

//     input:
//     tuple val(sample), val(single_end), path(scaffold) from ch_metaspades_plasmidid.filter { it.size() > 0 }
//     path fasta from ch_fasta

//     output:
//     path "$sample"

//     script:
//     """
//     plasmidID -d $fasta -s $sample -c $scaffold --only-reconstruct -C 47 -S 47 -i 60 --no-trim -o .
//     mv NO_GROUP/$sample ./$sample
//     """
// }

// /*
//  * STEP 6.3.4: Run Quast on MetaSPAdes de novo assembly
//  */
// process METASPADES_QUAST {
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/metaspades", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (!filename.endsWith(".tsv")) filename
//                 }

//     when:
//     !params.skip_assembly && 'metaspades' in assemblers && !single_end && !params.skip_assembly_quast

//     input:
//     path scaffolds from ch_metaspades_quast.collect{ it[2] }
//     path fasta from ch_fasta
//     path gff from ch_gff

//     output:
//     path "quast"
//     path "report.tsv" into ch_quast_metaspades_mqc

//     script:
//     features = params.gff ? "--features $gff" : ""
//     """
//     quast.py \\
//         --output-dir quast \\
//         -r $fasta \\
//         $features \\
//         --threads $task.cpus \\
//         ${scaffolds.join(' ')}
//     ln -s quast/report.tsv
//     """
// }

// /*
//  * STEP 6.3.5: Overlap scaffolds with Minimap2, induce and polish assembly, and call variants with seqwish and vg
//  */
// process METASPADES_VG {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/metaspades/variants", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".txt")) "bcftools_stats/$filename"
//                       else if (filename.endsWith(".png")) "bandage/$filename"
//                       else if (filename.endsWith(".svg")) "bandage/$filename"
//                       else filename
//                 }

//     when:
//     !params.skip_assembly && 'metaspades' in assemblers && !single_end && !params.skip_vg

//     input:
//     tuple val(sample), val(single_end), path(scaffolds) from ch_metaspades_vg
//     path fasta from ch_fasta

//     output:
//     tuple val(sample), val(single_end), path("${sample}.vcf.gz*") into ch_metaspades_vg_vcf
//     path "*.bcftools_stats.txt" into ch_metaspades_vg_bcftools_mqc
//     path "*.{gfa,png,svg}"

//     script:
//     """
//     minimap2 -c -t $task.cpus -x asm20 $fasta $scaffolds > ${sample}.paf

//     cat $scaffolds $fasta > ${sample}.withRef.fasta
//     seqwish --paf-alns ${sample}.paf --seqs ${sample}.withRef.fasta --gfa ${sample}.gfa --threads $task.cpus

//     vg view -Fv ${sample}.gfa --threads $task.cpus > ${sample}.vg
//     vg convert -x ${sample}.vg > ${sample}.xg

//     samtools faidx $fasta
//     vg snarls ${sample}.xg > ${sample}.snarls
//     for chrom in `cat ${fasta}.fai | cut -f1`
//     do
//         vg deconstruct -p \$chrom ${sample}.xg -r ${sample}.snarls --threads $task.cpus \\
//             | bcftools sort -O v -T ./ \\
//             | bgzip -c > ${sample}.\$chrom.vcf.gz
//     done
//     bcftools concat --output-type z --output ${sample}.vcf.gz *.vcf.gz
//     tabix -p vcf -f ${sample}.vcf.gz
//     bcftools stats ${sample}.vcf.gz > ${sample}.bcftools_stats.txt

//     if [ -s ${sample}.gfa ]
//     then
//         Bandage image ${sample}.gfa ${sample}.png --height 1000
//         Bandage image ${sample}.gfa ${sample}.svg --height 1000
//     fi
//     """
// }

// /*
//  * STEP 6.3.6: Variant annotation with SnpEff and SnpSift
//  */
// process METASPADES_SNPEFF {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/metaspades/variants/snpeff", mode: params.publish_dir_mode

//     when:
//     !params.skip_assembly && 'metaspades' in assemblers && !single_end && !params.skip_vg && params.gff && !params.skip_snpeff

//     input:
//     tuple val(sample), val(single_end), path(vcf) from ch_metaspades_vg_vcf
//     tuple file(db), file(config) from ch_snpeff_db_metaspades

//     output:
//     path "*.snpEff.csv" into ch_metaspades_snpeff_mqc
//     path "*.vcf.gz*"
//     path "*.{txt,html}"

//     script:
//     """
//     snpEff ${index_base} \\
//         -config $config \\
//         -dataDir $db \\
//         ${vcf[0]} \\
//         -csvStats ${sample}.snpEff.csv \\
//         | bgzip -c > ${sample}.snpEff.vcf.gz
//     tabix -p vcf -f ${sample}.snpEff.vcf.gz
//     mv snpEff_summary.html ${sample}.snpEff.summary.html

//     SnpSift extractFields -s "," \\
//         -e "." \\
//         ${sample}.snpEff.vcf.gz \\
//         CHROM POS REF ALT \\
//         "ANN[*].GENE" "ANN[*].GENEID" \\
//         "ANN[*].IMPACT" "ANN[*].EFFECT" \\
//         "ANN[*].FEATURE" "ANN[*].FEATUREID" \\
//         "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \\
//         "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \\
//         "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \\
//         "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \\
//         "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \\
//         > ${sample}.snpSift.table.txt
//     	"""
// }

// ////////////////////////////////////////////////////
// /* --               UNICYCLER                  -- */
// ////////////////////////////////////////////////////

// /*
//  * STEP 6.3: De novo assembly with Unicycler
//  */
// process UNICYCLER {
//     tag "$sample"
//     label 'process_high'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/unicycler", mode: params.publish_dir_mode,
//     saveAs: { filename ->
//                   if (filename.endsWith(".png")) "bandage/$filename"
//                   else if (filename.endsWith(".svg")) "bandage/$filename"
//                   else filename
//             }

//     when:
//     !params.skip_assembly && 'unicycler' in assemblers

//     input:
//     tuple val(sample), val(single_end), path(reads) from ch_kraken2_unicycler

//     output:
//     tuple val(sample), val(single_end), path("*scaffolds.fa") into ch_unicycler_blast,
//                                                                    ch_unicycler_abacas,
//                                                                    ch_unicycler_plasmidid,
//                                                                    ch_unicycler_quast,
//                                                                    ch_unicycler_vg
//     path "*assembly.{gfa,png,svg}"

//     script:
//     input_reads = single_end ? "-s $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
//     """
//     unicycler \\
//         --threads $task.cpus \\
//         $input_reads \\
//         --out ./
//     mv assembly.fasta ${sample}.scaffolds.fa
//     mv assembly.gfa ${sample}.assembly.gfa

//     if [ -s ${sample}.assembly.gfa ]
//     then
//         Bandage image ${sample}.assembly.gfa ${sample}.assembly.png --height 1000
//         Bandage image ${sample}.assembly.gfa ${sample}.assembly.svg --height 1000
//     fi
//     """
// }

// /*
//  * STEP 6.3.1: Run Blast on MetaSPAdes de novo assembly
//  */
// process UNICYCLER_BLAST {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/unicycler/blast", mode: params.publish_dir_mode

//     when:
//     !params.skip_assembly && 'unicycler' in assemblers && !params.skip_blast

//     input:
//     tuple val(sample), val(single_end), path(scaffold) from ch_unicycler_blast
//     path db from ch_blast_db
//     path header from ch_blast_outfmt6_header

//     output:
//     path "*.blast*"

//     script:
//     """
//     blastn \\
//         -num_threads $task.cpus \\
//         -db $db/$fasta_base \\
//         -query $scaffold \\
//         -outfmt \'6 stitle std slen qlen qcovs\' \\
//         -out ${sample}.blast.txt

//     awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"}{print \$0,\$5/\$15,\$5/\$14}' ${sample}.blast.txt | awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"} \$15 > 200 && \$17 > 0.7 && \$1 !~ /phage/ {print \$0}' > ${sample}.blast.filt.txt
//     cat $header ${sample}.blast.filt.txt > ${sample}.blast.filt.header.txt
//     """
// }

// /*
//  * STEP 6.3.2: Run ABACAS on Unicycler de novo assembly
//  */
// process UNICYCLER_ABACAS {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/unicycler/abacas", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.indexOf("nucmer") > 0) "nucmer/$filename"
//                       else filename
//                 }

//     when:
//     !params.skip_assembly && 'unicycler' in assemblers && !params.skip_abacas

//     input:
//     tuple val(sample), val(single_end), path(scaffold) from ch_unicycler_abacas
//     path fasta from ch_fasta

//     output:
//     path "*.abacas*"

//     script:
//     """
//     abacas.pl -r $fasta -q $scaffold -m -p nucmer -o ${sample}.abacas
//     mv nucmer.delta ${sample}.abacas.nucmer.delta
//     mv nucmer.filtered.delta ${sample}.abacas.nucmer.filtered.delta
//     mv nucmer.tiling ${sample}.abacas.nucmer.tiling
//     mv unused_contigs.out ${sample}.abacas.unused.contigs.out
//     """
// }

// /*
//  * STEP 6.3.3: Run PlasmidID on Unicycler de novo assembly
//  */
// process UNICYCLER_PLASMIDID {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/unicycler/plasmidid", mode: params.publish_dir_mode

//     when:
//     !params.skip_assembly && 'unicycler' in assemblers && !params.skip_plasmidid

//     input:
//     tuple val(sample), val(single_end), path(scaffold) from ch_unicycler_plasmidid.filter { it.size() > 0 }
//     path fasta from ch_fasta

//     output:
//     path "$sample"

//     script:
//     """
//     plasmidID -d $fasta -s $sample -c $scaffold --only-reconstruct -C 47 -S 47 -i 60 --no-trim -o .
//     mv NO_GROUP/$sample ./$sample
//     """
// }

// /*
//  * STEP 6.3.4: Run Quast on Unicycler de novo assembly
//  */
// process UNICYCLER_QUAST {
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/unicycler", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (!filename.endsWith(".tsv")) filename
//                 }

//     when:
//     !params.skip_assembly && 'unicycler' in assemblers && !params.skip_assembly_quast

//     input:
//     path scaffolds from ch_unicycler_quast.collect{ it[2] }
//     path fasta from ch_fasta
//     path gff from ch_gff

//     output:
//     path "quast"
//     path "report.tsv" into ch_quast_unicycler_mqc

//     script:
//     features = params.gff ? "--features $gff" : ""
//     """
//     quast.py \\
//         --output-dir quast \\
//         -r $fasta \\
//         $features \\
//         --threads $task.cpus \\
//         ${scaffolds.join(' ')}
//     ln -s quast/report.tsv
//     """
// }

// /*
//  * STEP 6.3.5: Overlap scaffolds with Minimap2, induce and polish assembly, and call variants with seqwish and vg
//  */
// process UNICYCLER_VG {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/unicycler/variants", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".txt")) "bcftools_stats/$filename"
//                       else if (filename.endsWith(".png")) "bandage/$filename"
//                       else if (filename.endsWith(".svg")) "bandage/$filename"
//                       else filename
//                 }

//     when:
//     !params.skip_assembly && 'unicycler' in assemblers && !params.skip_vg

//     input:
//     tuple val(sample), val(single_end), path(scaffolds) from ch_unicycler_vg
//     path fasta from ch_fasta

//     output:
//     tuple val(sample), val(single_end), path("${sample}.vcf.gz*") into ch_unicycler_vg_vcf
//     path "*.bcftools_stats.txt" into ch_unicycler_vg_bcftools_mqc
//     path "*.{gfa,png,svg}"

//     script:
//     """
//     minimap2 -c -t $task.cpus -x asm20 $fasta $scaffolds > ${sample}.paf

//     cat $scaffolds $fasta > ${sample}.withRef.fasta
//     seqwish --paf-alns ${sample}.paf --seqs ${sample}.withRef.fasta --gfa ${sample}.gfa --threads $task.cpus

//     vg view -Fv ${sample}.gfa --threads $task.cpus > ${sample}.vg
//     vg convert -x ${sample}.vg > ${sample}.xg

//     samtools faidx $fasta
//     vg snarls ${sample}.xg > ${sample}.snarls
//     for chrom in `cat ${fasta}.fai | cut -f1`
//     do
//         vg deconstruct -p \$chrom ${sample}.xg -r ${sample}.snarls --threads $task.cpus \\
//             | bcftools sort -O v -T ./ \\
//             | bgzip -c > ${sample}.\$chrom.vcf.gz
//     done
//     bcftools concat --output-type z --output ${sample}.vcf.gz *.vcf.gz
//     tabix -p vcf -f ${sample}.vcf.gz
//     bcftools stats ${sample}.vcf.gz > ${sample}.bcftools_stats.txt

//     if [ -s ${sample}.gfa ]
//     then
//         Bandage image ${sample}.gfa ${sample}.png --height 1000
//         Bandage image ${sample}.gfa ${sample}.svg --height 1000
//     fi
//     """
// }

// /*
//  * STEP 6.3.6: Variant annotation with SnpEff and SnpSift
//  */
// process UNICYCLER_SNPEFF {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/unicycler/variants/snpeff", mode: params.publish_dir_mode

//     when:
//     !params.skip_assembly && 'unicycler' in assemblers && !params.skip_vg && params.gff && !params.skip_snpeff

//     input:
//     tuple val(sample), val(single_end), path(vcf) from ch_unicycler_vg_vcf
//     tuple file(db), file(config) from ch_snpeff_db_unicycler

//     output:
//     path "*.snpEff.csv" into ch_unicycler_snpeff_mqc
//     path "*.vcf.gz*"
//     path "*.{txt,html}"

//     script:
//     """
//     snpEff ${index_base} \\
//         -config $config \\
//         -dataDir $db \\
//         ${vcf[0]} \\
//         -csvStats ${sample}.snpEff.csv \\
//         | bgzip -c > ${sample}.snpEff.vcf.gz
//     tabix -p vcf -f ${sample}.snpEff.vcf.gz
//     mv snpEff_summary.html ${sample}.snpEff.summary.html

//     SnpSift extractFields -s "," \\
//         -e "." \\
//         ${sample}.snpEff.vcf.gz \\
//         CHROM POS REF ALT \\
//         "ANN[*].GENE" "ANN[*].GENEID" \\
//         "ANN[*].IMPACT" "ANN[*].EFFECT" \\
//         "ANN[*].FEATURE" "ANN[*].FEATUREID" \\
//         "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \\
//         "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \\
//         "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \\
//         "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \\
//         "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \\
//         > ${sample}.snpSift.table.txt
//     	"""
// }

// ////////////////////////////////////////////////////
// /* --                MINIA                     -- */
// ////////////////////////////////////////////////////

// /*
//  * STEP 6.3: De novo assembly with minia
//  */
// process MINIA {
//     tag "$sample"
//     label 'process_high'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/minia/${params.minia_kmer}", mode: params.publish_dir_mode

//     when:
//     !params.skip_assembly && 'minia' in assemblers

//     input:
//     tuple val(sample), val(single_end), path(reads) from ch_kraken2_minia

//     output:
//     tuple val(sample), val(single_end), path("*scaffolds.fa") into ch_minia_vg,
//                                                                    ch_minia_blast,
//                                                                    ch_minia_abacas,
//                                                                    ch_minia_plasmidid,
//                                                                    ch_minia_quast

//     script:
//     """
//     echo "${reads.join("\n")}" > input_files.txt
//     minia \\
//         -kmer-size $params.minia_kmer \\
//         -abundance-min 20 \\
//         -nb-cores $task.cpus \\
//         -in input_files.txt \\
//         -out ${sample}.k${params.minia_kmer}.a20
//     mv ${sample}.k${params.minia_kmer}.a20.contigs.fa ${sample}.k${params.minia_kmer}.scaffolds.fa
//     """
// }

// /*
//  * STEP 6.3.1: Run Blast on minia de novo assembly
//  */
// process MINIA_BLAST {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/minia/${params.minia_kmer}/blast", mode: params.publish_dir_mode

//     when:
//     !params.skip_assembly && 'minia' in assemblers && !params.skip_blast

//     input:
//     tuple val(sample), val(single_end), path(scaffold) from ch_minia_blast
//     path db from ch_blast_db
//     path header from ch_blast_outfmt6_header

//     output:
//     path "*.blast*"

//     script:
//     """
//     blastn \\
//         -num_threads $task.cpus \\
//         -db $db/$fasta_base \\
//         -query $scaffold \\
//         -outfmt \'6 stitle std slen qlen qcovs\' \\
//         -out ${sample}.blast.txt

//     awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"}{print \$0,\$5/\$15,\$5/\$14}' ${sample}.blast.txt | awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"} \$15 > 200 && \$17 > 0.7 && \$1 !~ /phage/ {print \$0}' > ${sample}.blast.filt.txt
//     cat $header ${sample}.blast.filt.txt > ${sample}.blast.filt.header.txt
//     """
// }

// /*
//  * STEP 6.3.2: Run ABACAS on minia de novo assembly
//  */
// process MINIA_ABACAS {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/minia/${params.minia_kmer}/abacas", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.indexOf("nucmer") > 0) "nucmer/$filename"
//                       else filename
//                 }

//     when:
//     !params.skip_assembly && 'minia' in assemblers && !params.skip_abacas

//     input:
//     tuple val(sample), val(single_end), path(scaffold) from ch_minia_abacas
//     path fasta from ch_fasta

//     output:
//     path "*.abacas*"

//     script:
//     """
//     abacas.pl -r $fasta -q $scaffold -m -p nucmer -o ${sample}.abacas
//     mv nucmer.delta ${sample}.abacas.nucmer.delta
//     mv nucmer.filtered.delta ${sample}.abacas.nucmer.filtered.delta
//     mv nucmer.tiling ${sample}.abacas.nucmer.tiling
//     mv unused_contigs.out ${sample}.abacas.unused.contigs.out
//     """
// }

// /*
//  * STEP 6.3.3: Run PlasmidID on minia de novo assembly
//  */
// process MINIA_PLASMIDID {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/minia/${params.minia_kmer}/plasmidid", mode: params.publish_dir_mode

//     when:
//     !params.skip_assembly && 'minia' in assemblers && !params.skip_plasmidid

//     input:
//     tuple val(sample), val(single_end), path(scaffold) from ch_minia_plasmidid.filter { it.size() > 0 }
//     path fasta from ch_fasta

//     output:
//     path "$sample"

//     script:
//     """
//     plasmidID -d $fasta -s $sample -c $scaffold --only-reconstruct -C 47 -S 47 -i 60 --no-trim -o .
//     mv NO_GROUP/$sample ./$sample
//     """
// }

// /*
//  * STEP 6.3.4: Run Quast on minia de novo assembly
//  */
// process MINIA_QUAST {
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/minia/${params.minia_kmer}", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (!filename.endsWith(".tsv")) filename
//                 }

//     when:
//     !params.skip_assembly && 'minia' in assemblers && !params.skip_assembly_quast

//     input:
//     path scaffolds from ch_minia_quast.collect{ it[2] }
//     path fasta from ch_fasta
//     path gff from ch_gff

//     output:
//     path "quast"
//     path "report.tsv" into ch_quast_minia_mqc

//     script:
//     features = params.gff ? "--features $gff" : ""
//     """
//     quast.py \\
//         --output-dir quast \\
//         -r $fasta \\
//         $features \\
//         --threads $task.cpus \\
//         ${scaffolds.join(' ')}
//     ln -s quast/report.tsv
//     """
// }

// /*
//  * STEP 6.3.5: Overlap scaffolds with Minimap2, induce and polish assembly, and call variants with seqwish and vg
//  */
// process MINIA_VG {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/minia/${params.minia_kmer}/variants", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".txt")) "bcftools_stats/$filename"
//                       else if (filename.endsWith(".png")) "bandage/$filename"
//                       else if (filename.endsWith(".svg")) "bandage/$filename"
//                       else filename
//                 }

//     when:
//     !params.skip_assembly && 'minia' in assemblers && !params.skip_vg

//     input:
//     tuple val(sample), val(single_end), path(scaffolds) from ch_minia_vg
//     path fasta from ch_fasta

//     output:
//     tuple val(sample), val(single_end), path("${sample}.vcf.gz*") into ch_minia_vg_vcf
//     path "*.bcftools_stats.txt" into ch_minia_vg_bcftools_mqc
//     path "*.{gfa,png,svg}"

//     script:
//     """
//     minimap2 -c -t $task.cpus -x asm20 $fasta $scaffolds > ${sample}.paf

//     cat $scaffolds $fasta > ${sample}.withRef.fasta
//     seqwish --paf-alns ${sample}.paf --seqs ${sample}.withRef.fasta --gfa ${sample}.gfa --threads $task.cpus

//     vg view -Fv ${sample}.gfa --threads $task.cpus > ${sample}.vg
//     vg convert -x ${sample}.vg > ${sample}.xg

//     samtools faidx $fasta
//     vg snarls ${sample}.xg > ${sample}.snarls
//     for chrom in `cat ${fasta}.fai | cut -f1`
//     do
//         vg deconstruct -p \$chrom ${sample}.xg -r ${sample}.snarls --threads $task.cpus \\
//             | bcftools sort -O v -T ./ \\
//             | bgzip -c > ${sample}.\$chrom.vcf.gz
//     done
//     bcftools concat --output-type z --output ${sample}.vcf.gz *.vcf.gz
//     tabix -p vcf -f ${sample}.vcf.gz
//     bcftools stats ${sample}.vcf.gz > ${sample}.bcftools_stats.txt

//     if [ -s ${sample}.gfa ]
//     then
//         Bandage image ${sample}.gfa ${sample}.png --height 1000
//         Bandage image ${sample}.gfa ${sample}.svg --height 1000
//     fi
//     """
// }

// /*
//  * STEP 6.3.6: Variant annotation with SnpEff and SnpSift
//  */
// process MINIA_SNPEFF {
//     tag "$sample"
//     label 'process_medium'
//     label 'error_ignore'
//     publishDir "${params.outdir}/assembly/minia/${params.minia_kmer}/variants/snpeff", mode: params.publish_dir_mode

//     when:
//     !params.skip_assembly && 'minia' in assemblers && !params.skip_vg && params.gff && !params.skip_snpeff

//     input:
//     tuple val(sample), val(single_end), path(vcf) from ch_minia_vg_vcf
//     tuple file(db), file(config) from ch_snpeff_db_minia

//     output:
//     path "*.snpEff.csv" into ch_minia_snpeff_mqc
//     path "*.vcf.gz*"
//     path "*.{txt,html}"

//     script:
//     """
//     snpEff ${index_base} \\
//         -config $config \\
//         -dataDir $db \\
//         ${vcf[0]} \\
//         -csvStats ${sample}.snpEff.csv \\
//         | bgzip -c > ${sample}.snpEff.vcf.gz
//     tabix -p vcf -f ${sample}.snpEff.vcf.gz
//     mv snpEff_summary.html ${sample}.snpEff.summary.html

//     SnpSift extractFields -s "," \\
//         -e "." \\
//         ${sample}.snpEff.vcf.gz \\
//         CHROM POS REF ALT \\
//         "ANN[*].GENE" "ANN[*].GENEID" \\
//         "ANN[*].IMPACT" "ANN[*].EFFECT" \\
//         "ANN[*].FEATURE" "ANN[*].FEATUREID" \\
//         "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \\
//         "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \\
//         "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \\
//         "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \\
//         "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \\
//         > ${sample}.snpSift.table.txt
//     	"""
// }

// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                          MULTIQC                                    -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////

// /*
//  * Parse software version numbers
//  */
// process get_software_versions {
//     publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith(".csv")) filename
//                       else null
//                 }

//     output:
//     path "software_versions_mqc.yaml" into ch_software_versions_yaml
//     path "software_versions.csv"

//     script:
//     """
//     echo $workflow.manifest.version > v_pipeline.txt
//     echo $workflow.nextflow.version > v_nextflow.txt
//     parallel-fastq-dump --version > v_parallel_fastq_dump.txt
//     fastqc --version > v_fastqc.txt
//     fastp --version 2> v_fastp.txt
//     bowtie2 --version > v_bowtie2.txt
//     samtools --version > v_samtools.txt
//     bedtools --version > v_bedtools.txt
//     mosdepth --version > v_mosdepth.txt
//     picard CollectMultipleMetrics --version &> v_picard.txt || true
//     ivar -v > v_ivar.txt
//     echo \$(varscan 2>&1) > v_varscan.txt
//     bcftools -v > v_bcftools.txt
//     snpEff -version > v_snpeff.txt
//     echo \$(SnpSift 2>&1) > v_snpsift.txt
//     quast.py --version > v_quast.txt
//     cutadapt --version > v_cutadapt.txt
//     kraken2 --version > v_kraken2.txt
//     spades.py --version > v_spades.txt
//     unicycler --version > v_unicycler.txt
//     minia --version > v_minia.txt
//     blastn -version > v_blast.txt
//     abacas.pl -v &> v_abacas.txt || true
//     plasmidID -v > v_plasmidid.txt  || true
//     Bandage --version > v_bandage.txt
//     minimap2 --version > v_minimap2.txt
//     vg version > v_vg.txt
//     echo \$(R --version 2>&1) > v_R.txt
//     multiqc --version > v_multiqc.txt
//     scrape_software_versions.py &> software_versions_mqc.yaml
//     """
// }

// /*
//  * STEP 7: MultiQC
//  */
// process MULTIQC {
//     label 'process_medium'
//     publishDir "${params.outdir}", mode: params.publish_dir_mode,
//         saveAs: { filename ->
//                       if (filename.endsWith("assembly_metrics_mqc.tsv")) "assembly/$filename"
//                       else if (filename.endsWith("variants_metrics_mqc.tsv")) "variants/$filename"
//                       else "multiqc/$filename"
//                 }

//     when:
//     !params.skip_multiqc

//     input:
//     path (multiqc_config) from ch_multiqc_config
//     path (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
//     path ('fastqc/*') from ch_fastqc_raw_reports_mqc.collect().ifEmpty([])
//     path ('fastp/log/*') from ch_fastp_mqc.collect().ifEmpty([])
//     path ('fastp/fastqc/*') from ch_fastp_fastqc_mqc.collect().ifEmpty([])
//     path ('bowtie2/log/*') from ch_bowtie2_mqc.collect().ifEmpty([])
//     path ('bowtie2/flagstat/*') from ch_sort_bam_flagstat_mqc.collect().ifEmpty([])
//     path ('ivar/trim/flagstat/*') from ch_ivar_trim_flagstat_mqc.collect().ifEmpty([])
//     path ('ivar/trim/log/*') from ch_ivar_trim_log_mqc.collect().ifEmpty([])
//     path ('picard/markdup/*') from ch_markdup_bam_flagstat_mqc.collect().ifEmpty([])
//     path ('picard/metrics/*') from ch_markdup_bam_metrics_mqc.collect().ifEmpty([])
//     path ('picard/metrics/*') from ch_picard_metrics_mqc.collect().ifEmpty([])
//     path ('mosdepth/genome/*') from ch_mosdepth_genome_mqc.collect().ifEmpty([])
//     path ('varscan2/counts/lowfreq/*') from ch_varscan2_log_mqc.collect().ifEmpty([])
//     path ('varscan2/bcftools/highfreq/*') from ch_varscan2_bcftools_highfreq_mqc.collect().ifEmpty([])
//     path ('varscan2/snpeff/highfreq/*') from ch_varscan2_snpeff_highfreq_mqc.collect().ifEmpty([])
//     path ('varscan2/quast/highfreq/*') from ch_varscan2_quast_mqc.collect().ifEmpty([])
//     path ('ivar/variants/counts/lowfreq/*') from ch_ivar_count_mqc.collect().ifEmpty([])
//     path ('ivar/variants/bcftools/highfreq/*') from ch_ivar_bcftools_highfreq_mqc.collect().ifEmpty([])
//     path ('ivar/variants/snpeff/highfreq/*') from ch_ivar_snpeff_highfreq_mqc.collect().ifEmpty([])
//     path ('ivar/consensus/quast/highfreq/*') from ch_ivar_quast_mqc.collect().ifEmpty([])
//     path ('bcftools/variants/bcftools/*') from ch_bcftools_variants_mqc.collect().ifEmpty([])
//     path ('bcftools/variants/snpeff/*') from ch_bcftools_snpeff_mqc.collect().ifEmpty([])
//     path ('bcftools/consensus/quast/*') from ch_bcftools_quast_mqc.collect().ifEmpty([])
//     path ('cutadapt/log/*') from ch_cutadapt_mqc.collect().ifEmpty([])
//     path ('cutadapt/fastqc/*') from ch_cutadapt_fastqc_mqc.collect().ifEmpty([])
//     path ('kraken2/*') from ch_kraken2_report_mqc.collect().ifEmpty([])
//     path ('spades/bcftools/*') from ch_spades_vg_bcftools_mqc.collect().ifEmpty([])
//     path ('spades/snpeff/*') from ch_spades_snpeff_mqc.collect().ifEmpty([])
//     path ('spades/quast/*') from ch_quast_spades_mqc.collect().ifEmpty([])
//     path ('metaspades/bcftools/*') from ch_metaspades_vg_bcftools_mqc.collect().ifEmpty([])
//     path ('metaspades/snpeff/*') from ch_metaspades_snpeff_mqc.collect().ifEmpty([])
//     path ('metaspades/quast/*') from ch_quast_metaspades_mqc.collect().ifEmpty([])
//     path ('unicycler/bcftools/*') from ch_unicycler_vg_bcftools_mqc.collect().ifEmpty([])
//     path ('unicycler/snpeff/*') from ch_unicycler_snpeff_mqc.collect().ifEmpty([])
//     path ('unicycler/quast/*') from ch_quast_unicycler_mqc.collect().ifEmpty([])
//     path ('minia/bcftools/*') from ch_minia_vg_bcftools_mqc.collect().ifEmpty([])
//     path ('minia/snpeff/*') from ch_minia_snpeff_mqc.collect().ifEmpty([])
//     path ('minia/quast/*') from ch_quast_minia_mqc.collect().ifEmpty([])
//     path ('software_versions/*') from ch_software_versions_yaml.collect()
//     path workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

//     output:
//     path "*multiqc_report.html" into ch_multiqc_report
//     path "*_data"
//     path "*.tsv"

//     script:
//     rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
//     rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
//     custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
//     """
//     multiqc . -f $rtitle $rfilename $custom_config_file
//     multiqc_to_custom_tsv.py
//     multiqc . -f $rtitle $rfilename $custom_config_file
//     """
// }
