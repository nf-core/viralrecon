/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def valid_params = [
    protocols   : ['metagenomic', 'amplicon'],
    callers     : ['ivar', 'bcftools'],
    assemblers  : ['spades', 'unicycler', 'minia'],
    spades_modes: ['rnaviral', 'corona', 'metaviral', 'meta', 'metaplasmid', 'plasmid', 'isolate', 'rna', 'bio']
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowIllumina.initialise(params, log, valid_params)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.fasta, params.gff, params.bowtie2_index,
    params.kraken2_db, params.primer_bed, params.primer_fasta,
    params.blast_db, params.spades_hmm, params.multiqc_config
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

if (params.input)      { ch_input      = file(params.input)      } else { exit 1, 'Input samplesheet file not specified!' }
if (params.spades_hmm) { ch_spades_hmm = file(params.spades_hmm) } else { ch_spades_hmm = ch_dummy_file                   }

def assemblers = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []
def callers    = params.callers    ? params.callers.split(',').collect{ it.trim().toLowerCase() }    : []
if (!callers)  { callers = params.protocol == 'amplicon' ? ['ivar'] : ['bcftools'] }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config_illumina.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

// Header files
ch_blast_outfmt6_header     = file("$projectDir/assets/headers/blast_outfmt6_header.txt", checkIfExists: true)
ch_ivar_variants_header_mqc = file("$projectDir/assets/headers/ivar_variants_header_mqc.txt", checkIfExists: true)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def multiqc_options   = modules['illumina_multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

if (!params.skip_assembly) {
    multiqc_options.publish_files.put('assembly_metrics_mqc.csv','')
}
if (!params.skip_variants) {
    multiqc_options.publish_files.put('variants_metrics_mqc.csv','')
}

include { BCFTOOLS_ISEC         } from '../modules/local/bcftools_isec'         addParams( options: modules['illumina_bcftools_isec'] )
include { CUTADAPT              } from '../modules/local/cutadapt'              addParams( options: modules['illumina_cutadapt']      )
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files: ['tsv':'']]       )
include { MULTIQC               } from '../modules/local/multiqc_illumina'      addParams( options: multiqc_options                   )
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_GENOME      } from '../modules/local/plot_mosdepth_regions' addParams( options: modules['illumina_plot_mosdepth_regions_genome']   )
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_AMPLICON    } from '../modules/local/plot_mosdepth_regions' addParams( options: modules['illumina_plot_mosdepth_regions_amplicon'] )
include { MULTIQC_CUSTOM_TSV_FROM_STRING as MULTIQC_CUSTOM_TSV_FAIL_READS         } from '../modules/local/multiqc_custom_tsv_from_string' addParams( options: [publish_files: false] )
include { MULTIQC_CUSTOM_TSV_FROM_STRING as MULTIQC_CUSTOM_TSV_FAIL_MAPPED        } from '../modules/local/multiqc_custom_tsv_from_string' addParams( options: [publish_files: false] )
include { MULTIQC_CUSTOM_TSV_FROM_STRING as MULTIQC_CUSTOM_TSV_IVAR_NEXTCLADE     } from '../modules/local/multiqc_custom_tsv_from_string' addParams( options: [publish_files: false] )
include { MULTIQC_CUSTOM_TSV_FROM_STRING as MULTIQC_CUSTOM_TSV_BCFTOOLS_NEXTCLADE } from '../modules/local/multiqc_custom_tsv_from_string' addParams( options: [publish_files: false] )
include { MULTIQC_CUSTOM_CSV_FROM_MAP as MULTIQC_CUSTOM_CSV_IVAR_PANGOLIN     } from '../modules/local/multiqc_custom_csv_from_map'    addParams( options: [publish_files: false] )
include { MULTIQC_CUSTOM_CSV_FROM_MAP as MULTIQC_CUSTOM_CSV_BCFTOOLS_PANGOLIN } from '../modules/local/multiqc_custom_csv_from_map'    addParams( options: [publish_files: false] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
def publish_genome_options    = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options     = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]
def publish_db_options        = params.save_reference ? [publish_dir: 'genome/db']    : [publish_files: false]
def bedtools_getfasta_options = modules['illumina_bedtools_getfasta']
def bowtie2_build_options     = modules['illumina_bowtie2_build']
def snpeff_build_options      = modules['illumina_snpeff_build']
def makeblastdb_options       = modules['illumina_blast_makeblastdb']
def kraken2_build_options     = modules['illumina_kraken2_build']
def collapse_primers_options  = modules['illumina_collapse_primers_illumina']
if (!params.save_reference) {
    bedtools_getfasta_options['publish_files'] = false
    bowtie2_build_options['publish_files']     = false
    snpeff_build_options['publish_files']      = false
    makeblastdb_options['publish_files']       = false
    kraken2_build_options['publish_files']     = false
    collapse_primers_options['publish_files']  = false
}

def ivar_trim_options   = modules['illumina_ivar_trim']
ivar_trim_options.args += params.ivar_trim_noprimer ? '' : Utils.joinModuleArgs(['-e'])
ivar_trim_options.args += params.ivar_trim_offset   ? Utils.joinModuleArgs(["-x ${params.ivar_trim_offset}"]) : ''

def ivar_trim_sort_bam_options = modules['illumina_ivar_trim_sort_bam']
if (params.skip_markduplicates) {
    ivar_trim_sort_bam_options.publish_files.put('bam','')
    ivar_trim_sort_bam_options.publish_files.put('bai','')
}

def spades_options   = modules['illumina_spades']
spades_options.args += params.spades_mode ? Utils.joinModuleArgs(["--${params.spades_mode}"]) : ''

include { INPUT_CHECK        } from '../subworkflows/local/input_check'             addParams( options: [:] )
include { PREPARE_GENOME     } from '../subworkflows/local/prepare_genome_illumina' addParams( genome_options: publish_genome_options, index_options: publish_index_options, db_options: publish_db_options, bowtie2_build_options: bowtie2_build_options, bedtools_getfasta_options: bedtools_getfasta_options, collapse_primers_options: collapse_primers_options, snpeff_build_options: snpeff_build_options, makeblastdb_options: makeblastdb_options, kraken2_build_options: kraken2_build_options )
include { PRIMER_TRIM_IVAR   } from '../subworkflows/local/primer_trim_ivar'        addParams( ivar_trim_options: ivar_trim_options, samtools_options: ivar_trim_sort_bam_options )
include { VARIANTS_IVAR      } from '../subworkflows/local/variants_ivar'           addParams( ivar_variants_options: modules['illumina_ivar_variants'], ivar_variants_to_vcf_options: modules['illumina_ivar_variants_to_vcf'], tabix_bgzip_options: modules['illumina_ivar_tabix_bgzip'], tabix_tabix_options: modules['illumina_ivar_tabix_tabix'], bcftools_stats_options: modules['illumina_ivar_bcftools_stats'], ivar_consensus_options: modules['illumina_ivar_consensus'], consensus_plot_options: modules['illumina_ivar_consensus_plot'], quast_options: modules['illumina_ivar_quast'], snpeff_options: modules['illumina_ivar_snpeff'], snpsift_options: modules['illumina_ivar_snpsift'], snpeff_bgzip_options: modules['illumina_ivar_snpeff_bgzip'], snpeff_tabix_options: modules['illumina_ivar_snpeff_tabix'], snpeff_stats_options: modules['illumina_ivar_snpeff_stats'], pangolin_options: modules['illumina_ivar_pangolin'], nextclade_options: modules['illumina_ivar_nextclade'], asciigenome_options: modules['illumina_ivar_asciigenome'] )
include { VARIANTS_BCFTOOLS  } from '../subworkflows/local/variants_bcftools'       addParams( bcftools_mpileup_options: modules['illumina_bcftools_mpileup'], quast_options: modules['illumina_bcftools_quast'], consensus_genomecov_options: modules['illumina_bcftools_consensus_genomecov'], consensus_merge_options: modules['illumina_bcftools_consensus_merge'], consensus_mask_options: modules['illumina_bcftools_consensus_mask'], consensus_maskfasta_options: modules['illumina_bcftools_consensus_maskfasta'], consensus_bcftools_options: modules['illumina_bcftools_consensus_bcftools'], consensus_plot_options: modules['illumina_bcftools_consensus_plot'], snpeff_options: modules['illumina_bcftools_snpeff'], snpsift_options: modules['illumina_bcftools_snpsift'], snpeff_bgzip_options: modules['illumina_bcftools_snpeff_bgzip'], snpeff_tabix_options: modules['illumina_bcftools_snpeff_tabix'], snpeff_stats_options: modules['illumina_bcftools_snpeff_stats'], pangolin_options: modules['illumina_bcftools_pangolin'], nextclade_options: modules['illumina_bcftools_nextclade'], asciigenome_options: modules['illumina_bcftools_asciigenome'] )
include { ASSEMBLY_SPADES    } from '../subworkflows/local/assembly_spades'         addParams( spades_options: spades_options, bandage_options: modules['illumina_spades_bandage'], blastn_options: modules['illumina_spades_blastn'], blastn_filter_options: modules['illumina_spades_blastn_filter'], abacas_options: modules['illumina_spades_abacas'], plasmidid_options: modules['illumina_spades_plasmidid'], quast_options: modules['illumina_spades_quast'] )
include { ASSEMBLY_UNICYCLER } from '../subworkflows/local/assembly_unicycler'      addParams( unicycler_options: modules['illumina_unicycler'], bandage_options: modules['illumina_unicycler_bandage'], blastn_options: modules['illumina_unicycler_blastn'], blastn_filter_options: modules['illumina_unicycler_blastn_filter'], abacas_options: modules['illumina_unicycler_abacas'], plasmidid_options: modules['illumina_unicycler_plasmidid'], quast_options: modules['illumina_unicycler_quast'] )
include { ASSEMBLY_MINIA     } from '../subworkflows/local/assembly_minia'          addParams( minia_options: modules['illumina_minia'], blastn_options: modules['illumina_minia_blastn'], blastn_filter_options: modules['illumina_minia_blastn_filter'], abacas_options: modules['illumina_minia_abacas'], plasmidid_options: modules['illumina_minia_plasmidid'], quast_options: modules['illumina_minia_quast'] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                     } from '../modules/nf-core/modules/cat/fastq/main'                     addParams( options: modules['illumina_cat_fastq']                     )
include { FASTQC                        } from '../modules/nf-core/modules/fastqc/main'                        addParams( options: modules['illumina_cutadapt_fastqc']               )
include { KRAKEN2_KRAKEN2               } from '../modules/nf-core/modules/kraken2/kraken2/main'               addParams( options: modules['illumina_kraken2_kraken2']               )
include { PICARD_COLLECTMULTIPLEMETRICS } from '../modules/nf-core/modules/picard/collectmultiplemetrics/main' addParams( options: modules['illumina_picard_collectmultiplemetrics'] )
include { MOSDEPTH as MOSDEPTH_GENOME   } from '../modules/nf-core/modules/mosdepth/main'                      addParams( options: modules['illumina_mosdepth_genome']               )
include { MOSDEPTH as MOSDEPTH_AMPLICON } from '../modules/nf-core/modules/mosdepth/main'                      addParams( options: modules['illumina_mosdepth_amplicon']             )

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
def fastp_options = modules['illumina_fastp']
if (params.save_trimmed_fail) { fastp_options.publish_files.put('fail.fastq.gz','') }

def bowtie2_align_options = modules['illumina_bowtie2_align']
if (params.save_unaligned) { bowtie2_align_options.publish_files.put('fastq.gz','unmapped') }

def markduplicates_options   = modules['illumina_picard_markduplicates']
markduplicates_options.args += params.filter_duplicates ?  Utils.joinModuleArgs(['REMOVE_DUPLICATES=true']) : ''

include { FASTQC_FASTP           } from '../subworkflows/nf-core/fastqc_fastp'           addParams( fastqc_raw_options: modules['illumina_fastqc_raw'], fastqc_trim_options: modules['illumina_fastqc_trim'], fastp_options: fastp_options )
include { ALIGN_BOWTIE2          } from '../subworkflows/nf-core/align_bowtie2'          addParams( align_options: bowtie2_align_options, samtools_options: modules['illumina_bowtie2_sort_bam']                                           )
include { MARK_DUPLICATES_PICARD } from '../subworkflows/nf-core/mark_duplicates_picard' addParams( markduplicates_options: markduplicates_options, samtools_options: modules['illumina_picard_markduplicates_sort_bam']                   )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report    = []
def pass_mapped_reads = [:]
def fail_mapped_reads = [:]

workflow ILLUMINA {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME (
        ch_dummy_file
    )

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

        // Check if the primer BED file supplied to the pipeline is from the SWIFT/SNAP protocol
        if (!params.ivar_trim_offset) {
            PREPARE_GENOME
                .out
                .primer_bed
                .map { WorkflowIllumina.checkIfSwiftProtocol(it, 'covid19genome', log) }
        }
    }

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input,
        params.platform
    )
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ]
    }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_FASTP (
        ch_cat_fastq
    )
    ch_variants_fastq    = FASTQC_FASTP.out.reads
    ch_software_versions = ch_software_versions.mix(FASTQC_FASTP.out.fastqc_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_FASTP.out.fastp_version.first().ifEmpty(null))

    //
    // Filter empty FastQ files after adapter trimming
    //
    ch_fail_reads_multiqc = Channel.empty()
    if (!params.skip_fastp) {
        ch_variants_fastq
            .join(FASTQC_FASTP.out.trim_json)
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
            .set { ch_pass_fail_reads }

        MULTIQC_CUSTOM_TSV_FAIL_READS (
            ch_pass_fail_reads.collect(),
            'Sample\tReads before trimming',
            'fail_mapped_reads'
        )
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
            PREPARE_GENOME.out.kraken2_db
        )
        ch_kraken2_multiqc   = KRAKEN2_KRAKEN2.out.txt
        ch_software_versions = ch_software_versions.mix(KRAKEN2_KRAKEN2.out.version.first().ifEmpty(null))

        if (params.kraken2_variants_host_filter) {
            ch_variants_fastq = KRAKEN2_KRAKEN2.out.unclassified
        }

        if (params.kraken2_assembly_host_filter) {
            ch_assembly_fastq = KRAKEN2_KRAKEN2.out.unclassified
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
        ALIGN_BOWTIE2 (
            ch_variants_fastq,
            PREPARE_GENOME.out.bowtie2_index
        )
        ch_bam                      = ALIGN_BOWTIE2.out.bam
        ch_bai                      = ALIGN_BOWTIE2.out.bai
        ch_bowtie2_multiqc          = ALIGN_BOWTIE2.out.log_out
        ch_bowtie2_flagstat_multiqc = ALIGN_BOWTIE2.out.flagstat
        ch_software_versions = ch_software_versions.mix(ALIGN_BOWTIE2.out.bowtie2_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ALIGN_BOWTIE2.out.samtools_version.first().ifEmpty(null))
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

        MULTIQC_CUSTOM_TSV_FAIL_MAPPED (
            ch_pass_fail_mapped.fail.collect(),
            'Sample\tMapped reads',
            'fail_mapped_samples'
        )
        .set { ch_fail_mapping_multiqc }
    }

    //
    // SUBWORKFLOW: Trim primer sequences from reads with iVar
    //
    ch_ivar_trim_flagstat_multiqc = Channel.empty()
    if (!params.skip_variants && !params.skip_ivar_trim && params.protocol == 'amplicon') {
        PRIMER_TRIM_IVAR (
            ch_bam.join(ch_bai, by: [0]),
            PREPARE_GENOME.out.primer_bed
        )
        ch_bam                        = PRIMER_TRIM_IVAR.out.bam
        ch_bai                        = PRIMER_TRIM_IVAR.out.bai
        ch_ivar_trim_flagstat_multiqc = PRIMER_TRIM_IVAR.out.flagstat
        ch_software_versions  = ch_software_versions.mix(PRIMER_TRIM_IVAR.out.ivar_version.first().ifEmpty(null))
    }

    //
    // SUBWORKFLOW: Mark duplicate reads
    //
    ch_markduplicates_flagstat_multiqc = Channel.empty()
    if (!params.skip_variants && !params.skip_markduplicates) {
        MARK_DUPLICATES_PICARD (
            ch_bam
        )
        ch_bam                             = MARK_DUPLICATES_PICARD.out.bam
        ch_bai                             = MARK_DUPLICATES_PICARD.out.bai
        ch_markduplicates_flagstat_multiqc = MARK_DUPLICATES_PICARD.out.flagstat
        ch_software_versions = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.picard_version.first().ifEmpty(null))
    }

    //
    // MODULE: Picard metrics
    //
    if (!params.skip_variants && !params.skip_picard_metrics) {
        PICARD_COLLECTMULTIPLEMETRICS (
            ch_bam,
            PREPARE_GENOME.out.fasta
        )
        ch_software_versions = ch_software_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.version.first().ifEmpty(null))
    }

    //
    // MODULE: Genome-wide and amplicon-specific coverage QC plots
    //
    ch_mosdepth_multiqc         = Channel.empty()
    ch_amplicon_heatmap_multiqc = Channel.empty()
    if (!params.skip_variants && !params.skip_mosdepth) {

        MOSDEPTH_GENOME (
            ch_bam.join(ch_bai, by: [0]),
            ch_dummy_file,
            200
        )
        ch_mosdepth_multiqc  = MOSDEPTH_GENOME.out.global_txt
        ch_software_versions = ch_software_versions.mix(MOSDEPTH_GENOME.out.version.first().ifEmpty(null))

        PLOT_MOSDEPTH_REGIONS_GENOME (
            MOSDEPTH_GENOME.out.regions_bed.collect { it[1] }
        )

        if (params.protocol == 'amplicon') {
            MOSDEPTH_AMPLICON (
                ch_bam.join(ch_bai, by: [0]),
                PREPARE_GENOME.out.primer_collapsed_bed,
                0
            )

            PLOT_MOSDEPTH_REGIONS_AMPLICON (
                MOSDEPTH_AMPLICON.out.regions_bed.collect { it[1] }
            )
            ch_amplicon_heatmap_multiqc = PLOT_MOSDEPTH_REGIONS_AMPLICON.out.heatmap_tsv
        }
    }

    //
    // SUBWORKFLOW: Call variants with IVar
    //
    ch_ivar_vcf               = Channel.empty()
    ch_ivar_tbi               = Channel.empty()
    ch_ivar_counts_multiqc    = Channel.empty()
    ch_ivar_stats_multiqc     = Channel.empty()
    ch_ivar_snpeff_multiqc    = Channel.empty()
    ch_ivar_quast_multiqc     = Channel.empty()
    ch_ivar_pangolin_multiqc  = Channel.empty()
    ch_ivar_nextclade_multiqc = Channel.empty()
    if (!params.skip_variants && 'ivar' in callers) {
        VARIANTS_IVAR (
            ch_bam,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.chrom_sizes,
            params.gff ? PREPARE_GENOME.out.gff : [],
            (params.protocol == 'amplicon' && params.primer_bed) ? PREPARE_GENOME.out.primer_bed : [],
            PREPARE_GENOME.out.snpeff_db,
            PREPARE_GENOME.out.snpeff_config,
            ch_ivar_variants_header_mqc
        )
        ch_ivar_vcf              = VARIANTS_IVAR.out.vcf
        ch_ivar_tbi              = VARIANTS_IVAR.out.tbi
        ch_ivar_counts_multiqc   = VARIANTS_IVAR.out.multiqc_tsv
        ch_ivar_stats_multiqc    = VARIANTS_IVAR.out.stats
        ch_ivar_snpeff_multiqc   = VARIANTS_IVAR.out.snpeff_csv
        ch_ivar_quast_multiqc    = VARIANTS_IVAR.out.quast_tsv
        ch_ivar_pangolin_report  = VARIANTS_IVAR.out.pangolin_report
        ch_ivar_nextclade_report = VARIANTS_IVAR.out.nextclade_report
        ch_software_versions = ch_software_versions.mix(VARIANTS_IVAR.out.ivar_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(VARIANTS_IVAR.out.tabix_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(VARIANTS_IVAR.out.bcftools_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(VARIANTS_IVAR.out.quast_version.ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(VARIANTS_IVAR.out.snpeff_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(VARIANTS_IVAR.out.snpsift_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(VARIANTS_IVAR.out.pangolin_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(VARIANTS_IVAR.out.nextclade_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(VARIANTS_IVAR.out.asciigenome_version.first().ifEmpty(null))

        //
        // MODULE: Get Pangolin lineage information for MultiQC report
        //
        ch_ivar_pangolin_report
            .map { meta, report ->
                def fields = WorkflowCommons.getPangolinFieldMap(report)
                return [sample:meta.id] << fields
            }
            .set { ch_ivar_pangolin_multiqc }

        MULTIQC_CUSTOM_CSV_IVAR_PANGOLIN (
            ch_ivar_pangolin_multiqc.collect(),
            'ivar_pangolin_lineage'
        )
        .set { ch_ivar_pangolin_multiqc }

        //
        // MODULE: Get Nextclade clade information for MultiQC report
        //
        ch_ivar_nextclade_report
            .map { meta, csv ->
                def clade = WorkflowCommons.getNextcladeFieldMapFromCsv(csv)['clade']
                return [ "$meta.id\t$clade" ]
            }
            .set { ch_ivar_nextclade_multiqc }

        MULTIQC_CUSTOM_TSV_IVAR_NEXTCLADE (
            ch_ivar_nextclade_multiqc.collect(),
            'Sample\tclade',
            'ivar_nextclade_clade'
        )
        .set { ch_ivar_nextclade_multiqc }
    }

    //
    // SUBWORKFLOW: Call variants with BCFTools
    //
    ch_bcftools_vcf               = Channel.empty()
    ch_bcftools_tbi               = Channel.empty()
    ch_bcftools_stats_multiqc     = Channel.empty()
    ch_bcftools_snpeff_multiqc    = Channel.empty()
    ch_bcftools_quast_multiqc     = Channel.empty()
    ch_bcftools_pangolin_multiqc  = Channel.empty()
    ch_bcftools_nextclade_multiqc = Channel.empty()
    if (!params.skip_variants && 'bcftools' in callers) {
        VARIANTS_BCFTOOLS (
            ch_bam,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.chrom_sizes,
            params.gff ? PREPARE_GENOME.out.gff : [],
            (params.protocol == 'amplicon' && params.primer_bed) ? PREPARE_GENOME.out.primer_bed : [],
            PREPARE_GENOME.out.snpeff_db,
            PREPARE_GENOME.out.snpeff_config
        )
        ch_bcftools_vcf              = VARIANTS_BCFTOOLS.out.vcf
        ch_bcftools_tbi              = VARIANTS_BCFTOOLS.out.tbi
        ch_bcftools_stats_multiqc    = VARIANTS_BCFTOOLS.out.stats
        ch_bcftools_snpeff_multiqc   = VARIANTS_BCFTOOLS.out.snpeff_csv
        ch_bcftools_quast_multiqc    = VARIANTS_BCFTOOLS.out.quast_tsv
        ch_bcftools_pangolin_report  = VARIANTS_BCFTOOLS.out.pangolin_report
        ch_bcftools_nextclade_report = VARIANTS_BCFTOOLS.out.nextclade_report
        ch_software_versions = ch_software_versions.mix(VARIANTS_BCFTOOLS.out.bcftools_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(VARIANTS_BCFTOOLS.out.bedtools_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(VARIANTS_BCFTOOLS.out.quast_version.ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(VARIANTS_BCFTOOLS.out.snpeff_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(VARIANTS_BCFTOOLS.out.snpsift_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(VARIANTS_BCFTOOLS.out.pangolin_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(VARIANTS_BCFTOOLS.out.nextclade_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(VARIANTS_BCFTOOLS.out.asciigenome_version.first().ifEmpty(null))

        //
        // MODULE: Get Pangolin lineage information for MultiQC report
        //
        ch_bcftools_pangolin_report
            .map { meta, report ->
                def fields = WorkflowCommons.getPangolinFieldMap(report)
                return [sample:meta.id] << fields
            }
            .set { ch_bcftools_pangolin_multiqc }

        MULTIQC_CUSTOM_CSV_BCFTOOLS_PANGOLIN (
            ch_bcftools_pangolin_multiqc.collect(),
            'bcftools_pangolin_lineage'
        )
        .set { ch_bcftools_pangolin_multiqc }

        //
        // MODULE: Get Nextclade clade information for MultiQC report
        //
        ch_bcftools_nextclade_report
            .map { meta, csv ->
                def clade = WorkflowCommons.getNextcladeFieldMapFromCsv(csv)['clade']
                return [ "$meta.id\t$clade" ]
            }
            .set { ch_bcftools_nextclade_multiqc }

        MULTIQC_CUSTOM_TSV_BCFTOOLS_NEXTCLADE (
            ch_bcftools_nextclade_multiqc.collect(),
            'Sample\tclade',
            'bcftools_nextclade_clade'
        )
        .set { ch_bcftools_nextclade_multiqc }
    }

    //
    // MODULE: Intersect variants across callers
    //
    if (!params.skip_variants && callers.size() > 1) {
        BCFTOOLS_ISEC (
            ch_ivar_vcf
                .join(ch_ivar_tbi, by: [0])
                .join(ch_bcftools_vcf, by: [0])
                .join(ch_bcftools_tbi, by: [0])
        )
    }

    //
    // MODULE: Primer trimming with Cutadapt
    //
    ch_cutadapt_multiqc = Channel.empty()
    if (params.protocol == 'amplicon' && !params.skip_assembly && !params.skip_cutadapt) {
        CUTADAPT (
            ch_assembly_fastq,
            PREPARE_GENOME.out.primer_fasta
        )
        ch_assembly_fastq    = CUTADAPT.out.reads
        ch_cutadapt_multiqc  = CUTADAPT.out.log
        ch_software_versions = ch_software_versions.mix(CUTADAPT.out.version.first().ifEmpty(null))

        if (!params.skip_fastqc) {
            FASTQC (
                CUTADAPT.out.reads
            )
        }
    }

    //
    // SUBWORKFLOW: Run SPAdes assembly and downstream analysis
    //
    ch_spades_quast_multiqc = Channel.empty()
    if (!params.skip_assembly && 'spades' in assemblers) {
        ASSEMBLY_SPADES (
            ch_assembly_fastq,
            ch_spades_hmm,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gff,
            PREPARE_GENOME.out.blast_db,
            ch_blast_outfmt6_header
        )
        ch_spades_quast_multiqc = ASSEMBLY_SPADES.out.quast_tsv
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_SPADES.out.spades_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_SPADES.out.bandage_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_SPADES.out.blast_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_SPADES.out.quast_version.ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_SPADES.out.abacas_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_SPADES.out.plasmidid_version.first().ifEmpty(null))
    }

    //
    // SUBWORKFLOW: Run Unicycler assembly and downstream analysis
    //
    ch_unicycler_quast_multiqc = Channel.empty()
    if (!params.skip_assembly && 'unicycler' in assemblers) {
        ASSEMBLY_UNICYCLER (
            ch_assembly_fastq,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gff,
            PREPARE_GENOME.out.blast_db,
            ch_blast_outfmt6_header
        )
        ch_unicycler_quast_multiqc = ASSEMBLY_UNICYCLER.out.quast_tsv
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_UNICYCLER.out.unicycler_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_UNICYCLER.out.bandage_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_UNICYCLER.out.blast_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_UNICYCLER.out.quast_version.ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_UNICYCLER.out.abacas_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_UNICYCLER.out.plasmidid_version.first().ifEmpty(null))
    }

    //
    // SUBWORKFLOW: Run minia assembly and downstream analysis
    //
    ch_minia_quast_multiqc = Channel.empty()
    if (!params.skip_assembly && 'minia' in assemblers) {
        ASSEMBLY_MINIA (
            ch_assembly_fastq,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gff,
            PREPARE_GENOME.out.blast_db,
            ch_blast_outfmt6_header
        )
        ch_minia_quast_multiqc = ASSEMBLY_MINIA.out.quast_tsv
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_MINIA.out.minia_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_MINIA.out.blast_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_MINIA.out.quast_version.ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_MINIA.out.abacas_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ASSEMBLY_MINIA.out.plasmidid_version.first().ifEmpty(null))
    }

    //
    // MODULE: Pipeline reporting
    //
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

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowCommons.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            GET_SOFTWARE_VERSIONS.out.yaml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_fail_reads_multiqc.ifEmpty([]),
            ch_fail_mapping_multiqc.ifEmpty([]),
            ch_amplicon_heatmap_multiqc.ifEmpty([]),
            FASTQC_FASTP.out.fastqc_raw_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_FASTP.out.trim_json.collect{it[1]}.ifEmpty([]),
            ch_kraken2_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_flagstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_multiqc.collect{it[1]}.ifEmpty([]),
            ch_ivar_trim_flagstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_flagstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_mosdepth_multiqc.collect{it[1]}.ifEmpty([]),
            ch_ivar_counts_multiqc.collect{it[1]}.ifEmpty([]),
            ch_ivar_stats_multiqc.collect{it[1]}.ifEmpty([]),
            ch_ivar_snpeff_multiqc.collect{it[1]}.ifEmpty([]),
            ch_ivar_quast_multiqc.collect().ifEmpty([]),
            ch_ivar_pangolin_multiqc.collect().ifEmpty([]),
            ch_ivar_nextclade_multiqc.collect().ifEmpty([]),
            ch_bcftools_stats_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bcftools_snpeff_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bcftools_quast_multiqc.collect().ifEmpty([]),
            ch_bcftools_pangolin_multiqc.collect().ifEmpty([]),
            ch_bcftools_nextclade_multiqc.collect().ifEmpty([]),
            ch_cutadapt_multiqc.collect{it[1]}.ifEmpty([]),
            ch_spades_quast_multiqc.collect().ifEmpty([]),
            ch_unicycler_quast_multiqc.collect().ifEmpty([]),
            ch_minia_quast_multiqc.collect().ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report, fail_mapped_reads)
    }
    NfcoreTemplate.summary(workflow, params, log, fail_mapped_reads, pass_mapped_reads)
}

/*
========================================================================================
    THE END
========================================================================================
*/
