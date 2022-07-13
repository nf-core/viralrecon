/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def valid_params = [
    protocols         : ['metagenomic', 'amplicon'],
    variant_callers   : ['ivar', 'bcftools'],
    consensus_callers : ['ivar', 'bcftools'],
    assemblers        : ['spades', 'unicycler', 'minia'],
    spades_modes      : ['rnaviral', 'corona', 'metaviral', 'meta', 'metaplasmid', 'plasmid', 'isolate', 'rna', 'bio']
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

if (params.input)      { ch_input      = file(params.input)      } else { exit 1, 'Input samplesheet file not specified!' }
if (params.spades_hmm) { ch_spades_hmm = file(params.spades_hmm) } else { ch_spades_hmm = []                              }

def assemblers = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []

def variant_caller = params.variant_caller
if (!variant_caller) { variant_caller = params.protocol == 'amplicon' ? 'ivar' : 'bcftools' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config_illumina.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? file(params.multiqc_config) : []

// Header files
ch_blast_outfmt6_header     = file("$projectDir/assets/headers/blast_outfmt6_header.txt", checkIfExists: true)
ch_ivar_variants_header_mqc = file("$projectDir/assets/headers/ivar_variants_header_mqc.txt", checkIfExists: true)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Loaded from modules/local/
//
include { CUTADAPT } from '../modules/local/cutadapt'
include { MULTIQC  } from '../modules/local/multiqc_illumina'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_GENOME   } from '../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_AMPLICON } from '../modules/local/plot_mosdepth_regions'
include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_FAIL_READS         } from '../modules/local/multiqc_tsv_from_list'
include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_FAIL_MAPPED        } from '../modules/local/multiqc_tsv_from_list'
include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_NEXTCLADE          } from '../modules/local/multiqc_tsv_from_list'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK         } from '../subworkflows/local/input_check'
include { PREPARE_GENOME      } from '../subworkflows/local/prepare_genome_illumina'
include { VARIANTS_IVAR       } from '../subworkflows/local/variants_ivar'
include { VARIANTS_BCFTOOLS   } from '../subworkflows/local/variants_bcftools'
include { CONSENSUS_IVAR      } from '../subworkflows/local/consensus_ivar'
include { CONSENSUS_BCFTOOLS  } from '../subworkflows/local/consensus_bcftools'
include { VARIANTS_LONG_TABLE } from '../subworkflows/local/variants_long_table'
include { ASSEMBLY_SPADES     } from '../subworkflows/local/assembly_spades'
include { ASSEMBLY_UNICYCLER  } from '../subworkflows/local/assembly_unicycler'
include { ASSEMBLY_MINIA      } from '../subworkflows/local/assembly_minia'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                     } from '../modules/nf-core/modules/cat/fastq/main'
include { FASTQC                        } from '../modules/nf-core/modules/fastqc/main'
include { KRAKEN2_KRAKEN2               } from '../modules/nf-core/modules/kraken2/kraken2/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from '../modules/nf-core/modules/picard/collectmultiplemetrics/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { MOSDEPTH as MOSDEPTH_GENOME   } from '../modules/nf-core/modules/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_AMPLICON } from '../modules/nf-core/modules/mosdepth/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { FASTQC_FASTP           } from '../subworkflows/nf-core/fastqc_fastp'
include { ALIGN_BOWTIE2          } from '../subworkflows/nf-core/align_bowtie2'
include { PRIMER_TRIM_IVAR       } from '../subworkflows/nf-core/primer_trim_ivar'
include { MARK_DUPLICATES_PICARD } from '../subworkflows/nf-core/mark_duplicates_picard'

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

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME ()
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
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input,
        params.platform
    )
    .sample_info
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
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_FASTP (
        ch_cat_fastq,
        params.save_trimmed_fail,
        false
    )
    ch_variants_fastq = FASTQC_FASTP.out.reads
    ch_versions = ch_versions.mix(FASTQC_FASTP.out.versions)

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

        MULTIQC_TSV_FAIL_READS (
            ch_pass_fail_reads.collect(),
            ['Sample', 'Reads before trimming'],
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
            PREPARE_GENOME.out.kraken2_db,
            params.kraken2_variants_host_filter || params.kraken2_assembly_host_filter,
            params.kraken2_variants_host_filter || params.kraken2_assembly_host_filter
        )
        ch_kraken2_multiqc = KRAKEN2_KRAKEN2.out.report
        ch_versions        = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first().ifEmpty(null))

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
        ALIGN_BOWTIE2 (
            ch_variants_fastq,
            PREPARE_GENOME.out.bowtie2_index,
            params.save_unaligned,
            false
        )
        ch_bam                      = ALIGN_BOWTIE2.out.bam
        ch_bai                      = ALIGN_BOWTIE2.out.bai
        ch_bowtie2_multiqc          = ALIGN_BOWTIE2.out.log_out
        ch_bowtie2_flagstat_multiqc = ALIGN_BOWTIE2.out.flagstat
        ch_versions                 = ch_versions.mix(ALIGN_BOWTIE2.out.versions)
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

        MULTIQC_TSV_FAIL_MAPPED (
            ch_pass_fail_mapped.fail.collect(),
            ['Sample', 'Mapped reads'],
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
        ch_versions                   = ch_versions.mix(PRIMER_TRIM_IVAR.out.versions)
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
        ch_versions                        = ch_versions.mix(MARK_DUPLICATES_PICARD.out.versions)
    }

    //
    // MODULE: Picard metrics
    //
    if (!params.skip_variants && !params.skip_picard_metrics) {
        PICARD_COLLECTMULTIPLEMETRICS (
            ch_bam,
            PREPARE_GENOME.out.fasta,
            []
        )
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first().ifEmpty(null))
    }

    //
    // MODULE: Genome-wide and amplicon-specific coverage QC plots
    //
    ch_mosdepth_multiqc         = Channel.empty()
    ch_amplicon_heatmap_multiqc = Channel.empty()
    if (!params.skip_variants && !params.skip_mosdepth) {
        MOSDEPTH_GENOME (
            ch_bam.join(ch_bai, by: [0]),
            [],
            []
        )
        ch_mosdepth_multiqc = MOSDEPTH_GENOME.out.global_txt
        ch_versions         = ch_versions.mix(MOSDEPTH_GENOME.out.versions.first().ifEmpty(null))

        PLOT_MOSDEPTH_REGIONS_GENOME (
            MOSDEPTH_GENOME.out.regions_bed.collect { it[1] }
        )
        ch_versions = ch_versions.mix(PLOT_MOSDEPTH_REGIONS_GENOME.out.versions)

        if (params.protocol == 'amplicon') {
            MOSDEPTH_AMPLICON (
                ch_bam.join(ch_bai, by: [0]),
                PREPARE_GENOME.out.primer_collapsed_bed,
                []
            )
            ch_versions = ch_versions.mix(MOSDEPTH_AMPLICON.out.versions.first().ifEmpty(null))

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
            (params.protocol == 'amplicon' || !params.skip_asciigenome) ? PREPARE_GENOME.out.fai : [],
            (params.protocol == 'amplicon' || !params.skip_asciigenome) ? PREPARE_GENOME.out.chrom_sizes : [],
            PREPARE_GENOME.out.gff,
            (params.protocol == 'amplicon' && params.primer_bed) ? PREPARE_GENOME.out.primer_bed : [],
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
            (params.protocol == 'amplicon' || !params.skip_asciigenome) ? PREPARE_GENOME.out.chrom_sizes : [],
            PREPARE_GENOME.out.gff,
            (params.protocol == 'amplicon' && params.primer_bed) ? PREPARE_GENOME.out.primer_bed : [],
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
    // SUBWORKFLOW: Call consensus with iVar and downstream QC
    //
    ch_quast_multiqc    = Channel.empty()
    ch_pangolin_multiqc = Channel.empty()
    ch_nextclade_report = Channel.empty()
    if (!params.skip_consensus && params.consensus_caller == 'ivar') {
        CONSENSUS_IVAR (
            ch_bam,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gff,
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
            PREPARE_GENOME.out.gff,
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
            .set { ch_nextclade_multiqc }

        MULTIQC_TSV_NEXTCLADE (
            ch_nextclade_multiqc.collect(),
            ['Sample', 'clade'],
            'nextclade_clade'
        )
        .set { ch_nextclade_multiqc }
    }

    //
    // SUBWORKFLOW: Create variants long table report
    //
    if (!params.skip_variants && !params.skip_variants_long_table && params.gff && !params.skip_snpeff) {
        VARIANTS_LONG_TABLE (
            ch_vcf,
            ch_tbi,
            ch_snpsift_txt,
            ch_pangolin_multiqc
        )
        ch_versions = ch_versions.mix(VARIANTS_LONG_TABLE.out.versions)
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
        ch_assembly_fastq   = CUTADAPT.out.reads
        ch_cutadapt_multiqc = CUTADAPT.out.log
        ch_versions         = ch_versions.mix(CUTADAPT.out.versions.first().ifEmpty(null))

        if (!params.skip_fastqc) {
            FASTQC (
                CUTADAPT.out.reads
            )
            ch_versions = ch_versions.mix(FASTQC.out.versions.first().ifEmpty(null))
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
            PREPARE_GENOME.out.gff,
            PREPARE_GENOME.out.blast_db,
            ch_blast_outfmt6_header
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
            PREPARE_GENOME.out.gff,
            PREPARE_GENOME.out.blast_db,
            ch_blast_outfmt6_header
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
            PREPARE_GENOME.out.gff,
            PREPARE_GENOME.out.blast_db,
            ch_blast_outfmt6_header
        )
        ch_minia_quast_multiqc = ASSEMBLY_MINIA.out.quast_tsv
        ch_versions            = ch_versions.mix(ASSEMBLY_MINIA.out.versions)
    }

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowCommons.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config,
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
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
            ch_bcftools_stats_multiqc.collect{it[1]}.ifEmpty([]),
            ch_snpeff_multiqc.collect{it[1]}.ifEmpty([]),
            ch_quast_multiqc.collect().ifEmpty([]),
            ch_pangolin_multiqc.collect{it[1]}.ifEmpty([]),
            ch_nextclade_multiqc.collect().ifEmpty([]),
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
