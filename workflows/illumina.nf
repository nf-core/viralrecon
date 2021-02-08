////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check genome key exists if provided
Checks.genome_exists(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.fasta, params.gff, 
    params.bowtie2_index, params.primer_bed, params.primer_fasta,
    params.multiqc_config
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

// Check mandatory parameters
if (params.input)  { ch_input = file(params.input) } else { exit 1, 'Input samplesheet file not specified!' }
if (!params.fasta) { exit 1, 'Genome fasta file not specified!' }

def protocolList = ['metagenomic', 'amplicon']
if (!protocolList.contains(params.protocol)) {
    exit 1, "Invalid protocol option: ${params.protocol}. Valid options: ${protocolList.join(', ')}"
}

// Variant calling parameter validation
def callerList = ['varscan2', 'ivar', 'bcftools']
def callers = params.callers ? params.callers.split(',').collect{ it.trim().toLowerCase() } : []
if ((callerList + callers).unique().size() != callerList.size()) {
    exit 1, "Invalid variant calller option: ${params.callers}. Valid options: ${callerList.join(', ')}"
}

if (params.protocol == 'amplicon' && !params.skip_variants && !params.primer_bed) {
    exit 1, "To perform variant calling in 'amplicon' mode please provide a valid primer BED file!"
}

// Assembly parameter validation
def assemblerList = ['spades', 'unicycler', 'minia']
def assemblers = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []
if ((assemblerList + assemblers).unique().size() != assemblerList.size()) {
    exit 1, "Invalid assembler option: ${params.assemblers}. Valid options: ${assemblerList.join(', ')}"
}

def spadesModeList = ['rnaviral', 'corona', 'metaviral', 'meta', 'metaplasmid', 'plasmid', 'isolate', 'rna', 'bio']
if (!spadesModeList.contains(params.spades_mode)) {
    exit 1, "Invalid spades mode option: ${params.spades_mode}. Valid options: ${spadesModeList.join(', ')}"
}
if (params.spades_hmm) { ch_spades_hmm = file(params.spades_hmm) } else { ch_spades_hmm = ch_dummy_file }

if (params.enable_conda && assemblers.contains('spades')) {
    assemblers = Checks.ignore_spades(params, log)
}

// if (!params.skip_kraken2 && !params.kraken2_db) {
//     if (!params.kraken2_db_name) { exit 1, "Please specify a valid name to build Kraken2 database for host e.g. 'human'!" }

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

def cat_fastq_options = modules['cat_fastq']
if (!params.save_merged_fastq) { cat_fastq_options['publish_files'] = false }

def cutadapt_options = modules['cutadapt']
if (params.save_trimmed) { cutadapt_options.publish_files.put('fastq.gz','') }

def kraken2_run_options = modules['kraken2_run']
if (params.save_kraken2_fastq) { kraken2_run_options.publish_files.put('fastq.gz','') }

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

include { CAT_FASTQ                  } from '../modules/local/cat_fastq'                  addParams( options: cat_fastq_options                   ) 
include { MULTIQC_CUSTOM_FAIL_MAPPED } from '../modules/local/multiqc_custom_fail_mapped' addParams( options: [publish_files: false]              )
include { PICARD_COLLECTWGSMETRICS   } from '../modules/local/picard_collectwgsmetrics'   addParams( options: modules['picard_collectwgsmetrics'] )
include { BCFTOOLS_ISEC              } from '../modules/local/bcftools_isec'              addParams( options: modules['bcftools_isec']            ) 
include { CUTADAPT                   } from '../modules/local/cutadapt'                   addParams( options: cutadapt_options                    )
include { KRAKEN2_RUN                } from '../modules/local/kraken2_run'                addParams( options: kraken2_run_options                 ) 
include { GET_SOFTWARE_VERSIONS      } from '../modules/local/get_software_versions'      addParams( options: [publish_files: ['csv':'']]         )
include { MULTIQC                    } from '../modules/local/multiqc'                    addParams( options: multiqc_options                     )

include { MOSDEPTH as MOSDEPTH_GENOME                             } from '../modules/local/mosdepth'              addParams( options: modules['mosdepth_genome']                )
include { MOSDEPTH as MOSDEPTH_AMPLICON                           } from '../modules/local/mosdepth'              addParams( options: modules['mosdepth_amplicon']              )
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_GENOME   } from '../modules/local/plot_mosdepth_regions' addParams( options: modules['plot_mosdepth_regions_genome']   )
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_AMPLICON } from '../modules/local/plot_mosdepth_regions' addParams( options: modules['plot_mosdepth_regions_amplicon'] )

/*
 * SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
 */
def publish_genome_options    = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options     = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]
def publish_db_options        = params.save_reference ? [publish_dir: 'genome/db']    : [publish_files: false]
def bedtools_getfasta_options = modules['bedtools_getfasta']
def bowtie2_build_options     = modules['bowtie2_build']
def snpeff_build_options      = modules['snpeff_build']
def makeblastdb_options       = modules['blast_makeblastdb']
def kraken2_build_options     = modules['kraken2_build']
def collapse_primers_options  = modules['collapse_primers']
if (!params.save_reference) { 
    bedtools_getfasta_options['publish_files'] = false
    bowtie2_build_options['publish_files']     = false
    snpeff_build_options['publish_files']      = false
    makeblastdb_options['publish_files']       = false
    kraken2_build_options['publish_files']     = false
    collapse_primers_options['publish_files']  = false
}

def ivar_trim_options   = modules['ivar_trim']
ivar_trim_options.args += params.ivar_trim_noprimer ? "" : " -e"
if (params.save_align_intermeds) { ivar_trim_options.publish_files.put('bam','') }

def ivar_trim_sort_bam_options = modules['ivar_trim_sort_bam']
if (params.save_align_intermeds || params.skip_markduplicates) {
    ivar_trim_sort_bam_options.publish_files.put('bam','')
    ivar_trim_sort_bam_options.publish_files.put('bai','')
}

if (params.save_align_intermeds) { bowtie2_align_options.publish_files.put('bam','') }
if (params.save_unaligned)       { bowtie2_align_options.publish_files.put('fastq.gz','unmapped') }

def varscan_mpileup2cns_options   = modules['varscan_mpileup2cns']
varscan_mpileup2cns_options.args += " --min-var-freq $params.min_allele_freq"
varscan_mpileup2cns_options.args += (params.protocol != 'amplicon' && params.varscan2_strand_filter) ? " --strand-filter 1" : " --strand-filter 0"

def varscan_bcftools_options   = modules['varscan_bcftools_filter']
varscan_bcftools_options.args += " -i 'FORMAT/AD / (FORMAT/AD + FORMAT/RD) >= $params.max_allele_freq'"

def ivar_variants_options   = modules['ivar_variants']
ivar_variants_options.args += " -t $params.min_allele_freq"

def ivar_consensus_options   = modules['ivar_consensus']
ivar_consensus_options.args += " -t $params.max_allele_freq"

def ivar_variants_to_vcf_highfreq_options   = modules['ivar_variants_to_vcf_highfreq']
ivar_variants_to_vcf_highfreq_options.args += " --allele_freq_thresh $params.max_allele_freq"

def spades_options   = modules['spades']
spades_options.args += (params.spades_mode && params.spades_mode != 'corona') ? " --${params.spades_mode}" : ""

include { INPUT_CHECK        } from '../subworkflows/local/input_check'        addParams( options: [:] )
include { PREPARE_GENOME     } from '../subworkflows/local/prepare_genome'     addParams( genome_options: publish_genome_options, index_options: publish_index_options, db_options: publish_db_options, bowtie2_build_options: bowtie2_build_options, bedtools_getfasta_options: bedtools_getfasta_options, collapse_primers_options: collapse_primers_options, snpeff_build_options: snpeff_build_options, makeblastdb_options: makeblastdb_options, kraken2_build_options: kraken2_build_options )
include { PRIMER_TRIM_IVAR   } from '../subworkflows/local/primer_trim_ivar'   addParams( ivar_trim_options: ivar_trim_options, samtools_options: ivar_trim_sort_bam_options )
include { VARIANTS_VARSCAN   } from '../subworkflows/local/variants_varscan'   addParams( varscan_mpileup2cns_options: varscan_mpileup2cns_options, quast_options: modules['varscan_quast'], bcftools_filter_options: modules['varscan_bcftools_filter'], consensus_genomecov_options: modules['varscan_consensus_genomecov'], consensus_merge_options: modules['varscan_consensus_merge'], consensus_mask_options: modules['varscan_consensus_mask'], consensus_maskfasta_options: modules['varscan_consensus_maskfasta'], consensus_bcftools_options: modules['varscan_consensus_bcftools'], consensus_plot_options: modules['varscan_consensus_plot'], bgzip_options: modules['varscan_bgzip'], tabix_options: modules['varscan_tabix'], stats_options: modules['varscan_stats'], bcftools_filter_tabix_options: modules['varscan_bcftools_filter_tabix'], bcftools_filter_stats_options: modules['varscan_bcftools_filter_stats'], snpeff_lowfreq_options: modules['varscan_snpeff_lowfreq'], snpsift_lowfreq_options: modules['varscan_snpsift_lowfreq'], snpeff_lowfreq_bgzip_options: modules['varscan_snpeff_lowfreq_bgzip'], snpeff_lowfreq_tabix_options: modules['varscan_snpeff_lowfreq_tabix'], snpeff_lowfreq_stats_options: modules['varscan_snpeff_lowfreq_stats'], snpeff_highfreq_options: modules['varscan_snpeff_highfreq'], snpsift_highfreq_options: modules['varscan_snpsift_highfreq'], snpeff_highfreq_bgzip_options: modules['varscan_snpeff_highfreq_bgzip'], snpeff_highfreq_tabix_options: modules['varscan_snpeff_highfreq_tabix'], snpeff_highfreq_stats_options: modules['varscan_snpeff_highfreq_stats'] )
include { VARIANTS_IVAR      } from '../subworkflows/local/variants_ivar'      addParams( ivar_variants_options: ivar_variants_options, ivar_consensus_options: ivar_consensus_options, quast_options: modules['ivar_quast'], ivar_variants_to_vcf_lowfreq_options: modules['ivar_variants_to_vcf_lowfreq'], ivar_variants_to_vcf_highfreq_options: ivar_variants_to_vcf_highfreq_options, ivar_bgzip_lowfreq_options: modules['ivar_bgzip_lowfreq'], ivar_tabix_lowfreq_options: modules['ivar_tabix_lowfreq'], ivar_stats_lowfreq_options: modules['ivar_stats_lowfreq'], ivar_bgzip_highfreq_options: modules['ivar_bgzip_highfreq'], ivar_tabix_highfreq_options: modules['ivar_tabix_highfreq'], ivar_stats_highfreq_options: modules['ivar_stats_highfreq'], snpeff_lowfreq_options: modules['ivar_snpeff_lowfreq'], snpsift_lowfreq_options: modules['ivar_snpsift_lowfreq'], snpeff_lowfreq_bgzip_options: modules['ivar_snpeff_lowfreq_bgzip'], snpeff_lowfreq_tabix_options: modules['ivar_snpeff_lowfreq_tabix'], snpeff_lowfreq_stats_options: modules['ivar_snpeff_lowfreq_stats'], snpeff_highfreq_options: modules['ivar_snpeff_highfreq'], snpsift_highfreq_options: modules['ivar_snpsift_highfreq'], snpeff_highfreq_bgzip_options: modules['ivar_snpeff_highfreq_bgzip'], snpeff_highfreq_tabix_options: modules['ivar_snpeff_highfreq_tabix'], snpeff_highfreq_stats_options: modules['ivar_snpeff_highfreq_stats'] )
include { VARIANTS_BCFTOOLS  } from '../subworkflows/local/variants_bcftools'  addParams( bcftools_mpileup_options: modules['bcftools_mpileup'], quast_options: modules['bcftools_quast'], consensus_genomecov_options: modules['bcftools_consensus_genomecov'], consensus_merge_options: modules['bcftools_consensus_merge'], consensus_mask_options: modules['bcftools_consensus_mask'], consensus_maskfasta_options: modules['bcftools_consensus_maskfasta'], consensus_bcftools_options: modules['bcftools_consensus_bcftools'], consensus_plot_options: modules['bcftools_consensus_plot'], snpeff_options: modules['bcftools_snpeff'], snpsift_options: modules['bcftools_snpsift'], snpeff_bgzip_options: modules['bcftools_snpeff_bgzip'], snpeff_tabix_options: modules['bcftools_snpeff_tabix'], snpeff_stats_options: modules['bcftools_snpeff_stats'] )
include { ASSEMBLY_SPADES    } from '../subworkflows/local/assembly_spades'    addParams( spades_options: spades_options, bandage_options: modules['spades_bandage'], blastn_options: modules['spades_blastn'], abacas_options: modules['spades_abacas'], plasmidid_options: modules['spades_plasmidid'], quast_options: modules['spades_quast'], snpeff_options: modules['spades_snpeff'], snpeff_bgzip_options: modules['spades_snpeff_bgzip'], snpeff_tabix_options: modules['spades_snpeff_tabix'], snpeff_stats_options: modules['spades_snpeff_tabix'], snpsift_options: modules['spades_snpsift'] )
include { ASSEMBLY_UNICYCLER } from '../subworkflows/local/assembly_unicycler' addParams( unicycler_options: modules['unicycler'], bandage_options: modules['unicycler_bandage'], blastn_options: modules['unicycler_blastn'], abacas_options: modules['unicycler_abacas'], plasmidid_options: modules['unicycler_plasmidid'], quast_options: modules['unicycler_quast'], snpeff_options: modules['unicycler_snpeff'], snpeff_bgzip_options: modules['unicycler_snpeff_bgzip'], snpeff_tabix_options: modules['unicycler_snpeff_tabix'], snpeff_stats_options: modules['unicycler_snpeff_tabix'], snpsift_options: modules['unicycler_snpsift'] )
include { ASSEMBLY_MINIA     } from '../subworkflows/local/assembly_minia'     addParams( minia_options: modules['minia'], blastn_options: modules['minia_blastn'], abacas_options: modules['minia_abacas'], plasmidid_options: modules['minia_plasmidid'], quast_options: modules['minia_quast'], snpeff_options: modules['minia_snpeff'], snpeff_bgzip_options: modules['minia_snpeff_bgzip'], snpeff_tabix_options: modules['minia_snpeff_tabix'], snpeff_stats_options: modules['minia_snpeff_tabix'], snpsift_options: modules['minia_snpsift'] )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

/*
 * MODULE: Installed directly from nf-core/modules
 */
include { FASTQC                        } from '../modules/nf-core/software/fastqc/main'                        addParams( options: modules['cutadapt_fastqc']               )
include { PICARD_COLLECTMULTIPLEMETRICS } from '../modules/nf-core/software/picard/collectmultiplemetrics/main' addParams( options: modules['picard_collectmultiplemetrics'] )
include { SAMTOOLS_MPILEUP              } from '../modules/nf-core/software/samtools/mpileup/main'              addParams( options: modules['samtools_mpileup']              )

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */
def fastp_options = modules['fastp']
if (params.save_trimmed)      { fastp_options.publish_files.put('trim.fastq.gz','') }
if (params.save_trimmed_fail) { fastp_options.publish_files.put('fail.fastq.gz','') }

def bowtie2_align_options = modules['bowtie2_align']
if (params.save_align_intermeds) { bowtie2_align_options.publish_files.put('bam','') }
if (params.save_unaligned)       { bowtie2_align_options.publish_files.put('fastq.gz','unmapped') }

def bowtie2_sort_bam_options = modules['bowtie2_sort_bam']
if ((params.protocol != 'amplicon' && params.skip_markduplicates) || params.save_align_intermeds) {
    bowtie2_sort_bam_options.publish_files.put('bam','')
    bowtie2_sort_bam_options.publish_files.put('bai','')
}

def markduplicates_options   = modules['picard_markduplicates']
markduplicates_options.args += params.filter_duplicates ? " REMOVE_DUPLICATES=true" : ""

include { FASTQC_FASTP           } from '../subworkflows/nf-core/fastqc_fastp'           addParams( fastqc_raw_options: modules['fastqc_raw'], fastqc_trim_options: modules['fastqc_trim'], fastp_options: fastp_options )
include { ALIGN_BOWTIE2          } from '../subworkflows/nf-core/align_bowtie2'          addParams( align_options: bowtie2_align_options, samtools_options: bowtie2_sort_bam_options                                     )
include { MARK_DUPLICATES_PICARD } from '../subworkflows/nf-core/mark_duplicates_picard' addParams( markduplicates_options: markduplicates_options, samtools_options: modules['picard_markduplicates_sort_bam']          )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report    = []
def pass_mapped_reads = [:]
def fail_mapped_reads = [:]

workflow ILLUMINA {

    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Uncompress and prepare reference genome files
     */
    PREPARE_GENOME (
        ch_dummy_file
    )

    // Check genome fasta only contains a single contig
    Checks.is_multifasta(PREPARE_GENOME.out.fasta, log)
    
    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK ( 
        ch_input
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
    
    /*
     * MODULE: Concatenate FastQ files from same sample if required
     */
    CAT_FASTQ ( 
        ch_fastq.multiple
    )
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    
    /*
     * SUBWORKFLOW: Read QC and trim adapters
     */
    FASTQC_FASTP (
        ch_cat_fastq
    )
    ch_trim_fastq          = FASTQC_FASTP.out.reads
    ch_fastqc_raw_multiqc  = FASTQC_FASTP.out.fastqc_raw_zip
    ch_fastqc_trim_multiqc = FASTQC_FASTP.out.fastqc_trim_zip
    ch_fastp_multiqc       = FASTQC_FASTP.out.trim_json
    ch_software_versions   = ch_software_versions.mix(FASTQC_FASTP.out.fastqc_version.first().ifEmpty(null))
    ch_software_versions   = ch_software_versions.mix(FASTQC_FASTP.out.fastp_version.first().ifEmpty(null))
    
    /*
     * SUBWORKFLOW: Alignment with Bowtie2
     */
    ch_bam               = Channel.empty()
    ch_bai               = Channel.empty()
    ch_bowtie2_stats     = Channel.empty()
    ch_bowtie2_flagstat  = Channel.empty()
    ch_bowtie2_idxstats  = Channel.empty()
    ch_bowtie2_multiqc   = Channel.empty()
    if (!params.skip_variants) {
        ALIGN_BOWTIE2 (
            ch_trim_fastq,
            PREPARE_GENOME.out.bowtie2_index
        )
        ch_bam               = ALIGN_BOWTIE2.out.bam
        ch_bai               = ALIGN_BOWTIE2.out.bai
        ch_bowtie2_stats     = ALIGN_BOWTIE2.out.stats
        ch_bowtie2_flagstat  = ALIGN_BOWTIE2.out.flagstat
        ch_bowtie2_idxstats  = ALIGN_BOWTIE2.out.idxstats
        ch_bowtie2_multiqc   = ALIGN_BOWTIE2.out.log_out
        ch_software_versions = ch_software_versions.mix(ALIGN_BOWTIE2.out.bowtie2_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ALIGN_BOWTIE2.out.samtools_version.first().ifEmpty(null))
    }

    /*
     * Filter channels to get samples that passed Bowtie2 minimum mapped reads threshold
     */
    ch_fail_mapping_multiqc = Channel.empty()
    if (!params.skip_variants) {
        ch_bowtie2_flagstat
            .map { meta, flagstat -> [ meta ] + Checks.get_flagstat_mapped_reads(workflow, params, log, flagstat) }
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

        MULTIQC_CUSTOM_FAIL_MAPPED ( 
            ch_pass_fail_mapped.fail.collect()
        )
        .set { ch_fail_mapping_multiqc }
    }
    
    /*
     * SUBWORKFLOW: Trim primer sequences from reads with iVar
     */
    ch_ivar_trim_stats    = Channel.empty()
    ch_ivar_trim_flagstat = Channel.empty()
    ch_ivar_trim_idxstats = Channel.empty()
    ch_ivar_trim_multiqc  = Channel.empty()
    if (!params.skip_variants && params.protocol == 'amplicon') {
        PRIMER_TRIM_IVAR (
            ch_bam.join(ch_bai, by: [0]),
            PREPARE_GENOME.out.primer_bed
        )
        ch_bam                = PRIMER_TRIM_IVAR.out.bam
        ch_bai                = PRIMER_TRIM_IVAR.out.bai
        ch_ivar_trim_stats    = PRIMER_TRIM_IVAR.out.stats
        ch_ivar_trim_flagstat = PRIMER_TRIM_IVAR.out.flagstat
        ch_ivar_trim_idxstats = PRIMER_TRIM_IVAR.out.idxstats
        ch_ivar_trim_multiqc  = PRIMER_TRIM_IVAR.out.log_out
        ch_software_versions  = ch_software_versions.mix(PRIMER_TRIM_IVAR.out.ivar_version.first().ifEmpty(null))
    }

    /*
     * SUBWORKFLOW: Mark duplicate reads
     */
    ch_markduplicates_stats    = Channel.empty()
    ch_markduplicates_flagstat = Channel.empty()
    ch_markduplicates_idxstats = Channel.empty()
    ch_markduplicates_multiqc  = Channel.empty()
    if (!params.skip_variants && !params.skip_markduplicates) {
        MARK_DUPLICATES_PICARD (
            ch_bam
        )
        ch_bam                     = MARK_DUPLICATES_PICARD.out.bam
        ch_bai                     = MARK_DUPLICATES_PICARD.out.bai
        ch_markduplicates_stats    = MARK_DUPLICATES_PICARD.out.stats
        ch_markduplicates_flagstat = MARK_DUPLICATES_PICARD.out.flagstat
        ch_markduplicates_idxstats = MARK_DUPLICATES_PICARD.out.idxstats
        ch_markduplicates_multiqc  = MARK_DUPLICATES_PICARD.out.metrics
        ch_software_versions = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.picard_version.first().ifEmpty(null))
    }

    /*
     * MODULE: Picard metrics
     */
    ch_picard_collectwgsmetrics_multiqc      = Channel.empty()
    ch_picard_collectmultiplemetrics_multiqc = Channel.empty()
    if (!params.skip_variants && !params.skip_picard_metrics) {
        PICARD_COLLECTWGSMETRICS (
            ch_bam,
            PREPARE_GENOME.out.fasta
        )
        ch_picard_collectwgsmetrics_multiqc = PICARD_COLLECTWGSMETRICS.out.metrics
        ch_software_versions = ch_software_versions.mix(PICARD_COLLECTWGSMETRICS.out.version.first().ifEmpty(null))

        PICARD_COLLECTMULTIPLEMETRICS (
            ch_bam,
            PREPARE_GENOME.out.fasta
        )
        ch_picard_collectmultiplemetrics_multiqc = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
    }

    /*
     * MODULE: Genome-wide and amplicon-specific coverage QC plots
     */
    ch_mosdepth_multiqc = Channel.empty()
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
        }
    }

    /*
     * MODULE: Make mpileup to re-use across callers
     */
    if (!params.skip_variants) {
        SAMTOOLS_MPILEUP (
            ch_bam,
            PREPARE_GENOME.out.fasta
        )
    }

    // /*
    //  * SUBWORKFLOW: Call variants with VarScan2
    //  */
    // if (!params.skip_variants && 'varscan2' in callers) {
    //     VARIANTS_VARSCAN (
    //         SAMTOOLS_MPILEUP.out.mpileup,
    //         ch_bam,
    //         PREPARE_GENOME.out.fasta,
    //         PREPARE_GENOME.out.gff,
    //         PREPARE_GENOME.out.snpeff_db,
    //         PREPARE_GENOME.out.snpeff_config
    //     )
    // }

    // /*
    //  * SUBWORKFLOW: Call variants with IVar
    //  */
    // if (!params.skip_variants && 'ivar' in callers) {
    //     VARIANTS_IVAR (
    //         SAMTOOLS_MPILEUP.out.mpileup,
    //         PREPARE_GENOME.out.fasta,
    //         PREPARE_GENOME.out.gff,
    //         PREPARE_GENOME.out.snpeff_db,
    //         PREPARE_GENOME.out.snpeff_config,
    //         ch_ivar_variants_header_mqc
    //     )
    // }

    /*
     * SUBWORKFLOW: Call variants with BCFTools
     */
    if (!params.skip_variants && 'bcftools' in callers) {
        VARIANTS_BCFTOOLS (
            ch_bam,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gff,
            PREPARE_GENOME.out.snpeff_db,
            PREPARE_GENOME.out.snpeff_config
        )
    }

    // // /*
    // //  * MODULE: Intersect variants across callers
    // //  */
    // // if (!params.skip_variants && callers.size() > 2) {
    // //     BCFTOOLS_ISEC (
    // //         VARIANTS_VARSCAN.out.vcf
    // //             .join(VCF_TABIX_STATS.out.tbi, by: [0])
    // //             .join(VCF_BGZIP_TABIX_STATS_IVAR_HIGHFREQ.out.vcf, by: [0])
    // //             .join(VCF_BGZIP_TABIX_STATS_IVAR_HIGHFREQ.out.tbi, by: [0])
    // //             .join(BCFTOOLS_MPILEUP.out.vcf, by: [0])
    // //             .join(BCFTOOLS_MPILEUP.out.tbi, by: [0])
    // //     )
    // // }

    /*
     * MODULE: Primer trimming with Cutadapt
     */
    ch_cutadapt_multiqc        = Channel.empty()
    ch_cutadapt_fastqc_multiqc = Channel.empty()
    if (params.protocol == 'amplicon' && !params.skip_assembly && !params.skip_cutadapt) {
        CUTADAPT (
            ch_trim_fastq,
            PREPARE_GENOME.out.primer_fasta
        )
        ch_trim_fastq        = CUTADAPT.out.reads
        ch_cutadapt_multiqc  = CUTADAPT.out.log
        ch_software_versions = ch_software_versions.mix(CUTADAPT.out.version.first().ifEmpty(null))

        if (!params.skip_fastqc) {
            FASTQC ( 
                CUTADAPT.out.reads 
            )
            ch_cutadapt_fastqc_multiqc = FASTQC.out.zip
        }
    }

    /*
     * MODULE: Run Kraken2 for removal of host reads
     */
    ch_kraken2_multiqc = Channel.empty()
    if (!params.skip_assembly && !params.skip_kraken2) {
        KRAKEN2_RUN ( 
            ch_trim_fastq,
            PREPARE_GENOME.out.kraken2_db
        )
        ch_trim_fastq        = KRAKEN2_RUN.out.unclassified
        ch_kraken2_multiqc   = KRAKEN2_RUN.out.txt
        ch_software_versions = ch_software_versions.mix(KRAKEN2_RUN.out.version.first().ifEmpty(null))
    }

    /*
     * SUBWORKFLOW: Run SPAdes assembly and downstream analysis
     */
    ch_spades_quast_multiqc    = Channel.empty()
    // ch_spades_bcftools_multiqc = Channel.empty()
    // ch_spades_snpeff_multiqc   = Channel.empty()
    if (!params.skip_assembly && 'spades' in assemblers) {
        ASSEMBLY_SPADES (
            ch_trim_fastq,
            ch_spades_hmm,
            params.spades_mode == 'corona',
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gff,
            PREPARE_GENOME.out.blast_db,
            PREPARE_GENOME.out.snpeff_db,
            PREPARE_GENOME.out.snpeff_config
        )
        ch_spades_quast_multiqc    = ASSEMBLY_SPADES.out.quast_tsv
        // ch_spades_bcftools_multiqc = ASSEMBLY_SPADES.out.quast_tsv
        // ch_spades_snpeff_multiqc   = ASSEMBLY_SPADES.out.quast_tsv
        ch_software_versions       = ch_software_versions.mix(ASSEMBLY_SPADES.out.spades_version.first().ifEmpty(null))
        ch_software_versions       = ch_software_versions.mix(ASSEMBLY_SPADES.out.bandage_version.first().ifEmpty(null))
        ch_software_versions       = ch_software_versions.mix(ASSEMBLY_SPADES.out.blast_version.first().ifEmpty(null))
        ch_software_versions       = ch_software_versions.mix(ASSEMBLY_SPADES.out.quast_version.ifEmpty(null))
        ch_software_versions       = ch_software_versions.mix(ASSEMBLY_SPADES.out.abacas_version.first().ifEmpty(null))
        ch_software_versions       = ch_software_versions.mix(ASSEMBLY_SPADES.out.plasmidid_version.first().ifEmpty(null))
    }

    /*
     * SUBWORKFLOW: Run Unicycler assembly and downstream analysis
     */
    ch_unicycler_quast_multiqc    = Channel.empty()
    // ch_unicycler_bcftools_multiqc = Channel.empty()
    // ch_unicycler_snpeff_multiqc   = Channel.empty()
    if (!params.skip_assembly && 'unicycler' in assemblers) {
        ASSEMBLY_UNICYCLER (
            ch_trim_fastq,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gff,
            PREPARE_GENOME.out.blast_db,
            PREPARE_GENOME.out.snpeff_db,
            PREPARE_GENOME.out.snpeff_config
        )
        ch_unicycler_quast_multiqc    = ASSEMBLY_UNICYCLER.out.quast_tsv
        // ch_unicycler_bcftools_multiqc = ASSEMBLY_UNICYCLER.out.quast_tsv
        // ch_unicycler_snpeff_multiqc   = ASSEMBLY_UNICYCLER.out.quast_tsv
        ch_software_versions          = ch_software_versions.mix(ASSEMBLY_UNICYCLER.out.unicycler_version.first().ifEmpty(null))
        ch_software_versions          = ch_software_versions.mix(ASSEMBLY_UNICYCLER.out.bandage_version.first().ifEmpty(null))
        ch_software_versions          = ch_software_versions.mix(ASSEMBLY_UNICYCLER.out.blast_version.first().ifEmpty(null))
        ch_software_versions          = ch_software_versions.mix(ASSEMBLY_UNICYCLER.out.quast_version.ifEmpty(null))
        ch_software_versions          = ch_software_versions.mix(ASSEMBLY_UNICYCLER.out.abacas_version.first().ifEmpty(null))
        ch_software_versions          = ch_software_versions.mix(ASSEMBLY_UNICYCLER.out.plasmidid_version.first().ifEmpty(null))
    }

    /*
     * SUBWORKFLOW: Run minia assembly and downstream analysis
     */
    ch_minia_quast_multiqc    = Channel.empty()
    // ch_minia_bcftools_multiqc = Channel.empty()
    // ch_minia_snpeff_multiqc   = Channel.empty()
    if (!params.skip_assembly && 'minia' in assemblers) {
        ASSEMBLY_MINIA (
            ch_trim_fastq,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gff,
            PREPARE_GENOME.out.blast_db,
            PREPARE_GENOME.out.snpeff_db,
            PREPARE_GENOME.out.snpeff_config
        )
        ch_minia_quast_multiqc    = ASSEMBLY_MINIA.out.quast_tsv
        // ch_minia_bcftools_multiqc = ASSEMBLY_MINIA.out.quast_tsv
        // ch_minia_snpeff_multiqc   = ASSEMBLY_MINIA.out.quast_tsv
        ch_software_versions      = ch_software_versions.mix(ASSEMBLY_MINIA.out.minia_version.first().ifEmpty(null))
        ch_software_versions      = ch_software_versions.mix(ASSEMBLY_MINIA.out.blast_version.first().ifEmpty(null))
        ch_software_versions      = ch_software_versions.mix(ASSEMBLY_MINIA.out.quast_version.ifEmpty(null))
        ch_software_versions      = ch_software_versions.mix(ASSEMBLY_MINIA.out.abacas_version.first().ifEmpty(null))
        ch_software_versions      = ch_software_versions.mix(ASSEMBLY_MINIA.out.plasmidid_version.first().ifEmpty(null))
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

    // /*
    //  * MODULE: MultiQC
    //  */
    // if (!params.skip_multiqc) {
    //     workflow_summary    = Schema.params_summary_multiqc(workflow, params.summary_params)
    //     ch_workflow_summary = Channel.value(workflow_summary)

    // //     MULTIQC (
    // //         ch_multiqc_config,
    // //         ch_multiqc_custom_config.collect().ifEmpty([]),
    // //         GET_SOFTWARE_VERSIONS.out.yaml.collect(),
    // //         ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
    // //         ch_fail_mapping_multiqc.ifEmpty([]),
    // //         FASTQC_FASTP.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
    // //         FASTQC_FASTP.out.trim_zip.collect{it[1]}.ifEmpty([]),
    // //         FASTQC_FASTP.out.trim_log.collect{it[1]}.ifEmpty([]),
    // //         ch_bowtie2_multiqc.collect{it[1]}.ifEmpty([]),
    // //         ch_samtools_stats.collect{it[1]}.ifEmpty([]),
    // //         ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
    // //         ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
    // //         ch_markduplicates_multiqc.collect{it[1]}.ifEmpty([]),
    // //     )
    // //     multiqc_report = MULTIQC.out.report.toList()
    // }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report, fail_mapped_reads)
    Completion.summary(workflow, params, log, fail_mapped_reads, pass_mapped_reads)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////