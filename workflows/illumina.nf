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
    params.bowtie2_index, params.amplicon_fasta, params.amplicon_bed,
    params.multiqc_config
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet file not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Genome fasta file not specified!' }
if (params.gff)   { ch_gff   = file(params.gff)   }

// Check optional parameters
def platformList = ['illumina']
if (!platformList.contains(params.platform)) {
    exit 1, "Invalid platform option: ${params.platform}. Valid options: ${platformList.join(', ')}"
}

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
def assemblerList = [ 'spades', 'unicycler', 'minia', 'none' ]
def assemblers = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []
if ((assemblerList + assemblers).unique().size() != assemblerList.size()) {
    exit 1, "Invalid assembler option: ${params.assemblers}. Valid options: ${assemblerList.join(', ')}"
}

def spadesModeList = [ 'rnaviral', 'corona', 'metaviral', 'meta', 'metaplasmid', 'plasmid', 'isolate', 'rna', 'bio' ]
if (!spadesModeList.contains(params.spades_mode)) {
    exit 1, "Invalid spades mode option: ${params.spades_mode}. Valid options: ${spadesModeList.join(', ')}"
}
if (params.spades_hmm) { ch_spades_hmm = file(params.spades_hmm) } else { ch_spades_hmm = ch_dummy_file }

if (params.protocol == 'amplicon' && !params.skip_assembly && !params.amplicon_fasta) {
    exit 1, "To perform de novo assembly in 'amplicon' mode please provide a valid amplicon fasta file!"
}
if (params.amplicon_fasta) { ch_amplicon_fasta = file(params.amplicon_fasta) }

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

def publish_genome_options = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options  = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]

def cat_fastq_options          = modules['cat_fastq']
if (!params.save_merged_fastq) { cat_fastq_options['publish_files'] = false }

def samtools_mpileup_options    = modules['samtools_mpileup']
if (params.save_mpileup) { samtools_mpileup_options.publish_files.put('mpileup','') }

def varscan_mpileup2cns_options = modules['varscan_mpileup2cns']
varscan_mpileup2cns_options.args += " --min-var-freq $params.min_allele_freq"
varscan_mpileup2cns_options.args += (params.protocol != 'amplicon' && params.varscan2_strand_filter) ? " --strand-filter 1" : " --strand-filter 0"

def varscan_bcftools_options   = modules['varscan_bcftools_filter']
varscan_bcftools_options.args += " -i 'FORMAT/AD / (FORMAT/AD + FORMAT/RD) >= $params.max_allele_freq'"

def ivar_variants_options = modules['ivar_variants']
ivar_variants_options.args += " -t $params.min_allele_freq"

def ivar_consensus_options = modules['ivar_consensus']
ivar_consensus_options.args += " -t $params.max_allele_freq"

def ivar_variants_to_vcf_highfreq_options   = modules['ivar_variants_to_vcf_highfreq']
ivar_variants_to_vcf_highfreq_options.args += " --allele_freq_thresh $params.max_allele_freq"

def kraken2_run_options = modules['kraken2_run']
if (params.save_kraken2_fastq) { kraken2_run_options.publish_files.put('fastq.gz','') }

def multiqc_options         = modules['multiqc']
multiqc_options.args       += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''
// if (params.skip_alignment)  { multiqc_options['publish_dir'] = '' }

include { CAT_FASTQ                  } from '../modules/local/cat_fastq'                  addParams( options: cat_fastq_options                   ) 
include { MULTIQC_CUSTOM_FAIL_MAPPED } from '../modules/local/multiqc_custom_fail_mapped' addParams( options: [publish_files: false]              )
include { PICARD_COLLECTWGSMETRICS   } from '../modules/local/picard_collectwgsmetrics'   addParams( options: modules['picard_collectwgsmetrics'] )
include { COLLAPSE_AMPLICONS         } from '../modules/local/collapse_amplicons'         addParams( options: modules['collapse_amplicons']       )
include { SAMTOOLS_MPILEUP           } from '../modules/local/samtools_mpileup'           addParams( options: samtools_mpileup_options            )
include { VARSCAN_MPILEUP2CNS        } from '../modules/local/varscan_mpileup2cns'        addParams( options: varscan_mpileup2cns_options         )
include { BCFTOOLS_FILTER            } from '../modules/local/bcftools_filter'            addParams( options: modules['varscan_bcftools_filter']  )
include { IVAR_VARIANTS              } from '../modules/local/ivar_variants'              addParams( options: ivar_variants_options               )
include { IVAR_CONSENSUS             } from '../modules/local/ivar_consensus'             addParams( options: ivar_consensus_options              )
include { BCFTOOLS_MPILEUP           } from '../modules/local/bcftools_mpileup'           addParams( options: modules['bcftools_mpileup']         ) 
include { BCFTOOLS_ISEC              } from '../modules/local/bcftools_isec'              addParams( options: modules['bcftools_isec']            ) 
include { CUTADAPT                   } from '../modules/local/cutadapt'                   addParams( options: modules['cutadapt']                 )
include { KRAKEN2_RUN                } from '../modules/local/kraken2_run'                addParams( options: kraken2_run_options                 ) 
include { GET_SOFTWARE_VERSIONS      } from '../modules/local/get_software_versions'      addParams( options: [publish_files : ['csv':'']]        )
include { MULTIQC                    } from '../modules/local/multiqc'                    addParams( options: multiqc_options                     )

include { MOSDEPTH as MOSDEPTH_GENOME                             } from '../modules/local/mosdepth'              addParams( options: modules['mosdepth_genome']                )
include { MOSDEPTH as MOSDEPTH_AMPLICON                           } from '../modules/local/mosdepth'              addParams( options: modules['mosdepth_amplicon']              )
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_GENOME   } from '../modules/local/plot_mosdepth_regions' addParams( options: modules['plot_mosdepth_regions_genome']   )
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_AMPLICON } from '../modules/local/plot_mosdepth_regions' addParams( options: modules['plot_mosdepth_regions_amplicon'] )
include { IVAR_VARIANTS_TO_VCF as IVAR_VARIANTS_TO_VCF_LOWFREQ    } from '../modules/local/ivar_variants_to_vcf'  addParams( options: modules['ivar_variants_to_vcf_lowfreq']   )
include { IVAR_VARIANTS_TO_VCF as IVAR_VARIANTS_TO_VCF_HIGHFREQ   } from '../modules/local/ivar_variants_to_vcf'  addParams( options: ivar_variants_to_vcf_highfreq_options     )
include { QUAST as QUAST_VARSCAN                                  } from '../modules/local/quast'                 addParams( options: modules['varscan_quast']                  )
include { QUAST as QUAST_IVAR                                     } from '../modules/local/quast'                 addParams( options: modules['ivar_quast']                     )
include { QUAST as QUAST_BCFTOOLS                                 } from '../modules/local/quast'                 addParams( options: modules['bcftools_quast']                 )

/*
 * SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
 */
def fastp_options    = modules['fastp']
if (params.save_trimmed)      { fastp_options.publish_files.put('trim.fastq.gz','') }
if (params.save_trimmed_fail) { fastp_options.publish_files.put('fail.fastq.gz','') }

def bowtie2_build_options    = modules['bowtie2_build']
if (!params.save_reference) { bowtie2_build_options['publish_files'] = false }

def bowtie2_align_options         = modules['bowtie2_align']
if (params.save_align_intermeds) { bowtie2_align_options.publish_files.put('bam','') }
if (params.save_unaligned)       { bowtie2_align_options.publish_files.put('fastq.gz','unmapped') }
bowtie2_align_options.args2      += params.filter_unmapped ? "-F4" : ""
def bowtie2_sort_bam_options = modules['bowtie2_sort_bam']
if (params.save_align_intermeds || params.skip_markduplicates) {
    bowtie2_sort_bam_options.publish_files.put('bam','')
    bowtie2_sort_bam_options.publish_files.put('bai','')
}

def filter_bam_sort_bam_options = modules['filter_bam_sort_bam']

def ivar_trim_options   = modules['ivar_trim']
ivar_trim_options.args += params.ivar_trim_noprimer ? "" : " -e"

def ivar_trim_sort_bam_options = modules['ivar_trim_sort_bam']

def spades_options   = modules['spades']
spades_options.args += (params.spades_mode && params.spades_mode != 'corona') ? "--${params.spades_mode}" : ""

include { FASTQC_FASTP          } from '../subworkflows/local/fastqc_fastp'        addParams( fastqc_raw_options: modules['fastqc_raw'], fastqc_trim_options: modules['fastqc_trim'], fastp_options: fastp_options )
include { INPUT_CHECK           } from '../subworkflows/local/input_check'         addParams( options: [:] )
include { PREPARE_GENOME        } from '../subworkflows/local/prepare_genome'      addParams( genome_options: publish_genome_options, index_options: publish_index_options, bowtie2_index_options: bowtie2_build_options, makeblastdb_options: modules['blast_makeblastdb'], kraken2_build_options: modules['kraken2_build'])
include { ALIGN_BOWTIE2         } from '../subworkflows/local/align_bowtie2'       addParams( align_options: bowtie2_align_options, samtools_options: bowtie2_sort_bam_options )
include { FILTER_BAM_SAMTOOLS   } from '../subworkflows/local/filter_bam_samtools' addParams( samtools_view_options: modules['filter_bam'], samtools_index_options: filter_bam_sort_bam_options )
include { AMPLICON_TRIM_IVAR    } from '../subworkflows/local/amplicon_trim_ivar'  addParams( ivar_trim_options: ivar_trim_options, samtools_options: ivar_trim_sort_bam_options )
include { VCF_TABIX_STATS       } from '../subworkflows/local/vcf_tabix_stats'     addParams( tabix_options: modules['varscan_bcftools_filter_tabix'], stats_options: modules['varscan_bcftools_filter_stats'] )
include { ASSEMBLY_SPADES       } from '../subworkflows/local/assembly_spades'     addParams( spades_options: modules['spades'], bandage_options: modules['spades_bandage'], blastn_options: modules['spades_blastn'], abacas_options: modules['spades_abacas'], plasmidid_options: modules['spades_plasmidid'], quast_options: modules['spades_quast'], snpeff_options: modules['spades_snpeff'], snpeff_bgzip_options: modules['spades_snpeff_bgzip'], snpeff_tabix_options: modules['spades_snpeff_tabix'], snpeff_stats_options: modules['spades_snpeff_tabix'], snpsift_options: modules['spades_snpsift'] )
include { ASSEMBLY_UNICYCLER    } from '../subworkflows/local/assembly_unicycler'  addParams( unicycler_options: modules['unicycler'], bandage_options: modules['unicycler_bandage'], blastn_options: modules['unicycler_blastn'], abacas_options: modules['unicycler_abacas'], plasmidid_options: modules['unicycler_plasmidid'], quast_options: modules['unicycler_quast'], snpeff_options: modules['unicycler_snpeff'], snpeff_bgzip_options: modules['unicycler_snpeff_bgzip'], snpeff_tabix_options: modules['unicycler_snpeff_tabix'], snpeff_stats_options: modules['unicycler_snpeff_tabix'], snpsift_options: modules['unicycler_snpsift'] )
include { ASSEMBLY_MINIA        } from '../subworkflows/local/assembly_minia'      addParams( minia_options: modules['minia'], blastn_options: modules['minia_blastn'], abacas_options: modules['minia_abacas'], plasmidid_options: modules['minia_plasmidid'], quast_options: modules['minia_quast'], snpeff_options: modules['minia_snpeff'], snpeff_bgzip_options: modules['minia_snpeff_bgzip'], snpeff_tabix_options: modules['minia_snpeff_tabix'], snpeff_stats_options: modules['minia_snpeff_tabix'], snpsift_options: modules['minia_snpsift'] )

include { MAKE_CONSENSUS as MAKE_CONSENSUS_VARSCAN                     } from '../subworkflows/local/make_consensus'        addParams( genomecov_options: modules['varscan_consensus_genomecov'], merge_options: modules['varscan_consensus_merge'], mask_options: modules['varscan_consensus_mask'], maskfasta_options: modules['varscan_consensus_maskfasta'], bcftools_options: modules['varscan_consensus_bcftools'], plot_bases_options: modules['varscan_consensus_plot'] )
include { MAKE_CONSENSUS as MAKE_CONSENSUS_BCFTOOLS                    } from '../subworkflows/local/make_consensus'        addParams( genomecov_options: modules['bcftools_consensus_genomecov'], merge_options: modules['bcftools_consensus_merge'], mask_options: modules['bcftools_consensus_mask'], maskfasta_options: modules['bcftools_consensus_maskfasta'], bcftools_options: modules['bcftools_consensus_bcftools'], plot_bases_options: modules['bcftools_consensus_plot'] )
include { VCF_BGZIP_TABIX_STATS as VCF_BGZIP_TABIX_STATS_VARSCAN       } from '../subworkflows/local/vcf_bgzip_tabix_stats' addParams( bgzip_options: modules['varscan_bgzip'], tabix_options: modules['varscan_tabix'], stats_options: modules['varscan_stats'] )
include { VCF_BGZIP_TABIX_STATS as VCF_BGZIP_TABIX_STATS_IVAR_LOWFREQ  } from '../subworkflows/local/vcf_bgzip_tabix_stats' addParams( bgzip_options: modules['ivar_bgzip_lowfreq'], tabix_options: modules['ivar_tabix_lowfreq'], stats_options: modules['ivar_stats_lowfreq'] )
include { VCF_BGZIP_TABIX_STATS as VCF_BGZIP_TABIX_STATS_IVAR_HIGHFREQ } from '../subworkflows/local/vcf_bgzip_tabix_stats' addParams( bgzip_options: modules['ivar_bgzip_highfreq'], tabix_options: modules['ivar_tabix_highfreq'], stats_options: modules['ivar_stats_highfreq'] )
include { SNPEFF_SNPSIFT as SNPEFF_SNPSIFT_VARSCAN_LOWFREQ             } from '../subworkflows/local/snpeff_snpsift'        addParams( snpeff_options: modules['varscan_snpeff_lowfreq'], snpsift_options: modules['varscan_snpsift_lowfreq'], bgzip_options: modules['varscan_snpeff_lowfreq_bgzip'], tabix_options: modules['varscan_snpeff_lowfreq_tabix'], stats_options: modules['varscan_snpeff_lowfreq_stats'] )
include { SNPEFF_SNPSIFT as SNPEFF_SNPSIFT_VARSCAN_HIGHFREQ            } from '../subworkflows/local/snpeff_snpsift'        addParams( snpeff_options: modules['varscan_snpeff_highfreq'], snpsift_options: modules['varscan_snpsift_highfreq'], bgzip_options: modules['varscan_snpeff_highfreq_bgzip'], tabix_options: modules['varscan_snpeff_highfreq_tabix'], stats_options: modules['varscan_snpeff_highfreq_stats'] )
include { SNPEFF_SNPSIFT as SNPEFF_SNPSIFT_IVAR_LOWFREQ                } from '../subworkflows/local/snpeff_snpsift'        addParams( snpeff_options: modules['ivar_snpeff_lowfreq'], snpsift_options: modules['ivar_snpsift_lowfreq'], bgzip_options: modules['ivar_snpeff_lowfreq_bgzip'], tabix_options: modules['ivar_snpeff_lowfreq_tabix'], stats_options: modules['ivar_snpeff_lowfreq_stats'] )
include { SNPEFF_SNPSIFT as SNPEFF_SNPSIFT_IVAR_HIGHFREQ               } from '../subworkflows/local/snpeff_snpsift'        addParams( snpeff_options: modules['ivar_snpeff_highfreq'], snpsift_options: modules['ivar_snpsift_highfreq'], bgzip_options: modules['ivar_snpeff_highfreq_bgzip'], tabix_options: modules['ivar_snpeff_highfreq_tabix'], stats_options: modules['ivar_snpeff_highfreq_stats'] )
include { SNPEFF_SNPSIFT as SNPEFF_SNPSIFT_BCFTOOLS                    } from '../subworkflows/local/snpeff_snpsift'        addParams( snpeff_options: modules['bcftools_snpeff'], snpsift_options: modules['bcftools_snpsift'], bgzip_options: modules['bcftools_snpeff_bgzip'], tabix_options: modules['bcftools_snpeff_tabix'], stats_options: modules['bcftools_snpeff_stats'] )

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

/*
 * MODULE: Installed directly from nf-core/modules
 */
include { PICARD_COLLECTMULTIPLEMETRICS } from '../modules/nf-core/software/picard/collectmultiplemetrics/main' addParams( options: modules['picard_collectmultiplemetrics'] )
include { FASTQC                        } from '../modules/nf-core/software/fastqc/main'                        addParams( options: modules['cutadapt_fastqc']               )

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */
def markduplicates_options   = modules['picard_markduplicates']
markduplicates_options.args += params.filter_duplicates ? " REMOVE_DUPLICATES=true" : ""
include { MARK_DUPLICATES_PICARD } from '../subworkflows/nf-core/mark_duplicates_picard' addParams( markduplicates_options: markduplicates_options, samtools_options: modules['picard_markduplicates_sort_bam'] )

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
     * SUBWORKFLOW: Read QC, extract UMI and trim adapters
     */
    FASTQC_FASTP (
        ch_cat_fastq,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_trim_fastq        = FASTQC_FASTP.out.reads
    ch_software_versions = ch_software_versions.mix(FASTQC_FASTP.out.fastqc_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_FASTP.out.fastp_version.first().ifEmpty(null))
    
    /*
     * SUBWORKFLOW: Alignment with Bowtie2
     */
    ch_sorted_bam        = Channel.empty()
    ch_sorted_bai        = Channel.empty()
    ch_samtools_stats    = Channel.empty()
    ch_samtools_flagstat = Channel.empty()
    ch_samtools_idxstats = Channel.empty()
    ch_bowtie2_multiqc   = Channel.empty()
    if (!params.skip_variants) {
        ALIGN_BOWTIE2 (
            ch_trim_fastq,
            PREPARE_GENOME.out.bowtie2_index
        )
        ch_sorted_bam        = ALIGN_BOWTIE2.out.bam
        ch_sorted_bai        = ALIGN_BOWTIE2.out.bai
        ch_samtools_stats    = ALIGN_BOWTIE2.out.stats
        ch_samtools_flagstat = ALIGN_BOWTIE2.out.flagstat
        ch_samtools_idxstats = ALIGN_BOWTIE2.out.idxstats
        ch_bowtie2_multiqc   = ALIGN_BOWTIE2.out.log_out
        ch_software_versions = ch_software_versions.mix(ALIGN_BOWTIE2.out.bowtie2_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ALIGN_BOWTIE2.out.samtools_version.first().ifEmpty(null))
    }

    /*
     * Filter channels to get samples that passed Bowtie2 minimum mapped reads threshold
     */
    ch_fail_mapping_multiqc = Channel.empty()
    if (!params.skip_variants) {
        ch_samtools_flagstat
            .map { meta, flagstat -> [ meta ] + Checks.get_flagstat_mapped_reads(workflow, params, log, flagstat) }
            .set { ch_mapped_reads }

        ch_sorted_bam
            .join(ch_mapped_reads, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_sorted_bam }

        ch_sorted_bai
            .join(ch_mapped_reads, by: [0])
            .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
            .set { ch_sorted_bai }

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
     * SUBWORKFLOW: Filter unmapped reads from BAM and trim reads with iVar
     */
    ch_ivar_trim_multiqc = Channel.empty()
    if (!params.skip_variants && params.protocol == 'amplicon') {
        FILTER_BAM_SAMTOOLS (
            ch_sorted_bam
        )
        ch_sorted_bam = FILTER_BAM_SAMTOOLS.out.bam
        ch_sorted_bai = FILTER_BAM_SAMTOOLS.out.bai

        AMPLICON_TRIM_IVAR (
            ch_sorted_bam.join(ch_sorted_bai, by: [0]),
            ch_amplicon_bed
        )
        ch_sorted_bam             = AMPLICON_TRIM_IVAR.out.bam
        ch_sorted_bai             = AMPLICON_TRIM_IVAR.out.bai
        // ch_samtools_stats         = AMPLICON_TRIM_IVAR.out.stats
        // ch_samtools_flagstat      = AMPLICON_TRIM_IVAR.out.flagstat
        // ch_samtools_idxstats      = AMPLICON_TRIM_IVAR.out.idxstats
        ch_ivar_trim_multiqc      = AMPLICON_TRIM_IVAR.out.log_out
        ch_software_versions      = ch_software_versions.mix(AMPLICON_TRIM_IVAR.out.ivar_version.first().ifEmpty(null))
    }

    /*
     * SUBWORKFLOW: Mark duplicate reads
     */
    ch_markduplicates_multiqc = Channel.empty()
    if (!params.skip_variants && !params.skip_markduplicates) {
        MARK_DUPLICATES_PICARD (
            ch_sorted_bam
        )
        ch_sorted_bam             = MARK_DUPLICATES_PICARD.out.bam
        ch_sorted_bai             = MARK_DUPLICATES_PICARD.out.bai
        // ch_samtools_stats         = MARK_DUPLICATES_PICARD.out.stats
        // ch_samtools_flagstat      = MARK_DUPLICATES_PICARD.out.flagstat
        // ch_samtools_idxstats      = MARK_DUPLICATES_PICARD.out.idxstats
        ch_markduplicates_multiqc = MARK_DUPLICATES_PICARD.out.metrics
        ch_software_versions      = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.picard_version.first().ifEmpty(null))
    }

    /*
     * MODULE: Picard metrics
     */
    ch_picard_collectmultiplemetrics_multiqc = Channel.empty()
    ch_picard_collectwgsmetrics_multiqc      = Channel.empty()
    if (!params.skip_variants && !params.skip_picard_metrics) {
        PICARD_COLLECTMULTIPLEMETRICS (
            ch_sorted_bam,
            PREPARE_GENOME.out.fasta
        )
        ch_picard_collectmultiplemetrics_multiqc = PICARD_COLLECTMULTIPLEMETRICS.out.metrics

        PICARD_COLLECTWGSMETRICS (
            ch_sorted_bam,
            PREPARE_GENOME.out.fasta
        )
        ch_picard_collectwgsmetrics_multiqc = PICARD_COLLECTWGSMETRICS.out.metrics
    }

    /*
     * MODULE: Coverage QC plots
     */
    ch_mosdepth_multiqc = Channel.empty()
    if (!params.skip_variants && !params.skip_mosdepth) {

        MOSDEPTH_GENOME (
            ch_sorted_bam.join(ch_sorted_bai, by: [0]),
            ch_dummy_file,
            200
        )
        ch_mosdepth_multiqc = MOSDEPTH_GENOME.out.global_txt
        
        PLOT_MOSDEPTH_REGIONS_GENOME (
            MOSDEPTH_GENOME.out.regions_bed.collect { it[1] }
        )
        
        if (params.protocol == 'amplicon') {

            COLLAPSE_AMPLICONS (
                ch_amplicon_bed
            )

            MOSDEPTH_AMPLICON (
                ch_sorted_bam.join(ch_sorted_bai, by: [0]),
                COLLAPSE_AMPLICONS.out.bed,
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
            ch_sorted_bam,
            PREPARE_GENOME.out.fasta
        )
    }

    /*
     * SUBWORKFLOW: Call variants with VarScan2
     */
    if (!params.skip_variants && 'varscan2' in callers) {

        VARSCAN_MPILEUP2CNS (
            SAMTOOLS_MPILEUP.out.mpileup,
            PREPARE_GENOME.out.fasta
        )

        VCF_BGZIP_TABIX_STATS_VARSCAN (
            VARSCAN_MPILEUP2CNS.out.vcf
        )

        BCFTOOLS_FILTER (
            VCF_BGZIP_TABIX_STATS_VARSCAN.out.vcf
        )

        VCF_TABIX_STATS (
            BCFTOOLS_FILTER.out.vcf
        )

        if (!params.skip_consensus) {
            MAKE_CONSENSUS_VARSCAN (
                ch_sorted_bam
                    .join(BCFTOOLS_FILTER.out.vcf, by: [0])
                    .join(VCF_TABIX_STATS.out.tbi, by: [0]),
                PREPARE_GENOME.out.fasta
            )

            if (!params.skip_variants_quast) {
                QUAST_VARSCAN (
                    MAKE_CONSENSUS_VARSCAN.out.fasta.collect{ it[1] },
                    PREPARE_GENOME.out.fasta,
                    PREPARE_GENOME.out.gff
                )
            }
        }

        if (params.gff && !params.skip_variants_snpeff) {
            SNPEFF_SNPSIFT_VARSCAN_LOWFREQ (
                VCF_BGZIP_TABIX_STATS_VARSCAN.out.vcf,
                PREPARE_GENOME.out.snpeff_db,
                PREPARE_GENOME.out.snpeff_config,
                PREPARE_GENOME.out.fasta,
            )

            SNPEFF_SNPSIFT_VARSCAN_HIGHFREQ (
                BCFTOOLS_FILTER.out.vcf,
                PREPARE_GENOME.out.snpeff_db,
                PREPARE_GENOME.out.snpeff_config,
                PREPARE_GENOME.out.fasta
            )
        }
    }

    /*
     * SUBWORKFLOW: Call variants with IVar
     */
    if (!params.skip_variants && 'ivar' in callers) {

        IVAR_VARIANTS (
            SAMTOOLS_MPILEUP.out.mpileup,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gff
        )

        IVAR_VARIANTS_TO_VCF_LOWFREQ (
            IVAR_VARIANTS.out.tsv,
            ch_ivar_variants_header_mqc
        )

        VCF_BGZIP_TABIX_STATS_IVAR_LOWFREQ (
            IVAR_VARIANTS_TO_VCF_LOWFREQ.out.vcf
        )

        IVAR_VARIANTS_TO_VCF_HIGHFREQ (
            IVAR_VARIANTS.out.tsv,
            ch_ivar_variants_header_mqc
        )

        VCF_BGZIP_TABIX_STATS_IVAR_HIGHFREQ (
            IVAR_VARIANTS_TO_VCF_HIGHFREQ.out.vcf
        )

        if (!params.skip_consensus) {
            IVAR_CONSENSUS (
                SAMTOOLS_MPILEUP.out.mpileup
            )

            if (!params.skip_variants_quast) {
                QUAST_IVAR (
                    IVAR_CONSENSUS.out.fasta.collect{ it[1] },
                    PREPARE_GENOME.out.fasta,
                    PREPARE_GENOME.out.gff
                )
            }
        }

        if (params.gff && !params.skip_variants_snpeff) {
            SNPEFF_SNPSIFT_IVAR_LOWFREQ (
                VCF_BGZIP_TABIX_STATS_IVAR_LOWFREQ.out.vcf,
                PREPARE_GENOME.out.snpeff_db,
                PREPARE_GENOME.out.snpeff_config,
                PREPARE_GENOME.out.fasta,
            )

            SNPEFF_SNPSIFT_IVAR_HIGHFREQ (
                VCF_BGZIP_TABIX_STATS_IVAR_HIGHFREQ.out.vcf,
                PREPARE_GENOME.out.snpeff_db,
                PREPARE_GENOME.out.snpeff_config,
                PREPARE_GENOME.out.fasta
            )
        }
    }

    /*
     * SUBWORKFLOW: Call variants with BCFTools
     */
    if (!params.skip_variants && 'bcftools' in callers) {
        BCFTOOLS_MPILEUP (
            ch_sorted_bam,
            PREPARE_GENOME.out.fasta
        )

        if (!params.skip_consensus) {
            MAKE_CONSENSUS_BCFTOOLS (
                ch_sorted_bam
                    .join(BCFTOOLS_MPILEUP.out.vcf, by: [0])
                    .join(BCFTOOLS_MPILEUP.out.tbi, by: [0]),
                PREPARE_GENOME.out.fasta
            )

            if (!params.skip_variants_quast) {
                QUAST_BCFTOOLS (
                    MAKE_CONSENSUS_BCFTOOLS.out.fasta.collect{ it[1] },
                    PREPARE_GENOME.out.fasta,
                    PREPARE_GENOME.out.gff
                )
            }
        }

        if (params.gff && !params.skip_variants_snpeff) {
            SNPEFF_SNPSIFT_BCFTOOLS (
                BCFTOOLS_MPILEUP.out.vcf,
                PREPARE_GENOME.out.snpeff_db,
                PREPARE_GENOME.out.snpeff_config,
                PREPARE_GENOME.out.fasta,
            )
        }
    }

    /*
     * SUBWORKFLOW: Intersect variants across callers
     */
    if (!params.skip_variants && callers.size() > 2) {
        BCFTOOLS_ISEC (
            BCFTOOLS_FILTER.out.vcf
                .join(VCF_TABIX_STATS.out.tbi, by: [0])
                .join(VCF_BGZIP_TABIX_STATS_IVAR_HIGHFREQ.out.vcf, by: [0])
                .join(VCF_BGZIP_TABIX_STATS_IVAR_HIGHFREQ.out.tbi, by: [0])
                .join(BCFTOOLS_MPILEUP.out.vcf, by: [0])
                .join(BCFTOOLS_MPILEUP.out.tbi, by: [0])
        )
    }

    /*
     * MODULE: Amplicon trimming with Cutadapt
     */
    ch_cutadapt_multiqc        = Channel.empty()
    ch_cutadapt_fastqc_multiqc = Channel.empty()
    ch_cutadapt_version        = Channel.empty()
    if (params.protocol == 'amplicon' && !params.skip_assembly && !params.skip_amplicon_trimming) {
        CUTADAPT (
            ch_trim_fastq,
            PREPARE_GENOME.out.amplicon_fasta
        )
        ch_trim_fastq       = CUTADAPT.out.reads
        ch_cutadapt_multiqc = CUTADAPT.out.log
        ch_cutadapt_version = CUTADAPT.out.version

        if (!params.skip_fastqc) {
            FASTQC ( 
                CUTADAPT.out.reads 
            )
            ch_cutadapt_fastqc_multiqc = FASTQC.out.zip
        }
    }

    /*
     * MODULE: Run Kraken2
     */
    ch_kraken2_multiqc = Channel.empty()
    ch_kraken2_version = Channel.empty()
    if (!params.skip_assembly && !params.skip_kraken2) {
        KRAKEN2_RUN ( 
            ch_trim_fastq,
            PREPARE_GENOME.out.kraken2_db
        )
        ch_trim_fastq      = KRAKEN2_RUN.out.unclassified
        ch_kraken2_multiqc = KRAKEN2_RUN.out.txt
        ch_kraken2_version = KRAKEN2_RUN.out.version
    }

    /*
     * SUBWORKFLOW: Run SPAdes assembly and downstream analysis
     */
    ch_spades_version = Channel.empty()
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
        ch_spades_version = ASSEMBLY_SPADES.out.spades_version
    }

    /*
     * SUBWORKFLOW: Run Unicycler assembly and downstream analysis
     */
    ch_unicycler_version = Channel.empty()
    if (!params.skip_assembly && 'unicycler' in assemblers) {
        ASSEMBLY_UNICYCLER (
            ch_trim_fastq,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gff,
            PREPARE_GENOME.out.blast_db,
            PREPARE_GENOME.out.snpeff_db,
            PREPARE_GENOME.out.snpeff_config
        )
        ch_unicycler_version = ASSEMBLY_UNICYCLER.out.unicycler_version
    }

    /*
     * SUBWORKFLOW: Run minia assembly and downstream analysis
     */
    ch_minia_version = Channel.empty()
    if (!params.skip_assembly && 'minia' in assemblers) {
        ASSEMBLY_MINIA (
            ch_trim_fastq,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gff,
            PREPARE_GENOME.out.blast_db,
            PREPARE_GENOME.out.snpeff_db,
            PREPARE_GENOME.out.snpeff_config
        )
        ch_minia_version = ASSEMBLY_MINIA.out.minia_version
    }

    /*
     * MODULE: Pipeline reporting
     */
    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions.map { it }.collect()
    )

    /*
     * MultiQC
     */
    if (!params.skip_multiqc) {
        workflow_summary    = Schema.params_summary_multiqc(workflow, params.summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

    //     MULTIQC (
    //         ch_multiqc_config,
    //         ch_multiqc_custom_config.collect().ifEmpty([]),
    //         GET_SOFTWARE_VERSIONS.out.yaml.collect(),
    //         ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
    //         ch_fail_mapping_multiqc.ifEmpty([]),
    //         FASTQC_FASTP.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
    //         FASTQC_FASTP.out.trim_zip.collect{it[1]}.ifEmpty([]),
    //         FASTQC_FASTP.out.trim_log.collect{it[1]}.ifEmpty([]),
    //         ch_bowtie2_multiqc.collect{it[1]}.ifEmpty([]),
    //         ch_samtools_stats.collect{it[1]}.ifEmpty([]),
    //         ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
    //         ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
    //         ch_markduplicates_multiqc.collect{it[1]}.ifEmpty([]),
    //     )
    //     multiqc_report = MULTIQC.out.report.toList()
    }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

// workflow.onComplete {
//     Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report, fail_mapped_reads)
//     Completion.summary(workflow, params, log, fail_mapped_reads, pass_mapped_reads)
// }

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////

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

// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////
// /* --                                                                     -- */
// /* --                          MULTIQC                                    -- */
// /* --                                                                     -- */
// ///////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////

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

