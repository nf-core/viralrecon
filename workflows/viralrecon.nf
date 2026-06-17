/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog             } from 'plugin/nf-schema'
include { paramsSummaryMap             } from 'plugin/nf-schema'
include { paramsSummaryMultiqc         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText       } from '../subworkflows/local/utils_nfcore_viralrecon_pipeline'
include { getFlagstatMappedReads       } from '../subworkflows/local/utils_nfcore_viralrecon_pipeline'
include { multiqcTsvFromList           } from '../subworkflows/local/utils_nfcore_viralrecon_pipeline'
include { checkPrimerSuffixes          } from '../subworkflows/local/utils_nfcore_viralrecon_pipeline'
include { getColFromFile               } from '../subworkflows/local/utils_nfcore_viralrecon_pipeline'
include { checkContigsInBED            } from '../subworkflows/local/utils_nfcore_viralrecon_pipeline'
include { getNextcladeFieldMapFromCsv  } from '../subworkflows/local/utils_nfcore_viralrecon_pipeline'
include { isMultiFasta                 } from '../subworkflows/local/utils_nfcore_viralrecon_pipeline'
include { checkIfSwiftProtocol         } from '../subworkflows/local/utils_nfcore_viralrecon_pipeline'
include { getFastpReadsAfterFiltering  } from '../subworkflows/local/utils_nfcore_viralrecon_pipeline'
include { getFastpReadsBeforeFiltering } from '../subworkflows/local/utils_nfcore_viralrecon_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_GENOME   } from '../modules/local/plot_mosdepth_regions'
include { PLOT_MOSDEPTH_REGIONS as PLOT_MOSDEPTH_REGIONS_AMPLICON } from '../modules/local/plot_mosdepth_regions'
include { PREPARE_PRIMER_FASTA                                    } from '../modules/local/prepare_primer_fasta'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_GENOME_ILLUMINA } from '../subworkflows/local/prepare_genome_illumina'
include { PREPARE_GENOME_NANOPORE } from '../subworkflows/local/prepare_genome_nanopore'

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
include { SNPEFF_SNPSIFT          } from '../subworkflows/local/snpeff_snpsift'
include { FILTER_BAM_SAMTOOLS     } from '../subworkflows/local/filter_bam_samtools'
include { HIV_RESISTANCE          } from '../subworkflows/local/hiv_resitance_detection'
include { ARTIC_MINION_PROTOCOL   } from '../subworkflows/local/artic_minion_protocol'
include { MINIMAP2_MAPPING         } from '../subworkflows/local/minimap2_mapping'

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
include { PYCOQC                        } from '../modules/nf-core/pycoqc/main'
include { NANOPLOT                      } from '../modules/nf-core/nanoplot/main'
include { ARTIC_GUPPYPLEX               } from '../modules/nf-core/artic/guppyplex/main'
include { BCFTOOLS_STATS                } from '../modules/nf-core/bcftools/stats/main'
include { QUAST                         } from '../modules/nf-core/quast/main'
include { PANGOLIN_UPDATEDATA           } from '../modules/nf-core/pangolin/updatedata/main'
include { PANGOLIN_RUN                  } from '../modules/nf-core/pangolin/run/main'
include { NEXTCLADE_RUN                 } from '../modules/nf-core/nextclade/run/main'
include { MULTIQC                       } from '../modules/nf-core/multiqc/main'
include { UNTAR as UNTAR_PANGODB        } from '../modules/nf-core/untar/main'
include { GUNZIP as GUNZIP_GFF          } from '../modules/nf-core/gunzip/main'

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

workflow VIRALRECON {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir
    ch_genome_fasta
    ch_genome_gff
    ch_primer_bed
    ch_bowtie2_index
    ch_nextclade_dataset
    ch_nextclade_dataset_name
    ch_nextclade_dataset_tag
    ch_artic_scheme

    main:
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        VALIDATE INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def ch_multiqc_files = channel.empty()
    def ch_versions = channel.empty()


    def valid_params = [
        protocols            : ['metagenomic', 'amplicon'],
        variant_callers      : ['ivar', 'bcftools'],
        consensus_callers    : ['ivar', 'bcftools'],
        assemblers           : ['spades', 'unicycler', 'minia'],
        spades_modes         : ['rnaviral', 'corona', 'metaviral', 'meta', 'metaplasmid', 'plasmid', 'isolate', 'rna', 'bio'],
    ]

    def checkPathParamList = []
    def sequencing_summary = (params.sequencing_summary == false || params.sequencing_summary == 'false') ? null : params.sequencing_summary

    // Check input path parameters to see if they exist
    if (params.platform == 'illumina') {
        checkPathParamList = [
            params.input, params.fasta, ch_genome_gff, params.bowtie2_index,
            params.kraken2_db, params.primer_bed, params.primer_fasta,
            params.blast_db, params.spades_hmm, params.multiqc_config,
            params.freyja_barcodes, params.freyja_lineages_meta, params.freyja_lineages_topology, params.additional_annotation
        ]
    } else if (params.platform == 'nanopore') {
        checkPathParamList = [
            params.input, params.fastq_dir,
            sequencing_summary, ch_genome_gff,
            params.freyja_barcodes, params.freyja_lineages_meta, params.freyja_lineages_topology, params.additional_annotation,
            params.kraken2_db
        ]
    }

    checkPathParamList.each { param ->
        if (param) { file(param, checkIfExists: true) }
    }

    if (params.input)                 { ch_input          = file(params.input)                 } else { exit 1, 'Input samplesheet file not specified!' }
    if (params.spades_hmm)            { ch_spades_hmm     = file(params.spades_hmm)            } else { ch_spades_hmm = []                              }
    if (params.additional_annotation) { ch_additional_gtf = file(params.additional_annotation) } else { ch_additional_gtf = channel.empty()             }
    if (params.taxidlist)             { ch_taxidlist      = file(params.taxidlist)             } else { ch_taxidlist = []                               }

    // If protocol amplicon you must provide primer set information
    if (params.protocol == 'amplicon' && !params.skip_variants && !params.primer_bed) {
        error("To perform variant calling in amplicon mode please provide a valid primer BED file e.g. '--primer_bed primers.bed'.")
    }

    if (!params.fasta) {
        error("Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file.")
    }

    if (!params.skip_kraken2 && !params.kraken2_db) {
        if (!params.kraken2_db_name) {
            error("Please specify a valid name to build Kraken2 database for host e.g. '--kraken2_db_name human'.")
        }
    }

    def assemblers = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []

    def variant_caller = params.variant_caller
    if (!variant_caller) { variant_caller = params.trim_primers ? 'ivar' : 'bcftools' }

    if (sequencing_summary)             { ch_sequencing_summary = file(sequencing_summary)             } else { ch_sequencing_summary = [] }

    // Need to stage artic model properly depending on whether it is a string or a file
    ch_clair3_model = params.clair3_model_dir ?
        channel.value([ [:], file(params.clair3_model_dir, type: 'dir'), params.clair3_model ]) :
        channel.value([ [:], [], params.clair3_model ])

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CONFIG FILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Header files
    ch_blast_outfmt6_header          = file("$projectDir/assets/headers/blast_outfmt6_header.txt", checkIfExists: true)
    ch_blast_filtered_outfmt6_header = file("$projectDir/assets/headers/blast_filtered_outfmt6_header.txt", checkIfExists: true)
    ch_ivar_variants_header_mqc      = file("$projectDir/assets/headers/ivar_variants_header_mqc.txt", checkIfExists: true)

    // Info required for completion email and summary
    def pass_mapped_reads  = [:]
    def fail_mapped_reads  = [:]
    def pass_barcode_reads = [:]
    def fail_barcode_reads = [:]

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    if (params.platform == 'illumina') {
        genome = PREPARE_GENOME_ILLUMINA (
            ch_genome_fasta,
            ch_genome_gff,
            ch_primer_bed,
            ch_bowtie2_index,
            ch_nextclade_dataset,
            ch_nextclade_dataset_name,
            ch_nextclade_dataset_tag
        )
    } else if (params.platform == 'nanopore') {
        genome = PREPARE_GENOME_NANOPORE (
            ch_genome_fasta,
            ch_genome_gff,
            ch_primer_bed,
            ch_bowtie2_index,
            ch_nextclade_dataset,
            ch_nextclade_dataset_name,
            ch_nextclade_dataset_tag
        )
    }
    ch_versions = ch_versions.mix(genome.versions)

    if (params.platform == 'illumina') {
        //
        // ILLUMINA WORKFLOW
        //

        def ch_reference_fasta_fai = genome.fasta
            .combine(genome.fai)
            .map { fasta_file, fai_file -> [ [:], fasta_file, fai_file ] }
            .collect(flat: false)
            .map { fasta_fai -> fasta_fai[0] }

        // Check genome fasta only contains a single contig
        genome
            .fasta
            .map { isMultiFasta(it, log) }

        if (params.trim_primers && !params.skip_variants) {
            // Check primer BED file only contains suffixes provided --primer_left_suffix / --primer_right_suffix
            genome
                .primer_bed
                .map { checkPrimerSuffixes(it, params.primer_left_suffix, params.primer_right_suffix, log) }

            // Check whether the contigs in the primer BED file are present in the reference genome
            genome
                .primer_bed
                .map { [ getColFromFile(it, 0, true) ] }
                .set { ch_bed_contigs }

            genome
                .fai
                .map { [ getColFromFile(it, 0, true) ] }
                .concat(ch_bed_contigs)
                .collect()
                .map { fai, bed -> checkContigsInBED(fai, bed, log) }

            // Check whether the primer BED file supplied to the pipeline is from the SWIFT/SNAP protocol
            if (!params.ivar_trim_offset) {
                genome
                    .primer_bed
                    .map { checkIfSwiftProtocol(it, 'covid19genome', log) }
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
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TRIM_FASTP_FASTQC.out.fastqc_raw_zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TRIM_FASTP_FASTQC.out.trim_json.collect{it[1]}.ifEmpty([]))

        //
        // Filter empty FastQ files after adapter trimming
        //
        if (!params.skip_fastp) {
            ch_variants_fastq
                .join(FASTQ_TRIM_FASTP_FASTQC.out.trim_json)
                .map {
                    meta, reads, json ->
                        def pass = getFastpReadsAfterFiltering(json) > 0
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
                        def num_reads = getFastpReadsBeforeFiltering(json)
                        return [ "$meta.id\t$num_reads" ]
                    }
                }
                .collect()
                .map {
                    tsv_data ->
                        def header = ['Sample', 'Reads before trimming']
                        multiqcTsvFromList(tsv_data, header)
                }
                .collectFile(name: 'fail_mapped_reads_mqc.tsv')
                .ifEmpty([])
                .set { ch_fail_reads_multiqc }
            ch_multiqc_files = ch_multiqc_files.mix(ch_fail_reads_multiqc)
        }

        //
        // MODULE: Run Kraken2 for removal of host reads
        //
        ch_assembly_fastq  = ch_variants_fastq
        if (!params.skip_kraken2) {
            KRAKEN2_KRAKEN2 (
                ch_variants_fastq,
                genome.kraken2_db,
                params.kraken2_variants_host_filter || params.kraken2_assembly_host_filter,
                params.kraken2_variants_host_filter || params.kraken2_assembly_host_filter
            )
            ch_multiqc_files =  ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))

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
        ch_bam = channel.empty()
        ch_bai = channel.empty()
        if (!params.skip_variants) {
            FASTQ_ALIGN_BOWTIE2 (
                ch_variants_fastq,
                genome.bowtie2_index,
                params.save_unaligned,
                false,
                ch_reference_fasta_fai
            )
        ch_bam           = FASTQ_ALIGN_BOWTIE2.out.bam
        ch_bai           = FASTQ_ALIGN_BOWTIE2.out.index
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.log_out.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.flagstat.collect{it[1]}.ifEmpty([]))
        }

        //
        // Filter channels to get samples that passed Bowtie2 minimum mapped reads threshold
        //
        ch_fail_mapping_multiqc = channel.empty()
        if (!params.skip_variants) {
            FASTQ_ALIGN_BOWTIE2.out.flagstat
                .map { meta, flagstat -> [ meta ] + getFlagstatMappedReads(flagstat, params) }
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
                        multiqcTsvFromList(tsv_data, header)
                }
                .collectFile(name: 'fail_mapped_samples_mqc.tsv')
                .ifEmpty([])
                .set { ch_fail_mapping_multiqc }
            ch_multiqc_files = ch_multiqc_files.mix(ch_fail_mapping_multiqc)
        }

        //
        // SUBWORKFLOW: Trim primer sequences from reads with iVar
        //
        if (!params.skip_variants && !params.skip_ivar_trim && params.trim_primers) {
            BAM_TRIM_PRIMERS_IVAR (
                ch_bam.join(ch_bai, by: [0]),
                genome.primer_bed,
                ch_reference_fasta_fai
            )
            ch_bam           = BAM_TRIM_PRIMERS_IVAR.out.bam
            ch_bai           = BAM_TRIM_PRIMERS_IVAR.out.bai
            ch_multiqc_files = ch_multiqc_files.mix(BAM_TRIM_PRIMERS_IVAR.out.flagstat.collect{it[1]}.ifEmpty([]))
            ch_versions      = ch_versions.mix(BAM_TRIM_PRIMERS_IVAR.out.versions)
        }

        //
        // SUBWORKFLOW: Mark duplicate reads
        //
        if (!params.skip_variants && !params.skip_markduplicates) {
            BAM_MARKDUPLICATES_PICARD (
                ch_bam,
                ch_reference_fasta_fai
            )
            ch_bam           = BAM_MARKDUPLICATES_PICARD.out.bam
            ch_bai           = BAM_MARKDUPLICATES_PICARD.out.index
            ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect{it[1]}.ifEmpty([]))
        }

        //
        // MODULE: Picard metrics
        //
        if (!params.skip_variants && !params.skip_picard_metrics) {
            PICARD_COLLECTMULTIPLEMETRICS (
                ch_bam.join(ch_bai, by: [0]),
                genome.fasta.map { [ [:], it ] },
                genome.fai.map { [ [:], it ] }
            )
        }

        //
        // MODULE: Genome-wide and amplicon-specific coverage QC plots
        //
        if (!params.skip_variants && !params.skip_mosdepth) {
            MOSDEPTH_GENOME (
                ch_bam
                    .join(ch_bai, by: [0])
                    .map { meta, bam, bai -> [ meta, bam, bai, [] ] },
                [ [:], [] ],
                []
            )
            ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_GENOME.out.global_txt.collect{it[1]}.ifEmpty([]))
            PLOT_MOSDEPTH_REGIONS_GENOME (
                MOSDEPTH_GENOME.out.regions_bed.collect { it[1] }
            )

            if (params.trim_primers) {
                MOSDEPTH_AMPLICON (
                    ch_bam
                        .join(ch_bai, by: [0])
                        .combine(genome.primer_collapsed_bed),
                    [ [:], [] ],
                    []
                )

                PLOT_MOSDEPTH_REGIONS_AMPLICON (
                    MOSDEPTH_AMPLICON.out.regions_bed.collect { it[1] }
                )
                ch_multiqc_files = ch_multiqc_files.mix(PLOT_MOSDEPTH_REGIONS_AMPLICON.out.heatmap_tsv.collect().ifEmpty([]))
            }
        }

        //
        // SUBWORKFLOW: Call variants with IVar
        //
        ch_vcf                    = channel.empty()
        ch_tbi                    = channel.empty()
        ch_snpsift_txt            = channel.empty()
        if (!params.skip_variants && variant_caller == 'ivar') {
            VARIANTS_IVAR (
                ch_bam,
                genome.fasta,
                (params.trim_primers || !params.skip_markduplicates) ? genome.fai : [],
                (params.trim_primers || !params.skip_markduplicates) ? genome.chrom_sizes : [],
                ch_genome_gff ? genome.gff : [],
                (params.trim_primers && ch_primer_bed) ? genome.primer_bed : [],
                genome.snpeff_db,
                genome.snpeff_config,
                ch_ivar_variants_header_mqc
            )
            ch_vcf           = VARIANTS_IVAR.out.vcf
            ch_tbi           = VARIANTS_IVAR.out.tbi
            ch_snpsift_txt   = VARIANTS_IVAR.out.snpsift_txt
            ch_multiqc_files = ch_multiqc_files.mix(VARIANTS_IVAR.out.multiqc_tsv.collect{it[1]}.ifEmpty([]))
            ch_multiqc_files = ch_multiqc_files.mix(VARIANTS_IVAR.out.stats.collect{it[1]}.ifEmpty([]))
            ch_multiqc_files = ch_multiqc_files.mix(VARIANTS_IVAR.out.snpeff_csv.collect{it[1]}.ifEmpty([]))
            ch_versions      = ch_versions.mix(VARIANTS_IVAR.out.versions)
        }

        //
        // SUBWORKFLOW: Call variants with BCFTools
        //
        if (!params.skip_variants && variant_caller == 'bcftools') {
            VARIANTS_BCFTOOLS (
                ch_bam,
                ch_reference_fasta_fai,
                (params.trim_primers || !params.skip_markduplicates) ? genome.chrom_sizes : [],
                ch_genome_gff ? genome.gff : [],
                (params.trim_primers && ch_primer_bed) ? genome.primer_bed : [],
                genome.snpeff_db,
                genome.snpeff_config
            )
            ch_vcf           = VARIANTS_BCFTOOLS.out.vcf
            ch_tbi           = VARIANTS_BCFTOOLS.out.tbi
            ch_multiqc_files = ch_multiqc_files.mix(VARIANTS_BCFTOOLS.out.stats.collect{it[1]}.ifEmpty([]))
            ch_multiqc_files = ch_multiqc_files.mix(VARIANTS_BCFTOOLS.out.snpeff_csv.collect{it[1]}.ifEmpty([]))
            ch_snpsift_txt   = VARIANTS_BCFTOOLS.out.snpsift_txt
        }

        //
        // SUBWORKFLOW: Determine variants with Freyja
        //
        if (!params.skip_variants && !params.skip_freyja) {
            BAM_VARIANT_DEMIX_BOOT_FREYJA(
                ch_bam,
                genome.fasta,
                params.skip_freyja_boot,
                params.freyja_repeats,
                params.freyja_db_name,
                params.freyja_barcodes,
                params.freyja_lineages_meta,
                params.freyja_lineages_topology,
            )
            ch_multiqc_files  = ch_multiqc_files.mix(BAM_VARIANT_DEMIX_BOOT_FREYJA.out.demix.collect{it -> it[1]}.ifEmpty([]))
        }

        //
        // SUBWORKFLOW: Call consensus with iVar and downstream QC
        //
        ch_nextclade_report = channel.empty()
        ch_pangolin_report  = channel.empty()
        ch_consensus_genome = channel.empty()

        if (!params.skip_variants && !params.skip_consensus && params.consensus_caller == 'ivar') {
            CONSENSUS_IVAR (
                ch_bam,
                genome.fasta,
                ch_genome_gff ? genome.gff.map { [ [:], it ] } : [ [:], [] ],
                genome.nextclade_db
            )
            ch_nextclade_report = CONSENSUS_IVAR.out.nextclade_report
            ch_pangolin_report  = CONSENSUS_IVAR.out.pangolin_report
            ch_consensus_genome = CONSENSUS_IVAR.out.consensus
            ch_multiqc_files    = ch_multiqc_files.mix(ch_pangolin_report.collect{it[1]}.ifEmpty([]))
            ch_multiqc_files    = ch_multiqc_files.mix(CONSENSUS_IVAR.out.quast_results.collect{it[1]}.ifEmpty([]))
            ch_versions         = ch_versions.mix(CONSENSUS_IVAR.out.versions)
        }

        //
        // SUBWORKFLOW: Call consensus with BCFTools
        //
        if (!params.skip_variants && !params.skip_consensus && params.consensus_caller == 'bcftools' && variant_caller) {
            CONSENSUS_BCFTOOLS (
                ch_bam,
                ch_vcf,
                ch_tbi,
                genome.fasta,
                ch_genome_gff ? genome.gff.map { [ [:], it ] } : [ [:], [] ],
                genome.nextclade_db
            )

            ch_nextclade_report = CONSENSUS_BCFTOOLS.out.nextclade_report
            ch_pangolin_report  = CONSENSUS_BCFTOOLS.out.pangolin_report
            ch_consensus_genome = CONSENSUS_BCFTOOLS.out.consensus
            ch_multiqc_files    = ch_multiqc_files.mix(CONSENSUS_BCFTOOLS.out.quast_results.collect{it[1]}.ifEmpty([]))
            ch_multiqc_files    = ch_multiqc_files.mix(ch_pangolin_report.collect{it[1]}.ifEmpty([]))
            ch_versions         = ch_versions.mix(CONSENSUS_BCFTOOLS.out.versions)
        }

        //
        // MODULE: Get Nextclade clade information for MultiQC report
        //
        ch_nextclade_multiqc = channel.empty()
        if (!params.skip_variants && !params.skip_nextclade) {
            ch_nextclade_report
                .map { meta, csv ->
                    def clade = getNextcladeFieldMapFromCsv(csv)['clade']
                    return [ "$meta.id\t$clade" ]
                }
                .collect()
                .map {
                    tsv_data ->
                        def header = ['Sample', 'clade']
                        multiqcTsvFromList(tsv_data, header)
                }
                .collectFile(name: 'nextclade_clade_mqc.tsv')
                .ifEmpty([])
                .set { ch_nextclade_multiqc }
            ch_multiqc_files = ch_multiqc_files.mix(ch_nextclade_multiqc)
        }

        //
        // SUBWORKFLOW: Create variants long table report
        //
        if (!params.skip_variants && !params.skip_variants_long_table && ch_genome_gff && !params.skip_snpeff) {
            VARIANTS_LONG_TABLE (
                ch_vcf,
                ch_tbi,
                ch_snpsift_txt,
                ch_pangolin_report
            )
        }

        //
        // SUBWORKFLOW: Create variants long table report for additional annotation file
        //
        if (!params.skip_variants && params.additional_annotation) {
            ch_annot = channel.empty()
            //
            // Uncompress additional annotation file
            //
            if (params.additional_annotation.endsWith('.gz')) {
                GUNZIP_GFF (
                    [ [:], ch_additional_gtf ]
                )
                ch_annot       = GUNZIP_GFF.out.gunzip.map { it[1] }
            } else {
                ch_annot = ch_additional_gtf
            }

            ADDITIONAL_ANNOTATION (
                ch_vcf,
                ch_tbi,
                genome.fasta,
                ch_annot,
                ch_pangolin_report

            )
        }

        //
        // SUBWORKFLOW: HIV resistance detection
        //

        if (!params.skip_variants && params.perform_hiv_resistance) {
            HIV_RESISTANCE (
                ch_consensus_genome,
                ch_bam.join(ch_bai, by: [0]),
                genome.fasta,
                genome.gff,
                ch_vcf,
                ch_tbi,
                ch_pangolin_report,
                ch_nextclade_report
            )
        }

        //
        // MODULE: Primer trimming with Cutadapt
        //
        if (params.trim_primers && !params.skip_assembly && !params.skip_cutadapt) {
            ch_primers =  genome.primer_fasta
            if (!params.skip_noninternal_primers){
                PREPARE_PRIMER_FASTA(
                    genome.primer_fasta
                )
                ch_primers = PREPARE_PRIMER_FASTA.out.adapters
            }

            CUTADAPT (
                ch_assembly_fastq,
                ch_primers
            )
            ch_assembly_fastq   = CUTADAPT.out.reads
            ch_multiqc_files    = ch_multiqc_files.mix(CUTADAPT.out.log.collect{it[1]}.ifEmpty([]))

            if (!params.skip_fastqc) {
                FASTQC (
                    CUTADAPT.out.reads
                )
            }
        }

        //
        // SUBWORKFLOW: Run SPAdes assembly and downstream analysis
        //
        if (!params.skip_assembly && 'spades' in assemblers) {
            ASSEMBLY_SPADES (
                ch_assembly_fastq.map { meta, fastq -> [ meta, fastq, [], [] ] },
                params.spades_mode,
                ch_spades_hmm,
                genome.fasta,
                ch_genome_gff ? genome.gff.map { [ [:], it ] } : [ [:], [] ],
                genome.blast_db,
                ch_blast_outfmt6_header,
                ch_blast_filtered_outfmt6_header,
                ch_taxidlist
            )
            ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLY_SPADES.out.quast_results.collect{it[1]}.ifEmpty([]))
            ch_versions      = ch_versions.mix(ASSEMBLY_SPADES.out.versions)
        }

        //
        // SUBWORKFLOW: Run Unicycler assembly and downstream analysis
        //
        if (!params.skip_assembly && 'unicycler' in assemblers) {
            ASSEMBLY_UNICYCLER (
                ch_assembly_fastq.map { meta, fastq -> [ meta, fastq, [] ] },
                genome.fasta,
                ch_genome_gff ? genome.gff.map { [ [:], it ] } : [ [:], [] ],
                genome.blast_db,
                ch_blast_outfmt6_header,
                ch_blast_filtered_outfmt6_header,
                ch_taxidlist
            )
            ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLY_UNICYCLER.out.quast_results.collect{it[1]}.ifEmpty([]))
            ch_versions      = ch_versions.mix(ASSEMBLY_UNICYCLER.out.versions)
        }

        //
        // SUBWORKFLOW: Run minia assembly and downstream analysis
        //
        if (!params.skip_assembly && 'minia' in assemblers) {
            ASSEMBLY_MINIA (
                ch_assembly_fastq,
                genome.fasta,
                ch_genome_gff ? genome.gff.map { [ [:], it ] } : [ [:], [] ],
                genome.blast_db,
                ch_blast_outfmt6_header,
                ch_blast_filtered_outfmt6_header,
                ch_taxidlist
            )
            ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLY_MINIA.out.quast_results.collect{it[1]}.ifEmpty([]))
            ch_versions      = ch_versions.mix(ASSEMBLY_MINIA.out.versions)
        }

    } else if (params.platform == 'nanopore') {
        //
        // NANOPORE WORKFLOW
        //
        def ch_fasta_fai_nanopore = genome.fasta
            .combine(genome.fai)
            .map { fasta_file, fai_file -> [ [:], fasta_file, fai_file ] }
            .collect(flat: false)
            .map { fasta_fai -> fasta_fai[0] }
        def ch_fasta_primer_bed_nanopore = genome.fasta
            .combine(genome.primer_bed)
            .map { fasta_file, primer_bed -> [ [:], fasta_file, primer_bed ] }
            .collect(flat: false)
            .map { fasta_primer_bed -> fasta_primer_bed[0] }
        def ch_gff_tuple_nanopore = ch_genome_gff ? genome.gff.map { gff_file -> [ [:], gff_file ] } : [ [:], [] ]
        def min_barcode_reads = params.min_barcode_reads as Integer
        def min_guppyplex_reads = params.min_guppyplex_reads as Integer

        //
        // MODULE: PycoQC on sequencing summary file
        //
        if (sequencing_summary && !params.skip_pycoqc) {
            PYCOQC (
                channel.of(ch_sequencing_summary).map { [ [:], it ] }
            )
            ch_multiqc_files = ch_multiqc_files.mix(PYCOQC.out.json.collect{it[1]}.ifEmpty([]))
        }

        // Check primer BED file only contains suffixes provided --primer_left_suffix / --primer_right_suffix
        genome
            .primer_bed
            .map { checkPrimerSuffixes(it, params.primer_left_suffix, params.primer_right_suffix, log) }

        // Check whether the contigs in the primer BED file are present in the reference genome
        genome
            .primer_bed
            .map { [ getColFromFile(it, 0, true) ] }
            .set { ch_bed_contigs }

        genome
            .fai
            .map { [ getColFromFile(it, 0, true) ] }
            .concat(ch_bed_contigs)
            .collect()
            .map { fai, bed -> checkContigsInBED(fai, bed, log) }

        barcode_dirs       = file("${params.fastq_dir}/barcode*", type: 'dir' , maxdepth: 1)
        single_barcode_dir = file("${params.fastq_dir}/*.fastq" , type: 'file', maxdepth: 1)
        if (barcode_dirs) {
            channel
                .fromPath( barcode_dirs )
                .filter( ~/.*barcode[0-9]{1,4}$/ )
                .map { dir ->
                    def count = 0
                    dir.listFiles().each { x ->
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
                    .filter { it[-1] >= min_barcode_reads }
                    .map { it -> [ "${it[0]}\t${it[-1]}" ] }
                    .collect()
                    .map {
                        tsv_data ->
                            def header = ['Barcode', 'Read count']
                            multiqcTsvFromList(tsv_data, header)
                    }
                    .collectFile(name: 'fail_barcodes_no_sample_mqc.tsv')
                    .ifEmpty([])
                    .set { ch_custom_no_sample_name_multiqc }

                ch_multiqc_files = ch_multiqc_files.mix ( ch_custom_no_sample_name_multiqc )
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
                            multiqcTsvFromList(tsv_data, header)
                    }
                    .collectFile(name: 'fail_no_barcode_samples_mqc.tsv')
                    .ifEmpty([])
                    .set { ch_custom_no_barcodes_multiqc }

                ch_multiqc_files = ch_multiqc_files.mix ( ch_custom_no_barcodes_multiqc )

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
            channel
                .fromPath("${params.fastq_dir}", type: 'dir', maxDepth: 1)
                .map { it -> [ 'SAMPLE_1', 'single_barcode', it, 10000000 ] }
                .set{ ch_fastq_dirs }
        } else {
            error "Please specify a valid folder containing ONT basecalled, barcoded fastq files generated by guppy_barcoder or guppy_basecaller e.g. '--fastq_dir ./20191023_1522_MC-110615_0_FAO93606_12bf9b4f/fastq_pass/"
        }

        //
        // MODULE: Create custom content file for MultiQC to report samples with reads < params.min_barcode_reads
        //
        ch_fastq_dirs
            .branch { barcode, sample, dir, count  ->
                pass: count > min_barcode_reads
                    pass_barcode_reads[sample] = count
                    return [ "$sample\t$count" ]
                fail: count < min_barcode_reads
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
                    multiqcTsvFromList(tsv_data, header)
            }
            .collectFile(name: 'fail_barcode_count_samples_mqc.tsv')
            .ifEmpty([])
            .set { ch_custom_fail_barcodes_count_multiqc }

        ch_multiqc_files = ch_multiqc_files.mix(ch_custom_fail_barcodes_count_multiqc)

        // Re-arrange channels to have meta map of information for sample
        ch_fastq_dirs
            .filter { it[-1] > min_barcode_reads }
            .map { barcode, sample, dir, count -> [ [ id: sample, barcode:barcode ], dir ] }
            .set { ch_fastq_dirs }

        //
        // MODULE: Run Artic Guppyplex
        //
        ARTIC_GUPPYPLEX (
            ch_fastq_dirs
        )

        //
        // MODULE: Run Kraken2 for removal of host reads
        //
        ch_variants_fastq = ARTIC_GUPPYPLEX.out.fastq.map { meta, fastq ->
                    meta += [single_end: true]
                    return [meta, fastq]
                }
        ch_assembly_fastq  = ch_variants_fastq
        if (!params.skip_kraken2) {
            KRAKEN2_KRAKEN2 (
                ch_variants_fastq,
                genome.kraken2_db,
                params.kraken2_variants_host_filter || params.kraken2_assembly_host_filter,
                params.kraken2_variants_host_filter || params.kraken2_assembly_host_filter
            )
            ch_multiqc_files =  ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))

            if (params.kraken2_variants_host_filter) {
                ch_variants_fastq = KRAKEN2_KRAKEN2.out.unclassified_reads_fastq
            }

            if (params.kraken2_assembly_host_filter) {
                ch_assembly_fastq = KRAKEN2_KRAKEN2.out.unclassified_reads_fastq
            }
        }

        //
        // MODULE: Create custom content file for MultiQC to report samples with reads < params.min_guppyplex_reads
        //
        ARTIC_GUPPYPLEX
            .out
            .fastq
            .branch { meta, fastq  ->
                def count = fastq.countFastq()
                pass: count > min_guppyplex_reads
                    return [ "$meta.id\t$count" ]
                fail: count < min_guppyplex_reads
                    return [ "$meta.id\t$count" ]
            }
            .set { ch_pass_fail_guppyplex_count }

        ch_pass_fail_guppyplex_count
            .fail
            .collect()
            .map {
                tsv_data ->
                    def header = ['Sample', 'Read count']
                    multiqcTsvFromList(tsv_data, header)
            }
            .collectFile(name: 'fail_guppyplex_count_samples_mqc.tsv')
            .ifEmpty([])
            .set { ch_custom_fail_guppyplex_count_multiqc }

        ch_multiqc_files = ch_multiqc_files.mix(ch_custom_fail_guppyplex_count_multiqc)

        //
        // MODULE: Nanoplot QC for FastQ files
        //
        if (!params.skip_nanoplot) {
            NANOPLOT (
                ARTIC_GUPPYPLEX.out.fastq
            )
            ch_versions = ch_versions.mix(NANOPLOT.out.versions)
        }

        //
        // MODULE: Run Nanopore mapping, variant calling and consensus generation
        //
        ch_bam       = channel.empty()
        ch_bai       = channel.empty()
        ch_vcf       = channel.empty()
        ch_tbi       = channel.empty()
        ch_consensus = channel.empty()

        if (params.mapper_nanopore == 'artic') {

            ARTIC_MINION_PROTOCOL (
                ARTIC_GUPPYPLEX.out.fastq.filter { it[-1].countFastq() > min_guppyplex_reads },
                ch_clair3_model,
                ch_fasta_primer_bed_nanopore
            )

            ch_bam  = ARTIC_MINION_PROTOCOL.out.bam
            ch_bai  = ARTIC_MINION_PROTOCOL.out.bai

            ch_vcf  = ARTIC_MINION_PROTOCOL.out.vcf
            ch_tbi  = ARTIC_MINION_PROTOCOL.out.tbi

            ch_consensus = ARTIC_MINION_PROTOCOL.out.consensus

            ch_multiqc_files = ch_multiqc_files.mix(ARTIC_MINION_PROTOCOL.out.artic_minion_report.collect{it[1]}.ifEmpty([]))
            ch_versions      = ch_versions.mix(ARTIC_MINION_PROTOCOL.out.versions)

        } else if (params.mapper_nanopore == 'minimap2') {

            MINIMAP2_MAPPING(
                ARTIC_GUPPYPLEX.out.fastq.filter { it[-1].countFastq() > min_guppyplex_reads },
                genome.fasta,
                genome.fai,
                genome.primer_bed
            )

            ch_bam  = MINIMAP2_MAPPING.out.bam
            ch_bai  = MINIMAP2_MAPPING.out.bai

            ch_vcf  = MINIMAP2_MAPPING.out.vcf
            ch_tbi  = MINIMAP2_MAPPING.out.tbi

            ch_consensus = MINIMAP2_MAPPING.out.consensus

            ch_multiqc_files = ch_multiqc_files.mix(MINIMAP2_MAPPING.out.multiqc_files.collect{it[1]}.ifEmpty([]))

            ch_versions      = ch_versions.mix(MINIMAP2_MAPPING.out.versions)

        }

        ch_bam_bai = ch_bam.join(ch_bai, by: [0])
        ch_vcf_tbi = ch_vcf.join(ch_tbi, by: [0])

        //
        // MODULE: VCF stats with bcftools stats
        //
        BCFTOOLS_STATS (
            ch_vcf_tbi,
            [ [:], [] ],
            [ [:], [] ],
            [ [:], [] ],
            [ [:], [] ],
            [ [:], [] ]
        )
        ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS.out.stats.collect{it[1]}.ifEmpty([]))

        //
        // SUBWORKFLOW: Filter unmapped reads from BAM
        //
        FILTER_BAM_SAMTOOLS (
            ch_bam_bai,
            ch_fasta_fai_nanopore
        )
        ch_multiqc_files = ch_multiqc_files.mix(FILTER_BAM_SAMTOOLS.out.flagstat.collect{it[1]}.ifEmpty([]))

        //
        // Filter channels to get samples that passed minimum mapped reads threshold
        //
        ch_fail_mapping_multiqc_nanopore = channel.empty()
        FILTER_BAM_SAMTOOLS.out.flagstat
            .map { meta, flagstat ->
                def (mapped_reads, pass) = getFlagstatMappedReads(flagstat, params)
                [ meta, mapped_reads, pass ]
            }
            .set { ch_mapped_reads_nanopore }

        // Filter BAM files based on mapping threshold
        ch_bam
            .join(ch_mapped_reads_nanopore, by: [0])
            .map { meta, bam, mapped, pass ->
                if (pass) [ meta, bam ]
            }
            .set { ch_filtered_bam_nanopore }

        // Filter BAI files based on mapping threshold
        ch_bai
            .join(ch_mapped_reads_nanopore, by: [0])
            .map { meta, bai, mapped, pass ->
                if (pass) [ meta, bai ]
            }
            .set { ch_filtered_bai_nanopore }

        // Filter FASTA files based on mapping threshold
        ch_consensus
            .join(ch_mapped_reads_nanopore, by: [0])
            .map { meta, fasta, mapped_reads, pass ->
                if (pass) [ meta, fasta ]
            }
            .set { ch_filtered_fasta_nanopore }

        // Track passed/failed samples for MultiQC
        ch_mapped_reads_nanopore
            .branch { meta, mapped, pass ->
                pass: pass
                    pass_mapped_reads[meta.id] = mapped
                    return [ "$meta.id\t$mapped" ]
                fail: !pass
                    fail_mapped_reads[meta.id] = mapped
                    return [ "$meta.id\t$mapped" ]
            }
            .set { ch_pass_fail_mapped_nanopore }

        // Create MultiQC report for failed samples
        ch_pass_fail_mapped_nanopore
            .fail
            .collect()
            .map {
                tsv_data ->
                    def header = ['Sample', 'Mapped reads']
                    multiqcTsvFromList(tsv_data, header)
            }
            .collectFile(name: 'fail_mapped_samples_nanopore_mqc.tsv')
            .ifEmpty([])
            .set { ch_fail_mapping_multiqc_nanopore }
        ch_multiqc_files = ch_multiqc_files.mix(ch_fail_mapping_multiqc_nanopore)

        //
        // MODULE: Genome-wide and amplicon-specific coverage QC plots
        //
        if (!params.skip_mosdepth) {

            MOSDEPTH_GENOME (
                ch_filtered_bam_nanopore
                    .join(ch_filtered_bai_nanopore, by: [0])
                    .map { meta, bam, bai -> [ meta, bam, bai, [] ] },
                [ [:], [] ],
                []
            )
            ch_multiqc_files  = ch_multiqc_files.mix(MOSDEPTH_GENOME.out.global_txt.collect{it[1]}.ifEmpty([]))

            PLOT_MOSDEPTH_REGIONS_GENOME (
                MOSDEPTH_GENOME.out.regions_bed.collect { it[1] }
            )

            MOSDEPTH_AMPLICON (
                ch_filtered_bam_nanopore
                    .join(ch_filtered_bai_nanopore, by: [0])
                    .combine(genome.primer_collapsed_bed),
                [ [:], [] ],
                []
            )


            PLOT_MOSDEPTH_REGIONS_AMPLICON (
                MOSDEPTH_AMPLICON.out.regions_bed.collect { it[1] }
            )
            ch_multiqc_files = ch_multiqc_files.mix(PLOT_MOSDEPTH_REGIONS_AMPLICON.out.heatmap_tsv.collect{it[1]}.ifEmpty([]))
        }

        //
        // MODULE: Lineage analysis with Pangolin
        //
        ch_pango_database = channel.empty()
        ch_pangolin_report = channel.empty()

        if (!params.skip_pangolin) {
            if (!params.pango_database) {
                PANGOLIN_UPDATEDATA('pangolin_db')
                ch_pango_database = PANGOLIN_UPDATEDATA.out.db
                ch_versions       = ch_versions.mix(PANGOLIN_UPDATEDATA.out.versions)
            } else {
                if (params.pango_database.endsWith('.tar.gz')) {
                    UNTAR_PANGODB (
                        [ [:], params.pango_database ]
                    )
                    ch_pango_database = UNTAR_PANGODB.out.untar.map { it[1] }
                } else {
                    ch_pango_database = channel.value(file(params.pango_database, type: 'dir'))
                }
            }
            def ch_pango_database_for_run = ch_pango_database.collect(flat: false).map { dirs -> dirs[0] }

            PANGOLIN_RUN (
                ch_consensus,
                ch_pango_database_for_run
            )
            ch_pangolin_multiqc = PANGOLIN_RUN.out.report
            ch_multiqc_files    = ch_multiqc_files.mix(ch_pangolin_multiqc.collect{it[1]}.ifEmpty([]))
            ch_versions         = ch_versions.mix(PANGOLIN_RUN.out.versions)
        }

        //
        // MODULE: Clade assignment, mutation calling, and sequence quality checks with Nextclade
        //
        if (!params.skip_nextclade) {
            NEXTCLADE_RUN (
                ch_consensus,
                genome.nextclade_db
            )
            ch_versions = ch_versions.mix(NEXTCLADE_RUN.out.versions)

            //
            // MODULE: Get Nextclade clade information for MultiQC report
            //
            NEXTCLADE_RUN
                .out
                .csv
                .map {
                    meta, csv ->
                        def clade = getNextcladeFieldMapFromCsv(csv)['clade']
                        return [ "$meta.id\t$clade" ]
                }
                .collect()
                .map {
                    tsv_data ->
                        def header = ['Sample', 'clade']
                        multiqcTsvFromList(tsv_data, header)
                }
                .collectFile(name: 'nextclade_clade_mqc.tsv')
                .ifEmpty([])
                .set{ nextclade_clade_mqc }

            ch_multiqc_files = ch_multiqc_files.mix(nextclade_clade_mqc)
        }

        //
        // SUBWORKFLOW: Determine variants with Freyja
        //
        if (!params.skip_freyja) {
            BAM_VARIANT_DEMIX_BOOT_FREYJA(
                ch_filtered_bam_nanopore,
                genome.fasta,
                params.skip_freyja_boot,
                params.freyja_repeats,
                params.freyja_db_name,
                params.freyja_barcodes,
                params.freyja_lineages_meta,
                params.freyja_lineages_topology,
            )
            ch_multiqc_files  = ch_multiqc_files.mix(BAM_VARIANT_DEMIX_BOOT_FREYJA.out.demix.collect{it[1]}.ifEmpty([]))
        }

        //
        // MODULE: Consensus QC across all samples with QUAST
        //
        if (!params.skip_variants_quast) {
            ch_filtered_fasta_nanopore
                .collect{ it[1] }
                .map { consensus_collect -> tuple([id: "quast"], consensus_collect) }
                .set { ch_to_quast }
            QUAST (
                ch_to_quast,
                genome.fasta.map { fasta_file -> [ [:], fasta_file ] },
                ch_gff_tuple_nanopore,
            )
            ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.results.collect{it[1]}.ifEmpty([]))
        }

        //
        // SUBWORKFLOW: Annotate variants with snpEff
        //
        ch_snpsift_txt    = channel.empty()
        if (ch_genome_gff && !params.skip_snpeff) {
            SNPEFF_SNPSIFT (
                ch_vcf,
                genome.snpeff_db,
                genome.snpeff_config,
                genome.fasta
            )
            ch_multiqc_files  = ch_multiqc_files.mix(SNPEFF_SNPSIFT.out.csv.collect{it[1]}.ifEmpty([]))
            ch_snpsift_txt    = SNPEFF_SNPSIFT.out.snpsift_txt
        }

        //
        // SUBWORKFLOW: Create variants long table report
        //
        if (!params.skip_variants_long_table && ch_genome_gff && !params.skip_snpeff) {
            VARIANTS_LONG_TABLE (
                ch_vcf,
                ch_tbi,
                ch_snpsift_txt,
                ch_pangolin_multiqc
            )
        }

        //
        // SUBWORKFLOW: Create variants long table report for additional annotation file
        //
        if (params.additional_annotation) {
            ch_annot = channel.empty()
            //
            // Uncompress additional annotation file
            //
            if (params.additional_annotation.endsWith('.gz')) {
                GUNZIP_GFF (
                    [ [:], ch_additional_gtf ]
                )
                ch_annot       = GUNZIP_GFF.out.gunzip.map { it[1] }
            } else {
                ch_annot = ch_additional_gtf
            }

            ADDITIONAL_ANNOTATION (
                ch_vcf,
                ch_tbi,
                genome.fasta,
                ch_annot,
                ch_pangolin_multiqc

            )
        }
    }

    //
    // MODULE: Pipeline reporting
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'viralrecon_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    if (params.platform == 'illumina') {
        ch_multiqc_config        = [ file("$projectDir/assets/multiqc_config_illumina.yml", checkIfExists: true) ]
    } else if (params.platform == 'nanopore') {
        ch_multiqc_config        = [ file("$projectDir/assets/multiqc_config_nanopore.yml", checkIfExists: true) ]
    }
    ch_multiqc_custom_config = params.multiqc_config ? [ file(params.multiqc_config, checkIfExists: true) ] : []
    ch_multiqc_logo          = params.multiqc_logo ? [ file(params.multiqc_logo, checkIfExists: true) ] : []

    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    ch_multiqc_input = ch_multiqc_files
        .collect()
        .map { multiqc_files ->
            [ [:], multiqc_files, ch_multiqc_config, ch_multiqc_logo, [], [], ch_multiqc_custom_config ]
        }

    MULTIQC (
        ch_multiqc_input
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channeel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
