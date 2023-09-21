//
// This file holds several functions specific to the workflow/illumina.nf in the nf-core/viralrecon pipeline
//
import nextflow.Nextflow
import groovy.json.JsonSlurper

class WorkflowIllumina {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, valid_params) {
        WorkflowCommons.genomeExistsError(params, log)

        // Generic parameter validation
        if (!valid_params['protocols'].contains(params.protocol)) {
            Nextflow.error("Invalid option: '${params.protocol}'. Valid options for '--protocol': ${valid_params['protocols'].join(', ')}.")
        }

        if (!params.fasta) {
            Nextflow.error("Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file.")
        }

        if (!params.skip_kraken2 && !params.kraken2_db) {
            if (!params.kraken2_db_name) {
                Nextflow.error("Please specify a valid name to build Kraken2 database for host e.g. '--kraken2_db_name human'.")
            }
        }

        // Variant calling parameter validation
        if (params.variant_caller) {
            if (!valid_params['variant_callers'].contains(params.variant_caller)) {
                Nextflow.error("Invalid option: ${params.variant_caller}. Valid options for '--variant_caller': ${valid_params['variant_callers'].join(', ')}.")
            }
        }

        // Consensus calling parameter validation
        if (params.consensus_caller) {
            if (!valid_params['consensus_callers'].contains(params.consensus_caller)) {
                Nextflow.error("Invalid option: ${params.consensus_caller}. Valid options for '--consensus_caller': ${valid_params['consensus_callers'].join(', ')}.")
            }
        }

        if (params.protocol == 'amplicon' && !params.skip_variants && !params.primer_bed) {
            Nextflow.error("To perform variant calling in amplicon mode please provide a valid primer BED file e.g. '--primer_bed primers.bed'.")
        }

        // Assembly parameter validation
        def assemblers = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []
        if ((valid_params['assemblers'] + assemblers).unique().size() != valid_params['assemblers'].size()) {
            Nextflow.error("Invalid option: ${params.assemblers}. Valid options for '--assemblers': ${valid_params['assemblers'].join(', ')}.")
        }

        if (!valid_params['spades_modes'].contains(params.spades_mode)) {
            Nextflow.error("Invalid option: ${params.spades_mode}. Valid options for '--spades_modes': ${valid_params['spades_modes'].join(', ')}.")
        }
    }

    //
    // Print warning if genome fasta has more than one sequence
    //
    public static void isMultiFasta(fasta_file, log) {
        def count = 0
        def line  = null
        fasta_file.withReader { reader ->
            while (line = reader.readLine()) {
                if (line.contains('>')) {
                    count++
                    if (count > 1) {
                        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                            "  This pipeline does not officially support multi-fasta genome files!\n\n" +
                            "  The parameters and processes are tailored for viral genome analysis.\n" +
                            "  Please amend the '--fasta' parameter.\n" +
                            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                        break
                    }
                }
            }
        }
    }

    //
    // Function that parses and returns the number of mapped reasds from flagstat files
    //
    public static ArrayList getFlagstatMappedReads(flagstat_file, params) {
        def mapped_reads = 0
        flagstat_file.eachLine { line ->
            if (line.contains(' mapped (')) {
                mapped_reads = line.tokenize().first().toInteger()
            }
        }

        def pass = false
        def logname = flagstat_file.getBaseName() - 'flagstat'
        if (mapped_reads > params.min_mapped_reads.toInteger()) {
            pass = true
        }
        return [ mapped_reads, pass ]
    }

    //
    // Check if the primer BED file supplied to the pipeline is from the SWIFT/SNAP protocol
    //
    public static void checkIfSwiftProtocol(primer_bed_file, name_prefix, log) {
        def count = 0
        def line  = null
        primer_bed_file.withReader { reader ->
            while (line = reader.readLine()) {
                def name = line.split('\t')[3]
                if (name.contains(name_prefix)) {
                    count++
                    if (count > 1) {
                        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                            "  Found '${name_prefix}' in the name field of the primer BED file!\n" +
                            "  This suggests that you have used the SWIFT/SNAP protocol to prep your samples.\n" +
                            "  If so, please set '--ivar_trim_offset 5' as suggested in the issue below:\n" +
                            "  https://github.com/nf-core/viralrecon/issues/170\n" +
                            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                        break
                    }
                }
            }
        }
    }

    //
    // Function that parses fastp json output file to get total number of reads after trimming
    //
    public static Integer getFastpReadsAfterFiltering(json_file) {
        def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
        return json['after_filtering']['total_reads'].toInteger()
    }

    //
    // Function that parses fastp json output file to get total number of reads before trimming
    //
    public static Integer getFastpReadsBeforeFiltering(json_file) {
        def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
        return json['before_filtering']['total_reads'].toInteger()
    }
}
