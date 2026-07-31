//
// Subworkflow with functionality specific to the nf-core/viralrecon pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/viralrecon ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/viralrecon/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //
    if (params.platform == 'illumina') {
        channel
            .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2, _barcode->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map { sample_group -> validateInputSamplesheet(sample_group) }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }
    } else if (params.platform == 'nanopore') {
        ch_samplesheet = channel.empty()
        if (input){
            channel
                .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
                .map {
                    meta, _fastq_1, _fastq_2, barcode->
                        tuple( "barcode"+ String.format('%02d', barcode).toString(), meta.id)
                }
                .set { ch_samplesheet }
        }
    }

    emit:
    samplesheet = ch_samplesheet
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// Function to check whether primer BED file has the correct suffixes as provided to the pipeline
//
def checkPrimerSuffixes(primer_bed_file, primer_left_suffix, primer_right_suffix, log) {
    def total = 0
    def left  = 0
    def right = 0
    primer_bed_file.eachLine { line ->
        total += 1
        def name = line.split('\t')[3]
        if (name.contains(primer_left_suffix)) {
            left += 1
        } else if (name.contains(primer_right_suffix)) {
            right += 1
        }
    }
    if (total != (left + right)) {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Please check the name field (column 4) in the file supplied via --primer_bed.\n\n" +
            "  All of the values in that column do not end with those supplied by:\n" +
            "      --primer_left_suffix : $primer_left_suffix\n" +
            "      --primer_right_suffix: $primer_right_suffix\n\n" +
            "  This information is required to collapse the primer intervals into amplicons\n" +
            "  for the coverage plots generated by the pipeline.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    }
}

//
// Create MultiQC tsv custom content from a list of values
//
def multiqcTsvFromList(tsv_data, header) {
    def tsv_string = ""
    if (tsv_data.size() > 0) {
        tsv_string += "${header.join('\t')}\n"
        tsv_string += tsv_data.join('\n')
    }
    return tsv_string
}

//
// Function to get column entries from a file
//
def getColFromFile(input_file, col=0, uniqify=false, sep='\t') {
    def vals = []
    input_file.eachLine { line ->
        def val = line.split(sep)[col]
        if (uniqify) {
            if (!vals.contains(val)) {
                vals << val
            }
        } else {
            vals << val
        }
    }
    return vals
}

//
// Function that returns the number of lines in a file
//
def getNumLinesInFile(input_file) {
    def num_lines = 0
    input_file.eachLine { _line ->
        num_lines += 1
    }
    return num_lines
}

//
// Function to generate an error if contigs in BED file do not match those in reference genome
//
def checkContigsInBED(fai_contigs, bed_contigs, _log) {
    def intersect = bed_contigs.intersect(fai_contigs)
    if (intersect.size() != bed_contigs.size()) {
        def diff = bed_contigs.minus(intersect).sort()
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Contigs in primer BED file do not match those in the reference genome:\n\n" +
            "  ${diff.join('\n  ')}\n\n" +
            "  Please check:\n" +
            "    - Primer BED file supplied with --primer_bed\n" +
            "    - Genome FASTA file supplied with --fasta\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Function to read in all fields into a Groovy Map from Nextclade CSV output file
//
// See: https://stackoverflow.com/a/67766919
def getNextcladeFieldMapFromCsv(nextclade_report) {
    def headers   = []
    def field_map = [:]
    nextclade_report.readLines().eachWithIndex { row, row_index ->
        def vals = row.split(';')
        if (row_index == 0) {
            headers = vals
        } else {
            headers.eachWithIndex { header, header_index ->
                def val = (header_index <= vals.size()-1) ? vals[header_index] : ''
                field_map[header] = val ?: 'NA'
            }
        }
    }
    return field_map
}

//
// Function to get number of variants reported in BCFTools stats file
//
def getNumVariantsFromBCFToolsStats(bcftools_stats) {
    def num_vars = 0
    bcftools_stats.eachLine { line ->
        def matcher = line =~ /SN\s*0\s*number\sof\srecords:\s*([\d]+)/
        if (matcher) num_vars = matcher[0][1].toInteger()
    }
    return num_vars
}

//
// Function that parses and returns the number of mapped reads from flagstat files
//
def getFlagstatMappedReads(flagstat_file, params) {
    def mapped_reads = 0
    flagstat_file.eachLine { line ->
        if (line.contains(' mapped (')) {
            mapped_reads = line.tokenize().first().toInteger()
        }
    }

    def pass = false
    if (mapped_reads > params.min_mapped_reads.toInteger()) {
        pass = true
    }
    return [ mapped_reads, pass ]
}


/// FUNCTIONS WORKFLOW ILLUMINA

//
// Print warning if genome fasta has more than one sequence
//
def isMultiFasta(fasta_file, log) {
    def count = 0
    fasta_file.eachLine { line ->
        if (count > 1)
            return
        if (line.contains('>')) {
            count += 1
            if (count > 1) {
                log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  This pipeline does not officially support multi-fasta genome files!\n\n" +
                    "  The parameters and processes are tailored for viral genome analysis.\n" +
                    "  Please amend the '--fasta' parameter.\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            }
        }
    }
}

//
// Check if the primer BED file supplied to the pipeline is from the SWIFT/SNAP protocol
//
def checkIfSwiftProtocol(primer_bed_file, name_prefix, log) {
    def count = 0
    primer_bed_file.eachLine { line ->
        if (count > 1)
            return
        def name = line.split('\t')[3]
        if (name.contains(name_prefix)) {
            count += 1
            if (count > 1) {
                log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                    "  Found '${name_prefix}' in the name field of the primer BED file!\n" +
                    "  This suggests that you have used the SWIFT/SNAP protocol to prep your samples.\n" +
                    "  If so, please set '--ivar_trim_offset 5' as suggested in the issue below:\n" +
                    "  https://github.com/nf-core/viralrecon/issues/170\n" +
                    "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            }
        }
    }
}

//
// Function that parses fastp json output file to get total number of reads after trimming
//
def getFastpReadsAfterFiltering(json_file) {
    def Map json = new groovy.json.JsonSlurper().parseText(json_file.text).get('summary')
    return json['after_filtering']['total_reads'].toInteger()
}

//
// Function that parses fastp json output file to get total number of reads before trimming
//
def getFastpReadsBeforeFiltering(json_file) {
    def Map json = new groovy.json.JsonSlurper().parseText(json_file.text).get('summary')
    return json['before_filtering']['total_reads'].toInteger()
}
