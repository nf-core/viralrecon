// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process MULTIQC_CUSTOM_CSV_FROM_MAP {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    memory 100.MB

    input:
    val csv_data
    val out_prefix

    output:
    path "*.csv"

    exec:
    // Write to file
    def file = task.workDir.resolve("${out_prefix}_mqc.csv")
    file.write csv_data[0].keySet().join(",") + '\n'
    csv_data.each { data ->
        file.append(data.values().join(",") + '\n')
    }
}
