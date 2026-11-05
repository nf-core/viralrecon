process SIERRALOCAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sierra-local:0.4.3--py310hdfd78af_0' :
        'biocontainers/sierra-local:0.4.3--py310hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path hivdb_xml
    path apobec_drm
    path apobec_csv
    path unusual_csv
    path sdrms_csv
    path mutation_csv

    output:
    tuple val(meta), path("*_resistance.json") , emit: json
    tuple val("${task.process}"), val('sierra-local'), val('0.4.3'), emit: versions_sierralocal, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix              = task.ext.prefix     ?: "${meta.id}"
    def args                = task.ext.args       ?: ''
    def hivdb_xml_arg       = params.hivdb_xml    ? "-xml  ${params.hivdb_xml}"            : ''
    def apobec_drm_arg      = params.apobec_drm   ? "-json ${params.apobec_drm}"           : ''
    def apobec_csv_arg      = params.apobec_csv   ? "-apobec_csv ${params.apobec_csv}"     : ''
    def unusual_csv_arg     = params.unusual_csv  ? "-unusual_csv ${params.unusual_csv}"   : ''
    def sdrms_csv_arg       = params.sdrms_csv    ? "-sdrms_csv ${params.sdrms_csv}"       : ''
    def mutation_csv_arg    = params.mutation_csv ? "-mutation_csv ${params.mutation_csv}" : ''

    """
    sierralocal \\
        $args \\
        $hivdb_xml_arg \\
        $apobec_drm_arg \\
        $apobec_csv_arg \\
        $unusual_csv_arg \\
        $sdrms_csv_arg \\
        $mutation_csv_arg \\
        -o ${prefix}_resistance.json \\
        ${fasta}
    """
}
