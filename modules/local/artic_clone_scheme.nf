// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process ARTIC_CLONE_SCHEME {
    tag "$scheme_repo_url"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    val scheme
    val scheme_version
    val scheme_repo_url
    
    output:
    path "${scheme_dir}"                                              , emit: scheme
    path "${scheme_dir}/${scheme}/${scheme_version}/*.reference.fasta", emit: fasta
    path "${scheme_dir}/${scheme}/${scheme_version}/*.primer.bed"     , emit: bed

    script:
    def lastPath = scheme_repo_url.lastIndexOf(File.separator)
    def lastExt = scheme_repo_url.lastIndexOf(".")
    scheme_dir = scheme_repo_url.substring(lastPath+1,lastExt)
    """
    git clone $scheme_repo_url $scheme_dir
    """
}