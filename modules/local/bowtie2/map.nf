// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process BOWTIE2_MAP {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bowtie2=2.4.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/https://depot.galaxyproject.org/singularity/bowtie2:2.4.4--py39hbb4e92a_0"
    } else {
        container "quay.io/biocontainers/quay.io/biocontainers/bowtie2:2.4.4--py36hd4290be_0"
    }

    input:
    tuple val(meta), path(bam), path(reference)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def ref_prefix = ${meta.ref}

    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
    bowtie2 \\
        -x \$INDEX \\
        -U $reads \\
        --threads ${split_cpus} \\
        $unaligned \\
        $options.args \\
        2> ${prefix}.bowtie2.log \\
        | samtools view -@ ${split_cpus} $options.args2 -bhS -o ${prefix}.bam -

    echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' > ${software}.version.txt
    """
}
