// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process BOWTIE2_MAP {
    tag "${meta_reads.id}-${meta_ref.id}"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:$prefix, publish_by_meta:$prefix) }

    conda (params.enable_conda ? "bioconda::bowtie2=2.4.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0"
    }

    input:
    tuple val(meta_reads), path(reads), val(meta_ref), path(reference)

    output:
    tuple val(meta_reads), val(meta_ref), path("*.bam"), emit: bam
    tuple val(meta_reads), val(meta_ref), path('*.log'), emit: log
    path "*.version.txt"          , emit: version

    script:
    def split_cpus = Math.floor(task.cpus/2)
    def software = getSoftwareName(task.process)
    def prefix   = "${meta_reads.id}-${meta_ref.id}"

    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`

    ## by default unaligned reads only to keep BAMs small!
    bowtie2 \\
        -x \$INDEX \\
        -U $reads \\
        --threads ${split_cpus} \\
        $options.args \\
        2> ${prefix}.bowtie2.log \\
        | samtools view -@ ${split_cpus} -F 4 $options.args2 -bhS -o ${prefix}.bam -

    echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' > ${software}.version.txt
    """
}
