#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/archaeodiet
========================================================================================
 nf-core/archaeodiet Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/archaeodiet
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/archaeodiet --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --input [file]                Path to preprocessed FASTQ files (must be surrounded with quotes)
      -profile [str]                Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, test, awsbatch, <institute> and more

    References                        If not specified in the configuration file or you wish to overwrite any of the references
      --target_db [file/dir]          Path to target metagenomic screening reference (e.g. eukaryotic dietary species)
      --contaminant_db [file/dir]     Path to contaminants to screen against (e.g. microbial database)
      --ete3toolkit_db                         Path to ete3 toolkit taxa.sqlite database, if not in ~/.etetoolkit/

    Target Screening
      --target_taxonomic_level                Specify at which taxonomic level to screen for target taxa for. Options: 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'. Default: ${params.target_taxonomic_level}
      --target_tool                           Specify which classifier to use. Options: 'malt'. Default: '${params.target_tool}'
      --target_min_support_reads              Specify a minimum number of reads  a taxon of sample total is required to have to be retained. Not compatible with . Default: ${params.target_min_support_reads}
      --target_percent_identity               Percent identity value threshold for MALT. Default: ${params.target_percent_identity}
      --target_malt_mode                      Specify which alignment method to use for MALT. Options: 'Unknown', 'BlastN', 'BlastP', 'BlastX', 'Classifier'. Default: '${params.target_malt_mode}'
      --target_malt_alignment_mode            Specify alignment method for MALT. Options: 'Local', 'SemiGlobal'. Default: '${params.target_malt_alignment_mode}'
      --target_malt_top_percent               Specify the percent for LCA algorithm for MALT (see MEGAN6 CE manual). Default: ${params.target_malt_top_percent}
      --target_malt_min_support_mode          Specify whether to use percent or raw number of reads for minimum support required for taxon to be retained for MALT. Options: 'percent', 'reads'. Default: '${params.target_malt_min_support_mode}'
      --target_malt_min_support_percent       Specify the minimum percentage of reads a taxon of sample total is required to have to be retained for MALT. Default: Default: ${params.target_malt_min_support_percent}
      --target_malt_max_queries               Specify the maximium number of queries a read can have for MALT. Default: ${params.target_malt_max_queries}
      --target_malt_memory_mode               Specify the memory load method. Do not use 'map' with GTFS file system for MALT. Options: 'load', 'page', 'map'. Default: '${params.target_malt_memory_mode}'
      --target_malt_weighted_lca              Specify whether to use MALT's 'weighted' LCA algorithm
      --target_malt_weighted_lca_perc         Specify the weighted-LCA percentage of weight to cover. Default: ${params.target_malt_weighted_lca_perc}

    Contaminant Screening
      --contaminant_tool                           Specify which classifier to use. Options: 'malt', 'kraken'. Default: '${params.contaminant_tool}'
      --contaminant_min_support_reads              Specify a minimum number of reads  a taxon of sample total is required to have to be retained. Not compatible with . Default: ${params.contaminant_min_support_reads}
      --contaminant_percent_identity               Percent identity value threshold for MALT. Default: ${params.contaminant_percent_identity}
      --contaminant_malt_mode                      Specify which alignment method to use for MALT. Options: 'Unknown', 'BlastN', 'BlastP', 'BlastX', 'Classifier'. Default: '${params.contaminant_malt_mode}'
      --contaminant_malt_alignment_mode            Specify alignment method for MALT. Options: 'Local', 'SemiGlobal'. Default: '${params.contaminant_malt_alignment_mode}'
      --contaminant_malt_top_percent               Specify the percent for LCA algorithm for MALT (see MEGAN6 CE manual). Default: ${params.contaminant_malt_top_percent}
      --contaminant_malt_min_support_mode          Specify whether to use percent or raw number of reads for minimum support required for taxon to be retained for MALT. Options: 'percent', 'reads'. Default: '${params.contaminant_malt_min_support_mode}'
      --contaminant_malt_min_support_percent       Specify the minimum percentage of reads a taxon of sample total is required to have to be retained for MALT. Default: Default: ${params.contaminant_malt_min_support_percent}
      --contaminant_malt_max_queries               Specify the maximium number of queries a read can have for MALT. Default: ${params.contaminant_malt_max_queries}
      --contaminant_malt_memory_mode               Specify the memory load method. Do not use 'map' with GTFS file system for MALT. Options: 'load', 'page', 'map'. Default: '${params.contaminant_malt_memory_mode}'
      --contaminant_malt_weighted_lca              Specify whether to use MALT's 'weighted' LCA algorihthm
      --contaminant_malt_weighted_lca_per          Specify the weighted-LCA percentage of weight to cover. Default: ${params.contaminant_malt_weighted_lca_perc}


    Damage Profiling
      --damageprofiler_length       Specify length filter for DamageProfiler. Default: ${params.damageprofiler_length}
      --damageprofiler_threshold    Specify number of bases to consider for damageProfiler (e.g. on damage plot). Default: ${params.damageprofiler_threshold}
      --damageprofiler_yaxis        Specify the maximum misincorporation frequency that should be displayed on damage plot. Set to 0 to 'autoscale'. Default: ${params.damageprofiler_yaxis} 
      --pydamage_windowlength       Specify length of read to perform model fitting
      --pydamage_minreads           Specify minimum number of reads required for calculation
      --pydamage_coverage           Specify minimum coverage required for calculation
      --pydamage_plots              Specify whether to produce pydamage plot images
      --pydamage_alignments         Specify whether to produce pydamage alignment representations

    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      --max_multiqc_email_size [str]  Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

//log.info("workflow params")
//log.info("${params}")

/*
 * SET UP CONFIGURATION VARIABLES
 */

if (params.target_db == '' ) {
    exit 1, "[nf-core/archaeodiet] error: target database alignment requires a path to a database directory. Please specify one with --target_db '/path/to/database/'."
}

if (params.contaminant_db == '' ) {
    exit 1, "[nf-core/archaeodiet] error: contaminant database alignment requires a path to a database directory. Please specify one with --contaminant_db '/path/to/database/'."
}

// Input validation
def valid_tax_levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
if ( !params.target_taxonomic_level in valid_tax_levels) {
    exit 1, "[nf-core/archaeodiet] error: invalid taxonomic level specificed. Options: ${valid_tax_levels}. You gave: --target_taxonomic_level '${target_taxonomic_level}'."
}

  if (params.target_tool != 'malt') {
    exit 1, "[nf-core/archaeodiet] error: metagenomic classification against target database can currently only be run with 'malt'. Please check your classifer. You gave: --target_tool '${params.target_tool}'."
  }

  if (params.target_tool == 'malt' && params.target_malt_mode != 'BlastN' && params.target_malt_mode != 'BlastP' && params.target_malt_mode != 'BlastX') {
    exit 1, "[nf-core/archaeodiet] error: unknown MALT mode specified. Options: 'BlastN', 'BlastP', 'BlastX'. You gave: --target_malt_mode '${params.target_malt_mode}'."
  }

  if (params.target_tool == 'malt' && params.target_malt_alignment_mode != 'Local' && params.target_malt_alignment_mode != 'SemiGlobal') {
    exit 1, "[nf-core/archaeodiet] error: unknown MALT alignment mode specified. Options: 'Local', 'SemiGlobal'. You gave: --target_malt_alignment_mode '${params.target_malt_alignment_mode}'."
  }

  if (params.target_tool == 'malt' && params.target_malt_min_support_mode == 'percent' && params.target_min_support_reads != 1) {
    exit 1, "[nf-core/archaeodiet] error: incompatible MALT min support configuration. Percent can only be used with --target_malt_min_support_percent. You modified --target_min_support_reads."
  }

  if (params.target_tool == 'malt' && params.target_malt_min_support_mode == 'reads' && params.target_malt_min_support_percent != 0.01) {
    exit 1, "[nf-core/archaeodiet] error: incompatible MALT min support configuration. Reads can only be used with --target_target_malt_min_supportreads. You modified --target_malt_min_support_percent."
  }

  if (params.target_tool == 'malt' && params.target_malt_memory_mode != 'load' && params.target_malt_memory_mode != 'page' && params.target_malt_memory_mode != 'map') {
    exit 1, "[nf-core/archaeodiet] error: unknown MALT memory mode specified. Options: 'load', 'page', 'map'. You gave: --target_malt_memory_mode '${params.target_malt_memory_mode}'."
  }

  if (!params.target_min_support_reads.toString().isInteger()){
    exit 1, "[nf-core/archaeodiet] error: incompatible min_support_reads configuration. min_support_reads can only be used with integers. --target_min_support_reads You gave: ${params.target_min_support_reads}."
  }

    if (params.contaminant_tool != 'malt') {
    exit 1, "[nf-core/archaeodiet] error: metagenomic classification against contaminant database can currently only be run with 'malt'. Please check your classifer. You gave: --target_tool '${params.target_tool}'."
  }

  if (params.contaminant_tool == 'malt' && params.contaminant_malt_mode != 'BlastN' && params.contaminant_malt_mode != 'BlastP' && params.contaminant_malt_mode != 'BlastX') {
    exit 1, "[nf-core/archaeodiet] error: unknown MALT mode specified. Options: 'BlastN', 'BlastP', 'BlastX'. You gave: --contaminant_malt_mode '${params.contaminant_malt_mode}'."
  }

  if (params.contaminant_tool == 'malt' && params.contaminant_malt_alignment_mode != 'Local' && params.contaminant_malt_alignment_mode != 'SemiGlobal') {
    exit 1, "[nf-core/archaeodiet] error: unknown MALT alignment mode specified. Options: 'Local', 'SemiGlobal'. You gave: --contaminant_malt_alignment_mode '${params.contaminant_malt_alignment_mode}'."
  }

  if (params.contaminant_tool == 'malt' && params.contaminant_malt_min_support_mode == 'percent' && params.contaminant_min_support_reads != 1) {
    exit 1, "[nf-core/archaeodiet] error: incompatible MALT min support configuration. Percent can only be used with --contaminant_malt_min_support_percent. You modified --contaminant_min_support_reads."
  }

  if (params.contaminant_tool == 'malt' && params.contaminant_malt_min_support_mode == 'reads' && params.contaminant_malt_min_support_percent != 0.01) {
    exit 1, "[nf-core/archaeodiet] error: incompatible MALT min support configuration. Reads can only be used with --contaminant_contaminant_malt_min_supportreads. You modified --contaminant_malt_min_support_percent."
  }

  if (params.contaminant_tool == 'malt' && params.contaminant_malt_memory_mode != 'load' && params.contaminant_malt_memory_mode != 'page' && params.contaminant_malt_memory_mode != 'map') {
    exit 1, "[nf-core/archaeodiet] error: unknown MALT memory mode specified. Options: 'load', 'page', 'map'. You gave: --contaminant_malt_memory_mode '${params.contaminant_malt_memory_mode}'."
  }

  if (!params.contaminant_min_support_reads.toString().isInteger()){
    exit 1, "[nf-core/archaeodiet] error: incompatible min_support_reads configuration. min_support_reads can only be used with integers. --target_min_support_reads You gave: ${params.target_min_support_reads}."
  }


// Has the run name been specified by the user? this has the bonus effect of
// catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch related:
    // https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling
    // files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

/*
 * Input Loading
 */

 // TODO Fix shitty splitting at '.'
if (params.input) {
    Channel
        .fromPath(params.input, checkIfExists: true)
        .ifEmpty { exit 1, "[nf-core/archaeodiet] Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!" }
        .set { ch_input_for_targetalignment }
} else {
    exit 1, "[nf-core/archaeodiet] input was not supplied. PLease check input parameters"
}

if ( params.target_db != '' ) {
  ch_db_for_targetalignment = Channel
      .fromPath(params.target_db, checkIfExists: true, type: 'dir')
} else {
    exit 1, "[nf-core/archaeodiet] target database was not supplied. Please check input parameters."
}

if ( params.contaminant_db != '' ) {
  ch_db_for_contaminantalignment = Channel
      .fromPath(params.contaminant_db, checkIfExists: true, type: 'dir')
} else {
    exit 1, "[nf-core/archaeodiet] contaminant database was not supplied. Please check input parameters."
}

// Load dummy header for pydamage output
pydamage_header = file("$baseDir/assets/multiqc_pydamage_customcontent_header.txt")

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Input']            = params.input
summary['Target DB']        = params.target_db
summary['Contaminant DB']   = params.contaminant_db
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-archaeodiet-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/archaeodiet Workflow Summary'
    section_href: 'https://github.com/nf-core/archaeodiet'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */

process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo "${workflow.manifest.version}" > v_pipeline.txt
    echo "${workflow.nextflow.version}" > v_nextflow.txt

    malt-run --help \
    |& tail -n 3 \
    | head -n 1 \
    | cut -f 2 -d'(' \
    | cut -f 1 -d ',' &> v_malt.txt 2>&1 || true
    
    samtools --version &> v_samtools.txt 2>&1 || true   
    damageprofiler --version &> v_damageprofiler.txt 2>&1 || true
    picard FilterSamReads --version &> v_filtersamreads.txt || true
    pydamage --version | cut -d " " -f 3 &> v_pydamage.txt || true
    collapse_sam_taxonomy.py --version &> vcollapse_sam_taxonomy.txt || true
    
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * STEP 1- Target DB Screening
 */

process target_alignment_malt {
  label 'mc_small'
  publishDir "${params.outdir}/target_alignment", 
    mode: params.publish_dir_mode, 
    pattern: '*.log'

  when:
  params.target_tool == 'malt'

  input:
  path(fastqs) from ch_input_for_targetalignment.collect()
  path db from ch_db_for_targetalignment

  output:
  path "*.sam.gz" into ch_sam_for_targetsamtobam mode flatten
  path "target_malt.log" into ch_targetmalt_for_multiqc

  script:
  if ( "${params.target_malt_min_support_mode}" == "percent" ) {
    min_supp = "-supp ${params.target_malt_min_support_percent}" 
  } else if ( "${params.target_malt_min_support_mode}" == "reads" ) {
    min_supp = "-sup ${params.target_min_support_reads}"
  }
    
  if ( params.target_malt_weighted_lca ) {
    wlca = "-wLCA -wLCAP ${params.target_malt_weighted_lca_perc}"
  } else {
    wlca = ""
  } 
  """
  malt-run \
  -J-Xmx${task.memory.toGiga()}g \
  -t ${task.cpus} \
  -v \
  -o . \
  --alignments ./ \
  -d ${db} \
  -id ${params.target_percent_identity} \
  -m ${params.target_malt_mode} \
  -at ${params.target_malt_alignment_mode} \
  -top ${params.target_malt_top_percent} \
  ${min_supp} \
  -mq ${params.target_malt_max_queries} \
  --memoryMode ${params.target_malt_memory_mode} \
  ${wlca} -i ${fastqs.join(' ')} |&tee target_malt.log

  rm *.rma6
  """
}

process target_fileconversion {
  label 'mc_small'
  tag "${sam}"
  publishDir "${params.outdir}/target_alignment", mode:"copy"

  input:
  path(sam) from ch_sam_for_targetsamtobam

  output:
  path "*bam" into ch_bam_for_contaminationalignment, ch_bam_for_contaminantremoval
  path "*.fastq.gz" into ch_fastq_for_contaminantalignment

  script:
  """
  zcat $sam | samtools view -@ ${task.cpus} -b -o ${sam}_putativehits.bam
  samtools fastq \
  -@ ${task.cpus} \
  ${sam}_putativehits.bam \
  | pigz -p ${task.cpus} > ${sam}_putativehits.fastq.gz
  """
}

process contaminant_alignment_malt {
  label 'mc_small'

  publishDir "${params.outdir}/contaminant_alignment", 
    mode: params.publish_dir_mode, 
    pattern: '*.log'

  when:
  params.target_tool == 'malt'

  input:
  path(fastqs) from ch_fastq_for_contaminantalignment.collect()
  path db from ch_db_for_contaminantalignment

  output:
  path "*.sam.gz" into ch_sam_for_contaminantsamtobam mode flatten
  path "contaminant_malt.log" into ch_contaminantmalt_for_multiqc

  script:
  if ( "${params.contaminant_malt_min_support_mode}" == "percent" ) {
    min_supp = "-supp ${params.contaminant_malt_min_support_percent}" 
  } else if ( "${params.contaminant_malt_min_support_mode}" == "reads" ) {
    min_supp = "-sup ${params.contaminant_min_support_reads}"
  }
  if ( params.contaminant_malt_weighted_lca ) {
    wlca = "-wLCA -wLCAP ${params.contaminant_malt_weighted_lca_perc}"
  } else {
    wlca = ""
  }
  """
  malt-run \
  -J-Xmx${task.memory.toGiga()}g \
  -t ${task.cpus} \
  -v \
  -o . \
  --alignments ./ \
  -d ${db} \
  -id ${params.contaminant_percent_identity} \
  -m ${params.contaminant_malt_mode} \
  -at ${params.contaminant_malt_alignment_mode} \
  -top ${params.contaminant_malt_top_percent} \
  ${min_supp} \
  -mq ${params.contaminant_malt_max_queries} \
  --memoryMode ${params.contaminant_malt_memory_mode} \
  ${wlca} -i ${fastqs.join(' ')} |&tee contaminant_malt.log

  rm *.rma6
  """
}

process contaminant_fileconversion {
  label 'mc_small'
  tag "${sam}"  

  publishDir "${params.outdir}/contaminant_alignment", 
    mode: params.publish_dir_mode, 
    pattern: '*.bam'
  
  input:
  path(sam) from ch_sam_for_contaminantsamtobam

  output:
  path "*bam" into ch_bam_for_contaminantreadextraction

  script:
  """
  zcat $sam | samtools view -@ ${task.cpus} -b -o ${sam}_putativecontaminants.bam
  """
}

process contaminant_readextraction {
  label 'mc_small'
  publishDir "${params.outdir}/contaminant_alignment", 
    mode: params.publish_dir_mode, 
    pattern: '*.txt.gz'

  tag "${bam}"  

  input:
  path(bam) from ch_bam_for_contaminantreadextraction

  output:
  path "*txt.gz" into ch_readlist_for_contaminantremoval

  script:
  """
  samtools view \
  -@ ${task.cpus} \
  ${bam} \
  | cut -f 1 \
  | pigz -p ${task.cpus} > ${bam}_putatitivecontaminants_readids.txt.gz
  """
}

process contaminant_removal {
  label 'mc_small'
  tag "${bam}"
  publishDir "${params.outdir}/target_alignment_cleaned", 
    mode: params.publish_dir_mode, 
    pattern: '*_cleaned.bam'

  input:
  path(bam) from ch_bam_for_contaminantremoval
  path(txt) from ch_readlist_for_contaminantremoval

  output:
  path "*_cleaned.bam" into ch_bam_for_damageprep

  """
  picard FilterSamReads READ_LIST_FILE=${txt} FILTER=excludeReadList I=${bam} O=${bam}_cleaned.bam
  """
}

process target_taxonomy_collapsing {
  label 'mc_small'
  tag "${bam}"

  input:
  path(bam) from ch_bam_for_damageprep

  output:
  tuple path("*_tophits.bam"), path("results/*.txt") into ch_bam_for_damageauthentication

  script:
  """
  samtools view -b -F 256 ${bam} | samtools sort -@ ${task.cpus} > ${bam}_tophits.bam
  collapse_sam_taxonomy.py -i ${bam}_tophits.bam -t ${params.target_taxonomic_level} ${params.ete3toolkit_db}
  """
}

process target_bam_splitting {
  label 'mc_small'
  tag "${bam}"

  input:
  tuple path(bam), path(taxids) from ch_bam_for_damageauthentication

  output:
  path("results/*.bam") into ch_bam_for_damageprofiler
  path("results/*.bam") into ch_bam_for_pydamage

  script:
  samplename = bam.baseName
  """
  samtools index ${bam}
  mkdir results/
  rm unresolved_taxonomic_ids.txt
  for i in *.txt; do
    taxonname=\$(echo "\$i" | rev | cut -d '.' -f 2-999999999 | rev)
    samtools view -b ${bam} \$(cat \$i) | samtools sort -@ ${task.cpus} > results/${samplename}_"\$taxonname".bam
  done
  """
}

process target_damageprofiler {
  label 'mc_small'
  tag "${bam}"
  publishDir "${params.outdir}/damageprofiler/", 
    mode: params.publish_dir_mode

  input:
  path(bam) from ch_bam_for_damageprofiler.flatten()

  output:
  file "*/*.txt"
  file "*/*.log"
  file "*/*.pdf"
  file "*/*json" into ch_damageprofiler_for_multiqc

  script:
  samplename = bam.baseName
  """
  damageprofiler -i ${bam} -o . -yaxis_damageplot 0.30
  """
}

process target_pydamage {
  label 'mc_small'
  tag "${bam}"
  publishDir "${params.outdir}/pydamage/", 
    mode: params.publish_dir_mode

  input:
  file(bam) from ch_bam_for_pydamage.flatten()

  output:
  path "${samplename}/pydamage_results.csv" optional true
  path "${samplename}/plots" optional true
  path("*/*_pydamage_mqc.csv") optional true into ch_pydamage_for_multiqc

  script:
  samplename = bam.baseName
  def alignments = params.pydamage_alignments ? "-s" : ""
  def plots = params.pydamage_plots ? "-pl" : ""
  """
  samtools index ${bam}
  pydamage -o ${samplename} -w ${params.pydamage_windowlength} -m ${params.pydamage_minreads} -c ${params.pydamage_coverage} ${plots} ${alignments} -p ${task.cpus} ${bam}
  
  if [[ -e  ${samplename}/pydamage_results.csv ]]; then
    ## make new sample column
    ref_len=\$(wc -l ${samplename}/pydamage_results.csv | cut -f 1 -d ' ')
    samples=\$(expr \$ref_len - 1)
    filename=\$(dirname ${samplename}/pydamage_results.csv)

    paste <(echo "sample_name" && yes "\$filename" | head -n \${samples}) ${samplename}/pydamage_results.csv  -d ',' > ${samplename}/tmp.csv
    ## merge header
    cat ${pydamage_header} ${samplename}/tmp.csv > ${samplename}/${samplename}_pydamage_mqc.csv
  fi
  """
}

/*
 * STEP 2 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: params.publish_dir_mode

    input:
    file (multiqc_config) from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    // TODO nf-core: Add in log files from your new processes for MultiQC to
    // find!
    file ('software_versions/*') from ch_software_versions_yaml.collect()
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
    file ('target_alignment/*') from ch_targetmalt_for_multiqc.collect().ifEmpty([])
    file ('contaminant_alignment/*') from ch_contaminantmalt_for_multiqc.collect().ifEmpty([])
    file ('damageprofiler*/') from ch_damageprofiler_for_multiqc.collect().ifEmpty([])
    file ('pydamage/*') from ch_pydamage_for_multiqc.collect().ifEmpty([])

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f $rtitle $rfilename $multiqc_config $custom_config_file . -s ## add s as test TODO to remove and fix properly
    """
}

/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/archaeodiet] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/archaeodiet] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/archaeodiet] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/archaeodiet] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/archaeodiet] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, email_address ].execute() << email_txt
            log.info "[nf-core/archaeodiet] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/archaeodiet]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/archaeodiet]${c_red} Pipeline completed with errors${c_reset}-"
    }

}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/archaeodiet v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
