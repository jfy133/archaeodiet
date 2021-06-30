/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def valid_params = [
    bt2_alignmode   : ['local', 'end-to-end'],
    bt2_sensitivity : ['no-preset', 'very-fast', 'fast', 'sensitive', 'very-sensitive']
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowArchaeodiet.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.multiqc_config,
    params.references, params.db_kraken, params.db_ete
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, '[archaeodiet] ERROR: Input samplesheet not specified!' }
if (params.input) { ch_input = file(params.references) } else { exit 1, '[archaeodiet] ERROR: Reference directory not specified!' }
if (params.input) { ch_input = file(params.db_kraken) } else { exit 1, '[archaeodiet] ERROR: Kraken2 database not specified!' }
// if (params.input) { ch_input = file(params.db_ete) } else { exit 1, '[archaeodiet] ERROR: ETE toolkit database not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { BOWTIE2_BUILD } from '../modules/nf-core/software/bowtie2/build/main' addParams( options: bowtie2_build_options   )
include { BOWTIE2_ALIGN } from '../modules/nf-core/software/bowtie2/align/main' addParams( options: bowtie2_align_options   )
include { KRAKEN2_RUN } from '../modules/nf-core/software/kraken/run/main' addParams( options: kraken2_run_options   )
include { SAMTOOLS_MERGE } from '../modules/nf-core/software/samtools/merge/main' addParams( options: samtools_merge_options   )
include { DAMAGEPROFILER } from '../modules/nf-core/software/damageprofiler/main' addParams( options: damageprofiler_options )

include { MULTIQC } from '../modules/nf-core/software/multiqc/main' addParams( options: multiqc_options   )


/*
========================================================================================
    PREP CHANNELS
========================================================================================
*/

Channel.value(file( "${params.references}/*.fna.gz" )).set { ch_references }
Channel.value(file( "${params.db_kraken}" )).set { ch_db_kraken }
// ch_db_ete = Channel.value(file( "${params.db_ete}" ))

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ARCHAEODIET {

    ch_software_versions = Channel.empty()


    SUBWORKFLOW: Read in samplesheet, validate and stage input files

    INPUT_CHECK (
        ch_input
    )

    BOWTIE2_INDEXING (
        ch_references
    )
    ch_software_versions = ch_software_versions.mix(BOWTIE2_INDEXING.out.bowtie2_version.first().ifEmpty(null))

    // ISSUE NEED TO USE .combine OPERATOR TO GET PAIRWISE COMBS OF READS WITH REFERENCES
    // HOW TO SUPPLY TO BT2_MAPPING?

    // BOWTIE2_MAPPING (
    //     INPUTCHECK.out, ch_references
    // )

    // KRAKEN_CLEANING (

    // )


    //
    // MODULE: Run FastQC
    //
    // FASTQC (
    //     INPUT_CHECK.out.reads
    // )
    // ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowArchaeodiet.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
