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
if (params.references) { ch_references = file(params.references) } else { exit 1, '[archaeodiet] ERROR: Reference directory not specified!' }
//if (params.input) { ch_input = file(params.db_kraken) } else { exit 1, '[archaeodiet] ERROR: Kraken2 database not specified!' }
//if (params.input) { ch_input = file(params.db_ete) } else { exit 1, '[archaeodiet] ERROR: ETE toolkit database not specified!' }

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

// TODO JAFY MODIFY ALL MODULE ARGS/OPTIONS HERE: https://youtu.be/ggGGhTMgyHI?t=1542

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include { BOWTIE2_MAP } from '../modules/local/bowtie2/map'

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
include { BOWTIE2_BUILD } from '../modules/nf-core/modules/bowtie2/build/main'
include { KRAKEN2_KRAKEN2 } from '../modules/nf-core/modules/kraken2/kraken2/main'
include { SAMTOOLS_MERGE } from '../modules/nf-core/modules/samtools/merge/main'
include { DAMAGEPROFILER } from '../modules/nf-core/modules/damageprofiler/main'

include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )


/*
========================================================================================
    PREP CHANNELS
========================================================================================
*/

Channel.value(file( "${params.references}/*.fasta.gz" )).flatten().dump(tag: "Ref. pickup").set { ch_references }
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

    //
    //SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    INPUT_CHECK (
        ch_input
    )

    BOWTIE2_BUILD (
        ch_references
    )

    // Convert reference channel to id + file
    BOWTIE2_BUILD.out.index
        .dump(tag: "output indices")
        .map{
            def reference = it
            def new_file = reference + '/*.bt2'
            def ind_path = file(new_file)[1].getBaseName()
            def id = ind_path.split("\\.")[0]

            [id, [reference]]
        }
        .dump(tag: "output indices")
        .set { ch_indexed_refs }


    INPUT_CHECK.out.reads.dump(tag: "inputcheck").cross(ch_indexed_refs).set{ ch_mapping_input }

    // TODO ISSUE NEED TO USE .combine OPERATOR TO GET PAIRWISE COMBS OF READS WITH REFERENCES
    // HOW TO SUPPLY TO BT2_MAPPING? .cross()
    // TODO Does not align because I _think_ bt2 is auto detects paired or single-end reads,
    // and takes it from auto metadata propagation in check_samplesheet.py. Need to back and
    // zre-add fastq_2 stuff to check_samplesheet.py AND sample sheet itself
    BOWTIE2_MAP (
        ch_mapping_input
    )

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
    ch_software_versions = ch_software_versions.mix(BOWTIE2_BUILD.out.version.first().ifEmpty(null),
                                                    BOWTIE2_MAP.out.version.first().ifEmpty(null)
                                                    )

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
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

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
