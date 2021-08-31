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
if (params.db_kraken) { ch_db_kraken = file(params.db_kraken) } else { exit 1, '[archaeodiet] ERROR: Kraken2 database not specified!' }
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

// TODO JAFY MODIFY ALL MODULE ARGS/OPTIONS HERE

// Load defaults from config
def bowtie2_map_options = modules['bowtie2']

// Modify defaults
bowtie2_map_options.args += params.bt2_n != 0 ? "-N ${params.bt2_n} " : ""
bowtie2_map_options.args += params.bt2_l != 0 ? "-L ${params.bt2_L} " : ""


if ( "${params.bt2_alignmode}" == "end-to-end"  ) {
    switch ( "${params.bt2_sensitivity}" ) {
        case "no-preset":
        sensitivity = ""; break
        case "very-fast":
        sensitivity = "--very-fast"; break
        case "fast":
        sensitivity = "--fast"; break
        case "sensitive":
        sensitivity = "--sensitive"; break
        case "very-sensitive":
        sensitivity = "--very-sensitive"; break
        default:
        sensitivity = ""; break
        }
} else if ("${params.bt2_alignmode}" == "local") {
    switch ( "${params.bt2_sensitivity}" ) {
        case "no-preset":
        sensitivity = ""; break
        case "very-fast":
        sensitivity = "--very-fast-local"; break
        case "fast":
        sensitivity = "--fast-local"; break
        case "sensitive":
        sensitivity = "--sensitive-local"; break
        case "very-sensitive":
        sensitivity = "--very-sensitive-local"; break
        default:
        sensitivity = ""; break
        }
}

bowtie2_map_options.args += sensitivity

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include { BOWTIE2_MAP } from '../modules/local/bowtie2/map' addParams( options: bowtie2_map_options )
include { EXTRACTID } from '../modules/local/extractid'
include { EXTRACTBAMHEADER } from '../modules/local/extractbamheader'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )
include { BAM_SORT_SAMTOOLS } from '../subworkflows/local/bam_sort_samtools' addParams( options: [:] )  addParams( options: [suffix : '.sorted'  ]   )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

// TODO JAFY MODIFY ALL MODULE ARGS/OPTIONS HERE
def kraken2_options = modules['kraken2']

kraken2_options.args += params.kraken2_confidence != 0.0 ? "--confidence ${params.kraken2_confidence}" : ""
kraken2_options.args += params.kraken2_min_base_qual != 0 ? "--minimum_base_quality ${params.kraken2_minimum_min_base_qual}" : ""
kraken2_options.args += params.kraken2_memory_mapping ? "--memory-mapping" : ""
kraken2_options.args += params.kraken2_min_hit_groups != 2 ? "--minimum-hit-groups ${params.kraken2_min_hit_groups}" : ""

//
// MODULE: Installed directly from nf-core/modules
//
include { BOWTIE2_BUILD } from '../modules/nf-core/modules/bowtie2/build/main'
include { SAMTOOLS_FLAGSTAT } from '../modules/nf-core/modules/samtools/flagstat/main'
include { SAMTOOLS_MERGE } from '../modules/nf-core/modules/samtools/merge/main'
include { SAMTOOLS_FASTQ } from '../modules/nf-core/modules/samtools/fastq/main'
include { UNTAR }  from '../modules/nf-core/modules/untar/main'
include { KRAKEN2_KRAKEN2 } from '../modules/nf-core/modules/kraken2/kraken2/main'
include { PICARD_FILTERSAMREADS } from '../modules/nf-core/modules/picard/filtersamreads/main'  addParams( options: [suffix : '.filtered'  ]   )
include { DAMAGEPROFILER } from '../modules/nf-core/modules/damageprofiler/main'

include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )

/*
========================================================================================
    PREP CHANNELS
========================================================================================
*/

Channel.value(file( "${params.references}/*.fasta.gz" )).flatten().set { ch_references }
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
        .map{
            def meta_ref = [:]
            def reference = it
            def new_file = reference + '/*.bt2'
            def ind_path = file(new_file)[1].getBaseName()
            meta_ref.id = ind_path.split("\\.")[0]

            [meta_ref, [reference]]
        }
        .set { ch_indexed_refs }

    // Generate pairwise combination all reads and refs
    INPUT_CHECK.out.reads
        .combine(ch_indexed_refs)
        .dump(tag: "input_check_combine")
        .set{ ch_mapping_input }

    BOWTIE2_MAP (
        ch_mapping_input
    )

    // TODO FIX HERE
    BOWTIE2_MAP.out.bam
        .dump(tag: "bt2_out")
        .map {

            def meta = [:]
            def meta_reads = it[0]
            def meta_ref = it[1]
            def bam = it[2]

            meta.id = meta_reads.id
            meta.ref = meta_ref.id
            meta.longname = "${meta.id}-${meta.ref}"

            [ meta, bam ]
        }
        .set { ch_bt2_for_flagstat }

    BAM_SORT_SAMTOOLS (
        ch_bt2_for_flagstat.dump(tag: "input to stats")
    )

    // Combine for merging
    BOWTIE2_MAP.out.bam
        .map {
            def meta = [:]
            meta.id = it[0].id
            meta.single_end = true
            bams = it[2]

            array = [ meta, bams ]

            return array
        }
        .groupTuple()
        .dump(tag: "input_to_merge")
        .set{ch_bam_for_merge}

    SAMTOOLS_MERGE (
        ch_bam_for_merge
    )

    SAMTOOLS_FASTQ (
        SAMTOOLS_MERGE.out.bam
    )

    // TODO subworkflow to allow either tar/gzipped OR uncompressed folder
    UNTAR (
        ch_db_kraken
    )

    KRAKEN2_KRAKEN2 (
        SAMTOOLS_FASTQ.out.fastq, UNTAR.out.untar
    )

    EXTRACTID (
        KRAKEN2_KRAKEN2.out.unclassified
    )

    SAMTOOLS_MERGE.out.bam
        .combine(EXTRACTID.out.readlist, by: 0)
        .set { ch_input_for_filtersamreads }

    PICARD_FILTERSAMREADS (
        ch_input_for_filtersamreads, 'excludeReadList'
    )

    // TODO replace with bamAlignerHeader
    EXTRACTBAMHEADER (
        PICARD_FILTERSAMREADS.out.bam
    )

    SAMTOOLS_MERGE.out.bam
        .combine( EXTRACTBAMHEADER.out.reflist, by: 0)
        .set { ch_input_for_damageprofiling }

    //DAMAGEPROFILER (
    //
    //)

    // TODO
    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions = ch_software_versions.mix(BOWTIE2_BUILD.out.version.first().ifEmpty(null),
                                                    BOWTIE2_MAP.out.version.first().ifEmpty(null),
                                                    SAMTOOLS_MERGE.out.version.first().ifEmpty(null),
                                                    SAMTOOLS_FASTQ.out.version.first().ifEmpty(null),
                                                    UNTAR.out.version.first().ifEmpty(null),
                                                    KRAKEN2_KRAKEN2.out.version.first().ifEmpty(null),
                                                    EXTRACTID.out.version.first().ifEmpty(null),
                                                    PICARD_FILTERSAMREADS.out.version.first().ifEmpty(null),
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
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_MAP.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_SAMTOOLS.out.idxstats.dump(tag: "idxstats").collect{it[1]}.ifEmpty([]))

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
