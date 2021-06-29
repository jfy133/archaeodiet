#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/archaeodiet
========================================================================================
    Github : https://github.com/nf-core/archaeodiet
    Website: https://nf-co.re/archaeodiet
    Slack  : https://nfcore.slack.com/channels/archaeodiet
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

//log.info("workflow params")
//log.info("${params}")

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { ARCHAEODIET } from './workflows/archaeodiet'

//
// WORKFLOW: Run main nf-core/archaeodiet analysis pipeline
//
workflow NFCORE_ARCHAEODIET {
    ARCHAEODIET ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_ARCHAEODIET ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
