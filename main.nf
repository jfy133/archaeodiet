#!/usr/bin/env nextflow
/*
========================================================================================
 archaeodiet Analysis Pipeline.
 Pipeline built based on the nf-core iniative 
 #### Homepage / Documentation
 https://github.com/jfy133/archaeodiet
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

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
// WORKFLOW: Run main archaeodiet analysis pipeline
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
