#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    singleron-RD/scrna
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/singleron-RD/scrna
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

 include { SCSNP  } from './workflows/scsnp'
 include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_scrna_pipeline'
 include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_scrna_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow SINGLERONRD_SCSNP {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    SCSNP(
        samplesheet
    )

    //
    // WORKFLOW: Run pipeline
    //
    
    //emit:
    //multiqc_report = SCSNP.out.multiqc_report // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    main:

//    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
     PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )
    //PIPELINE_INITIALISATION.out.samplesheet.view()
    SINGLERONRD_SCSNP (
        PIPELINE_INITIALISATION.out.samplesheet
    )
     PIPELINE_INITIALISATION.out.samplesheet.view()

    //
    // WORKFLOW: Run main workflow
    //
    //
    // SUBWORKFLOW: Run completion tasks
    //
//    PIPELINE_COMPLETION (
//        params.email,
//        params.email_on_fail,
 //       params.plaintext_email,
 //       params.outdir,
 //       params.monochrome_logs,
 //       params.hook_url,
 //       SINGLERONRD_SCSNP.out.multiqc_report
 //   )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
