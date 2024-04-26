/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { MULTIQC                } from '../modules/local/multiqc_sgr/main'
// include { paramsSummaryMap       } from 'plugin/nf-validation'
// include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
// include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_scrna_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// test

params.bam = ""
params.sample = ""
params.outdir = ""
params.match_dir = ""
params.gene_list = ""

process TARGET_METRICS {
    tag "TARGET_METRICS"
    label "snp_target"

    input:
    path bam
    val sample
    path outdir
    path match_dir
    path gene_list

    output:
    path "$outdir/${sample}_filtered_sorted.bam"

    script:
    """
    python /SGRNJ06/randd/USER/liuzihao/work/scsnp/bin/target_metrics.py \\
    --input_bam $bam \\
    --sample $sample \\
    --outdir $outdir \\
    --gene_list $gene_list \\
    --match_dir $match_dir  \\
    --add_RG
    """
}

workflow {

    def ch_bam = Channel.fromPath(params.bam)
    def ch_outdir = Channel.fromPath(params.outdir)
    def ch_match_dir = Channel.fromPath(params.match_dir)
    def ch_gene_list = Channel.fromPath(params.gene_list)

    TARGET_METRICS(
        ch_bam, params.sample, ch_outdir, ch_match_dir, ch_gene_list
    )

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
