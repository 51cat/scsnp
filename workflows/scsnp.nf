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

//

process TARGET_METRICS {
    tag "TARGET_METRICS"
    label "snp_target"

    // conda
    // container 

    input:
    path        bam 
    val         meta
    path        outdir
    path        match_dir
    path        gene_list

    output:
    path "$outdir/${sample}_filtered_sorted.bam" ,emit: target_bam

    script:
    """
    target_metrics.py \\
    --input_bam $bam \\
    --sample $sample \\
    --outdir $outdir \\
    --gene_list $gene_list \\
    --match_dir $match_dir  \\
    --add_RG
    """
}

process CALLING_PREPROCESS {
    tag "CALLING_PREPROCESS"
    label "snp_calling_preprocess"

    // conda
    // container 

    input:
    path        bam 
    val         sample
    path        outdir
    path        fasta

    output:
    path "$outdir/${sample}_splitN.bam" ,emit: exon_bam

    """
    calling_preprocess.py \
    --input_bam $bam\
    --sample $sample\
    --outdir $outdir\
    --fasta $fasta
    """
}

process CALLING {
    tag "CALLING"
    label "snp_calling"

    // conda
    // container 

    input:
    path        bam 
    val         sample
    path        outdir
    path        fasta
    path        bed_file
    val         thread

    output:
    path "$outdir/${sample}_raw.bcf" ,emit: raw_bcf
    path "$outdir/${sample}_raw.vcf" ,emit: raw_vcf
    path "$outdir/${sample}_fixed.vcf" ,emit: fixed_vcf
    path "$outdir/${sample}_norm.vcf" ,emit: norm_vcf

    """
    calling.py\
    --input_bam $bam\
    --sample $sample\
    --outdir $outdir\
    --fasta $fasta \
    --thread $thread
    """
}

process FILTER_SNP {
    tag "FILTER_SNP"
    label "snp_filter"

    // conda
    // container 

    input:
    path        vcf 
    val         sample
    path        outdir
    val         ref_threshold_method
    val         alt_threshold_method
    val         VAF
    val         ref_min_support_read
    val         alt_min_support_read

    output:
    path "$outdir/${sample}_filtered.vcf" ,emit: filter_vcf

    """
    filter_snp.py\
    --vcf $vcf\
    --sample $sample\
    --outdir $outdir \\
    --ref_threshold_method $ref_threshold_method \\
    --alt_threshold_method $alt_threshold_method \\
    --VAF $VAF \\
    --ref_min_support_read $ref_min_support_read \\
    --alt_min_support_read $alt_min_support_read
    """

}

workflow SCSNP {
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
