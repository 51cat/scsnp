process CUTADAPTE {
    tag "cutadapt_test"
    //label 'process_single'
    
    input:
        path    input_fq
        val     sample
        path     outdir
        val     thread
        val     minimum_length
        val     nextseq_trim
        val     overlap
        //val     cutadapt_param

    output:
    path "${outdir}/${sample}_clean_2.fq" , emit: clean_r2

    script:
    
   // def other_args = ''
   // if (cutadapt_param != '') other_args += cutadapt_param

    """
    python /workspaces/scsnp/bin/cutadapt.py \
    --sample ${sample} \
    --outdir ${outdir}\
    --thread ${thread}\
    --input_fq ${input_fq} \
    --minimum_length ${minimum_length} \
    --nextseq_trim ${nextseq_trim} \
    --overlap ${overlap} 
    """
}


workflow  {

    main:
        CUTADAPTE(
           input_fq = Channel.fromPath( './out/barcode/test_2_2.fq' ),
            sample = "test_2",
            outdir = Channel.fromPath( "./out/cutadapt/"),
            thread = 4,
            minimum_length = params.minimum_length,
            nextseq_trim= params.nextseq_trim,
            overlap = params.overlap
        )
    
}
