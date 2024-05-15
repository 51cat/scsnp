process BARCODE {
    tag "barcode_test"
    label 'process_single'
    
    input:
    path    fq1
    path    fq2
    path    outdir
    val     sample
    // chemistry
    val     chemistry
    val     pattern
    val     whitelist
    val     linker
    // qc
    val     lowQual
    val     lowNum
    // true store
    val     nopolyT
    val     noLinker
    val     filterNoPolyT
    val     allowNoLinker
    val     output_R1

    output:
    path "${outdir}/${sample}_2.fq" , emit: barcode_r2

    script:
    def ture_args = ""
    if (nopolyT == true) ture_args += " --nopolyT "
    if (noLinker == true) ture_args += " --noLinker "
    if (filterNoPolyT == true) ture_args += " --filterNoPolyT "
    if (allowNoLinker == true) ture_args += " --allowNoLinker "
    if (output_R1 == true) ture_args += " --output_R1 "
    def chem_args = ""
    if (pattern != '') chem_args += " --pattern ${pattern} "
    if (linker != '') chem_args += " --linker ${linker} "
    if (whitelist != '') chem_args += " --pattern ${whitelist} "
    if (chemistry == "auto") chem_args += " --chemistry auto "

    
    """
    python /workspaces/scsnp/bin/barcode.py \
    --fq1 ${fq1} \
    --fq2 ${fq2}\
    --sample ${sample}\
    --outdir ${outdir} \
    --lowQual ${lowQual} \
    --lowNum ${lowNum} \
    ${ture_args} \
    ${chem_args}
    """
}


workflow  {

    main:
        BARCODE(
            fq1 = in_fq1 = Channel.fromPath( '/workspaces/demo/data/celescope_test_data/snp/fastqs/snp1_1.fq.gz' ),
            fq2 = in_fq2 = Channel.fromPath( '/workspaces/demo/data/celescope_test_data/snp/fastqs/snp1_2.fq.gz' ),
            outdir = Channel.fromPath( './out/barcode/' ),
            sample = 'test_2',
            
            params.chemistry,
            params.pattern,
            params.whitelist,
            params.linker,

            params.lowQual,
            params.lowNum,
            
            params.nopolyT,
            params.noLinker,
            params.filterNoPolyT,
            params.allowNoLinker,
            params.output_R1
        )
    
}
