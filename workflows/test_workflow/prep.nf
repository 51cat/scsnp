/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN PREP WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FASTQC {
    tag "$meta.id"
    label 'process_single'

    //conda "bioconda::fastqc=0.12.1"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
    //    'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads), path(matchdirs)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Make list of old name and new name pairs to use for renaming in the bash while loop
    def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[ reads, "${prefix}.${reads}" ]] : reads.withIndex().collect { entry, index -> [ entry, "${prefix}.${entry}" ] }
    def rename_to = old_new_pairs*.join(' ').join(' ')
    def renamed_files = old_new_pairs.collect{ old_name, new_name -> new_name }.join(' ')
    
    """
    printf "%s %s\\n" $rename_to | while read old_name new_name; do
        [ -f "\${new_name}" ] || ln -s \$old_name \$new_name
    done

    fastqc \\
        $args \\
        --threads 4 \\
        $renamed_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
    
}

process FILTER_GTF {
    tag "$gtf"
    label 'process_single'

    //conda 'conda-forge::python==3.12'
    //container "biocontainers/python:3.12"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    path gtf
    val attributes

    output:
    path "*.filtered.gtf", emit: filtered_gtf
    path "gtf_filter.log", emit: log_file

    script:
    def args = task.ext.args ?: ''

    """
    filter_gtf.py ${gtf} \"${attributes}\"
    """
}

process STAR_GENOME {
    tag "$genome_name"
    label 'process_medium'

    //conda "bioconda::star==2.7.11b"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_0' :
    //    'biocontainers/star:2.7.11b--h43eeafb_0' }"

    input:
    path fasta
    path gtf
    val genome_name

    output:
    path "$genome_name"            , emit: index
    path "versions.yml"            , emit: versions

    script:
    def args        = task.ext.args ?: ''
    def memory      = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    def include_gtf = gtf ? "--sjdbGTFfile $gtf" : ''
    def fasta_sa = ( Math.log(fasta.size()) / Math.log(2) ) / 2 - 1
    def sa = Math.floor( Math.min(14, fasta_sa) )
    """
    mkdir ${genome_name}
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir ${genome_name}/ \\
        --genomeFastaFiles $fasta \\
        $include_gtf \\
        --runThreadN $task.cpus \\
        --genomeSAindexNbases ${sa} \\
        $memory \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}

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



process FEATURECOUNT {
    tag "FEATURECOUNT"
    label 'featureCount'

    input:
    path    input_bam
    path    gene_list
    val     sample
    val     outdir
    val     thread
    path    gtf
    path    gtf_type

    output:
    path "$outdir/${sample}_aligned_posSorted_addTag.bam" ,emit: featurecount_bam

    script:
    """
    featureCount.py \
    --input_bam ${input_bam}\
    --sample ${sample}\
    --outdir ${outdir} \
    --gtf ${gtf} \
    --gtf_type ${gtf_type}\
    --thread ${thread}
    """
}

workflow SNP_PREP {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
 
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    
    // filter gtf
    FILTER_GTF(
        channel.fromPath(params.gtf),
       params.keep_attributes
    )
    // star
    

    // cutadapt


    // STAR_GENOME

    // featurecount





}