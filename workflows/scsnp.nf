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
process BARCODE {
    tag "barcode_test"
    label 'process_singel'
    
    input:
    tuple val(meta), path(reads), path(matchdirs)
    val    outdir
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
    path "${outdir}/${meta.id}_2.fq" , emit: barcode_r2

    script:
    def (fq1, fq2) = reads.collate(2).transpose()
    def fq1_in = fq1[0]
    def fq2_in = fq2[0]
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
    --fq1 ${fq1_in} \
    --fq2 ${fq2_in}\
    --sample ${meta.id}\
    --outdir ${outdir} \
    --lowQual ${lowQual} \
    --lowNum ${lowNum} \
    ${ture_args} \
    ${chem_args}
    """
}

process CUTADAPTE {
    tag "cutadapt_test"
    label 'process_single'
    
    input:
    tuple val(meta), path(reads), path(matchdirs)
    path    input_fq
    val     outdir
    val     minimum_length
    val     nextseq_trim
    val     overlap
    //val     cutadapt_param

    output:
    path "${outdir}/${meta.id}_clean_2.fq" , emit: clean_r2

    script:
    
   // def other_args = ''
   // if (cutadapt_param != '') other_args += cutadapt_param

    """
    python /workspaces/scsnp/bin/cutadapt.py \
    --sample ${meta.id} \
    --outdir ${outdir}\
    --thread ${task.cpus}\
    --input_fq ${input_fq} \
    --minimum_length ${minimum_length} \
    --nextseq_trim ${nextseq_trim} \
    --overlap ${overlap} 
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

process STAR {
    tag "$genome_name"
    label 'process_high'

    input:
    tuple val(meta), path(reads), path(matchdirs)
    val  input_fq
    val  outdir
    path genomedir
    val  outFilterMatchNmin
    val  out_unmapped
    val  outFilterMultimapNmax
    val  consensus_fq

    output:
    path "${outdir}/${meta.id}_Aligned.out.bam", emit: star_bam

    script:
    def ture_args = ""
    if (consensus_fq == true) ture_args += " --consensus_fq "
    if (out_unmapped == true) ture_args += " --out_unmapped "

    """
    python /workspaces/scsnp/bin/star.py \
    --sample ${sample} \
    --outdir ${outdir}\
    --thread ${task.cpus}\
    --input_fq ${input_fq} \
    --genomeDir ${genomedir} \
    --outFilterMatchNmin ${outFilterMatchNmin} \
    --outFilterMultimapNmax ${outFilterMultimapNmax} \
    --starMem ${memory} \
    ${ture_args}
    """
}

process FEATURECOUNT {
    tag "$genome_name"
    label 'process_high'

    input:
    tuple val(meta), path(reads), path(matchdirs)
    path input_bam
    val  outdir
    path gtf
    val  gtf_type

    output:
    path "${outdir}/${sample}_nameSorted.bam", emit: featurecount_bam

    script:
    """
    python /workspaces/scsnp/bin/featureCount.py \
    --input_bam ${input_bam} \
    --sample ${meta.id} \
    --outdir ${outdir} \
    --thread ${task.cpus} \
    --gtf ${gtf} \
    --gtf_type ${gtf_type}
    """
}

process TARGET_METRICS {
    tag "target_meterics"
    label "process_single"

    // conda 'pysam=0.22.0 samtools'
    // container 

    input:
    tuple       val(meta), path(reads), path(matchdirs)
    path        bam 
    path        outdir
    path        gene_list

    output:
    path "$outdir/${sample}_filtered_sorted.bam" ,emit: target_bam

    script:

    """
    python /workspaces/scsnp/bin/target_metrics.py \\
    --input_bam $bam \\
    --sample $meta.id \\
    --outdir $outdir \\
    --gene_list $gene_list \\
    --match_dir $matchdirs  \\
    --add_RG
    """
}

process CALLING_PREPROCESS {
    tag "CALLING_PREPROCESS"
    label "process_medium"

    // conda
    // container 

    input:
    tuple       val(meta), path(reads), path(matchdirs)
    path        bam 
    path        outdir
    path        fasta

    output:
    path "$outdir/${sample}_splitN.bam" ,emit: exon_bam

    """
    python /workspaces/scsnp/bin/calling_preprocess.py \
    --input_bam $bam\
    --sample $meta.id\
    --outdir $outdir\
    --fasta $fasta
    """
}

process CALLING {
    tag "CALLING"
    label "process_high"

    // conda
    // container 

    input:
    tuple       val(meta), path(reads), path(matchdirs)
    path        bam 
    path        outdir
    path        fasta
    path        bed_file

    output:
    path "$outdir/${meta.id}_raw.bcf" ,emit: raw_bcf
    path "$outdir/${meta.id}_raw.vcf" ,emit: raw_vcf
    path "$outdir/${meta.id}_fixed.vcf" ,emit: fixed_vcf
    path "$outdir/${meta.id}_norm.vcf" ,emit: norm_vcf

    """
    python /workspaces/scsnp/bin/calling.py\
    --input_bam $bam\
    --sample $meta.id\
    --outdir $outdir\
    --fasta $fasta \
    --thread $task.cpus
    """
}

process FILTER_SNP {
    tag "FILTER_SNP"
    label "snp_filter"

    // conda
    // container 

    input:
    tuple       val(meta), path(reads), path(matchdirs)
    path        vcf
    path        outdir
    val         ref_threshold_method
    val         alt_threshold_method
    val         VAF
    val         ref_min_support_read
    val         alt_min_support_read

    output:
    path "$outdir/${meta.id}_filtered.vcf" ,emit: filtered_vcf

    """
    python /workspaces/scsnp/bin/filter_snp.py\
    --vcf $vcf\
    --sample $meta.id\
    --outdir $outdir \\
    --ref_threshold_method $ref_threshold_method \\
    --alt_threshold_method $alt_threshold_method \\
    --VAF $VAF \\
    --ref_min_support_read $ref_min_support_read \\
    --alt_min_support_read $alt_min_support_read
    """

}

process ANALYSIS {
    tag "ANALYSIS_SNP"
    label "process_medium"

    // conda
    // container 

    input:
    tuple       val(meta), path(reads), path(matchdirs)
    path        vcf 
    path        outdir
    path        gene_list
    val         database

    output:
    path        "${outdir}/${meta.id}_variant_table.csv", emit: variant_table

    script:
    """
    python /workspaces/scsnp/bin/analysis.py \
    --input_vcf ${vcf}\
    --sample ${meta.id}\
    --outdir ${outdir}\
    --gene_list ${gene_list}\
    --database ${database}
    """
}

workflow SCSNP {
    take:
        ch_samplesheet

    main:
        // fastqc 
        FASTQC(
            ch_samplesheet
        )

        // barcode
        BARCODE(
            ch_samplesheet,
            params.outdir,
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
            params.output_R1,
        )   
        r2 = BARCODE.out.barcode_r2
        
        // cutadapt
        CUTADAPTE(
            ch_samplesheet,
            r2,
            params.outdir,
            params.minimum_length,
            params.nextseq_trim,
            params.overlap
        )
        clean_r2 = CUTADAPTE.out.clean_r2

        // filter gtf
        FILTER_GTF(
            Channel.fromPath(params.gtf),
            params.keep_attributes
        )
        gtf_use = FILTER_GTF.out.filtered_gtf

        // index Gemome
        if (params.genomedir == null){
            STAR_GENOME(
                Channel.fromPath(params.fasta),
                Channel.fromPath(params.gtf),
                params.genome_name
            )
            genomeDir = STAR_GENOME.out.index
        }else {
            genomeDir = Channel.fromPath(genomedir)
        }

        // STAR
        STAR(
            ch_samplesheet,
            clean_r2,
            params.outdir,
            genomeDir,
            params.outFilterMatchNmin,
            params.out_unmapped,
            params.outFilterMultimapNmax,
            params.consensus_fq
        )
        star_bam = STAR.out.star_bam
        
        // featurecount
        FEATURECOUNT(
            ch_samplesheet,
            star_bam,
            params.outdir,
            params.gtf,
            params.gtf_type,
        )
        featurecount_bam = FEATURECOUNT.out.featurecount_bam

        // target meterics
        TARGET_METRICS(
            ch_samplesheet,
            featurecount_bam,
            params.outdir,
            params.gene_list
        )
        target_bam = TARGET_METRICS.out.target_bam

        //calling

        CALLING_PREPROCESS(
            ch_samplesheet,
            target_bam,
            params.outdir,
            params.fasta
        )
        exon_bam = CALLING_PREPROCESS.out.exon_bam
        
        CALLING(
            ch_samplesheet,
            exon_bam ,
            params.outdir,
            params.fasta,
            params.bed_file
        )
        norm_vcf = CALLING.out.norm_vcf
    
    // filter vcf

    FILTER_SNP(
            ch_samplesheet,
            norm_vcf,
            params.outdir,
            params.ref_threshold_method,
            params.alt_threshold_method,
            params.VAF,
            params.ref_min_support_read,
            params.alt_min_support_read
    )
    filtered_vcf = FILTER_SNP.out.filtered_vcf

    // analysis vcf
    ANALYSIS(
        ch_samplesheet,
            filtered_vcf,
            params.outdir,
            params.gene_list,
            params.database,
    )

    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
