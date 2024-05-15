python ../barcode.py \
    --fq1 ./fastqs/snp1_1.fq.gz\
    --fq2  ./fastqs/snp1_2.fq.gz\
    --sample test_1\
    --outdir ./out/barcode/

python ../cutadapt.py \
    --input_fq ./out/barcode/test_1_2.fq\
    --sample test_1\
    --outdir ./out/cutadapt/

python ../star.py \
   --fq ./out/cutadapt/test_1_clean_2.fq\
   --sample test_1\
   --outdir ./out/star/\
   --genomeDir /SGRNJ06/randd/USER/zhouyiqi/genome/rna/celescope2.0.0/hs

python ../featureCount.py \
     --input_bam ./out/star/test_1_Aligned.out.bam\
     --sample test_1\
     --outdir ./out/featureCount/\
     --gtf /SGRNJ06/randd/USER/zhouyiqi/genome/rna/celescope2.0.0/hs/Homo_sapiens.GRCh38.99.filter.gtf \
     --gtf_type gene\
     --thread 4

python ../target_metrics.py \
    --input_bam ./out/featureCount/test_1_aligned_posSorted_addTag.bam\
    --sample test_1\
    --outdir ./out/target_metric_out/\
    --gene_list ./fastqs/gene_list.tsv\
    --match_dir ./fastqs/fake_match_dir \
    --add_RG

python ../calling_preprocess.py \
    --input_bam ./out/target_metric_out/test_1_filtered_sorted.bam\
    --sample test_1\
    --outdir ./out/calling_out\
    --fasta /SGRNJ06/randd/USER/zhouyiqi/genome/rna/celescope2.0.0/hs/Homo_sapiens.GRCh38.dna.primary_assembly.fa

python ../calling.py\
    --input_bam ./out/calling_out/test_1_splitN.bam\
    --sample test_1\
    --outdir ./out/calling_out\
    --fasta /SGRNJ06/randd/USER/zhouyiqi/genome/rna/celescope2.0.0/hs/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --thread 4

python ../filter_snp.py\
    --vcf ./out/calling_out/test_1_norm.vcf\
    --sample test_1\
    --outdir ./out/calling_out

python ../analysis.py \
    --input_vcf ./out/calling_out/test_1_filtered.vcf\
    --sample test_1\
    --outdir ./out/analysis_snp/\
    --gene_list ./fastqs/gene_list.tsv\
    --database GRCh38.mane.1.0.ensembl