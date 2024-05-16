#!/usr/bin/env python
import os
import pathlib
from collections import defaultdict
from itertools import groupby
import subprocess
import shutil
import pysam
import argparse
from pymodules.minana import MinAna
from pymodules.gtf_parser import GtfParser
from pymodules.utils import samtools

GTF_TYPES = ['exon','gene']
TAG_BAM_SUFFIX = 'aligned_posSorted_addTag.bam'

def genDict(dim=3, valType=int):
    if dim == 1:
        return defaultdict(valType)
    else:
        return defaultdict(lambda: genDict(dim - 1, valType=valType))

class FeatureCounts(MinAna):

    """
    ## Features
    - Assigning uniquely mapped reads to genomic features with FeatureCounts.
    ## Output
    - `{sample}` Numbers of reads assigned to features (or meta-features).
    - `{sample}_summary` Stat info for the overall summrization results, including number of 
    successfully assigned reads and number of reads that failed to be assigned due to 
    various reasons (these reasons are included in the stat info).
    - `{sample}_aligned_sortedByCoord_addTag.bam` featureCounts output BAM, 
    sorted by coordinates
    """

    def __init__(self,args):
        super().__init__(args.sample, args.outdir)

        self.input_bam = args.input_bam
        self.thread = args.thread
        self.gtf = args.gtf
        self.gtf_type = args.gtf_type
        gp = GtfParser(self.gtf)
        self.id_name = gp.get_id_name()
        self.intron_dict = {}

        # stats
        self.exon = self.intron = self.intergenic = self.ambiguity = 0


        # temp file
        self.tmp_dir = f'{self.outdir}/tmp/'
        self.add_tag_bam = f'{self.outdir}/{self.sample}_addTag.bam'
        input_basename = os.path.basename(self.input_bam)
        self.exon_bam = f'{self.tmp_dir}/exon/{input_basename}.featureCounts.bam'
        self.intron_bam = f'{self.tmp_dir}/intron/{input_basename}.featureCounts.bam'

        # out
        self.count_detail_file = f'{self.outdir}/{self.sample}_count_detail.txt'
        self.nameSorted_bam = f'{self.outdir}/{self.sample}_nameSorted.bam'
        self.out_bam = f'{self.outdir}/{self.sample}_{TAG_BAM_SUFFIX}'

    def add_tag(self, seg, id_name):
        """
        Add intron reads and tag

        Args:
            seg: pysam bam segment
            id_name: {gene_id: gene_name}

        Returns:
            seg with tag added

        Tags:
            CB: cell barcode
            UB: error-corrected UMI
            UR: original UMI
            GN: gene name
            GX: gene_id

        """
        attr = seg.query_name.split(':')
        barcode = attr[0]
        ur = ub = attr[1]

        # assign to some gene
        xs = seg.get_tag('XS')
        if xs == 'Assigned':
            gene_id = seg.get_tag('XT')
            gene_name = id_name[gene_id]
            seg.set_tag(tag='GN', value=gene_name, value_type='Z')
            seg.set_tag(tag='GX', value=gene_id, value_type='Z')
            seg.set_tag(tag='RE', value='E', value_type='Z')
            self.exon += 1
        else:
            if self.intron_dict and seg.query_name in self.intron_dict:
                gene_id = self.intron_dict[seg.query_name]
                gene_name = id_name[gene_id]
                seg.set_tag(tag='GN', value=gene_name, value_type='Z')
                seg.set_tag(tag='GX', value=gene_id, value_type='Z')
                seg.set_tag(tag='RE', value='N', value_type='Z')
                seg.set_tag(tag='XT', value=gene_id, value_type='Z')
                self.intron += 1
            elif xs == 'Unassigned_NoFeatures':
                seg.set_tag(tag='RE', value='I', value_type='Z')
                self.intergenic += 1
            elif xs == 'Unassigned_Ambiguity':
                seg.set_tag(tag='RE', value='A', value_type='Z')
                self.ambiguity += 1
        
        seg.set_tag(tag='CB', value=barcode, value_type='Z')
        seg.set_tag(tag='UB', value=ub, value_type='Z')
        seg.set_tag(tag='UR', value=ur, value_type='Z')

        return seg

    def get_intron_dict(self):
        with pysam.AlignmentFile(self.intron_bam, "rb") as in_bam:
            for seg in in_bam:
                if seg.has_tag('XT'):
                    self.intron_dict[seg.query_name] = seg.get_tag('XT')

    def run_featureCounts(self, outdir, gtf_type):
        '''
        allow multimapping with -M; but each multi-mapped reads only have one alignment because of --outSAMmultNmax 1
        '''
        cmd = (
            'featureCounts '
            f'-s 1 '
            f'--largestOverlap '
            f'-M '
            f'-a {self.gtf} '
            f'-o {outdir}/{self.sample} '  
            '-R BAM '
            f'-T {self.thread} '
            f'-t {gtf_type} '
            f'{self.input_bam} '
            '2>&1 '
        )
        #if self.args.featureCounts_param:
        #    cmd += (" " + self.args.featureCounts_param)
        
        subprocess.check_call(cmd, shell=True)

    def run_exon_intron(self):
        tmp_dir = f'{self.outdir}/tmp/'
        for gtf_type in ['exon', 'intron']:
            outdir = f'{tmp_dir}/{gtf_type}'
            pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)            
            self.run_featureCounts(outdir, gtf_type)

    def remove_temp_file(self):
        shutil.rmtree(self.tmp_dir)
        os.remove(self.add_tag_bam)

    def run(self):
        self.run_exon_intron()
        self.get_intron_dict()

        samtools(
            self.exon_bam,
            "sort",
            self.nameSorted_bam,
            'name'
            )
        
        self.get_count_detail_add_tag()
        self.add_metrics()
        
        samtools(
            inputbam=self.add_tag_bam,
            do="sort",
            outputbam=self.out_bam)
        
        self.remove_temp_file()

        self.write_log()

    def add_metrics(self):
        total = self.exon + self.intron + self.intergenic + self.ambiguity
        self.add_log_record(
            (
                f"Feature Type: "
                f"{self.gtf_type.capitalize()} "
            )
        )
        self.add_log_record(
            (
                f"Reads Assigned To Exonic Regions: "
                f"{self.exon}({round(100*self.exon/total, 2)} %)"
            )
        )

        self.add_log_record(
            (
                f"Reads Assigned To Intronic Regions: "
                f"{self.intron}({round(100*self.intron/total, 2)}%)"
            )
        )

        self.add_log_record(
            (
                f"Reads Assigned To Intergenic Regions: "
                f"{self.intergenic}({round(100*self.intergenic/total, 2)}%)"
            )
        )

        self.add_log_record(
            (
                f"Reads Unassigned Ambiguity: "
                f"{self.ambiguity}({round(100*self.ambiguity/total, 2)}%)"
            )
        )

    def get_count_detail_add_tag(self):
        """
        bam to detail table
        must be used on name_sorted bam
        Output file:
            - count_detail_file
            - bam with tag(remain name sorted)
        """
        save = pysam.set_verbosity(0)
        inputFile = pysam.AlignmentFile(self.nameSorted_bam, "rb")
        outputFile = pysam.AlignmentFile(self.add_tag_bam, 'wb', header=inputFile.header)
        pysam.set_verbosity(save)

        with open(self.count_detail_file, 'wt') as fh1:
            fh1.write('\t'.join(['Barcode', 'geneID', 'UMI', 'read', 'unique', 'PCR_duplicate']) + '\n')

            def keyfunc(x):
                return x.query_name.split(':', 1)[0]
            for _, g in groupby(inputFile, keyfunc):
                gene_umi_dict = defaultdict(lambda: defaultdict(int))
                gene_umi_pos = genDict(dim=3, valType=int)
                for seg in g:
                    seg = self.add_tag(seg, self.id_name)
                    outputFile.write(seg)
                    barcode, umi = seg.get_tag('CB'), seg.get_tag('UB')
                    if not seg.has_tag('GX'):
                        continue
                    gene_id = seg.get_tag('GX')
                    gene_umi_dict[gene_id][umi] += 1
                    gene_umi_pos[gene_id][umi][seg.reference_start] += 1                

                # output
                for gene_id in gene_umi_dict:
                    n_umi = len(gene_umi_dict[gene_id])
                    n_read = 0
                    unique = dup = 0
                    for umi in gene_umi_dict[gene_id]:
                        read_count = gene_umi_dict[gene_id][umi] 
                        n_read += read_count
                        if read_count == 1:
                            # unique
                            unique += 1
                        else:
                            # only add postion duplicate read number
                            for pos in gene_umi_pos[gene_id][umi]:
                                if gene_umi_pos[gene_id][umi][pos] > 1:
                                    dup += gene_umi_pos[gene_id][umi][pos]
                    fh1.write(f'{barcode}\t{gene_id}\t{n_umi}\t{n_read}\t{unique}\t{dup}\n')

        inputFile.close()
        outputFile.close()

def main():
    parser = argparse.ArgumentParser(description='FeatureCount')
    parser.add_argument('--input_bam', help='bamfile')
    parser.add_argument('--sample', help='sample name')
    parser.add_argument('--outdir', help='outdir')
    parser.add_argument('--thread', help='', type=int, default=4)
    parser.add_argument('--gtf', help='gtf file')

    parser.add_argument('--gtf_type',
                        help='Specify feature type in GTF annotation', 
                        default='gene',choices=['exon', 'gene'])


    args = parser.parse_args()
    runner = FeatureCounts(args)
    runner.run()

if __name__ == '__main__':
    main()