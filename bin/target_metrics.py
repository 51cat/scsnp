#!/usr/bin/env python
import pysam
import pymodules.utils as utils
from pymodules.minana import MinAna
import argparse

class TargetMetrics(MinAna):
    """
    ## Features
    - Filter bam file
        - Filter reads that are not cell-associated.
        - Filter reads that are not mapped to target genes. 

    - Collect enrichment metrics.

    ## Output
    - `filtered.bam` BAM file after filtering. Reads that are not cell-associated or not mapped to target genes are filtered.
    """

    def __init__(self, args):
        super().__init__(args.sample, args.outdir)

        self.input_bam = args.input_bam
        self.add_RG = args.add_RG

        self.match_barcode_list, self.n_cell = utils.get_barcode_from_matrix_dir(args.match_dir)
        self.match_barcode = set(self.match_barcode_list)
        self.gene_list = utils.read_one_col(args.gene_list)
        self.n_gene = len(self.gene_list)

        self.cell = set()
        self.cell_valid = set()
        self.enrich_reads_incell = 0
        self.enrich_reads = 0

        # out file
        self.out_bam_file = f'{self.outdir}/{self.sample}_filtered.bam'
        self.out_bam_file_sorted = f'{self.outdir}/{self.sample}_filtered_sorted.bam'

    def read_bam_write_filtered(self):
        sam_temp = f'{self.out_bam_file}.temp'
        with pysam.AlignmentFile(self.input_bam, "rb") as reader:
            header = reader.header.to_dict()
            # add RG to header
            if self.add_RG:
                header['RG'] = []
                for barcode in self.match_barcode_list:
                    header['RG'].append({
                        'ID': barcode,
                        'SM': barcode,
                    })
            with pysam.AlignmentFile(sam_temp, "w", header=header) as writer:
                for record in reader:
                    try:
                        gene_name = record.get_tag('GN')
                    except KeyError:
                        continue
                    # compatible with tag bam
                    try:
                        barcode = record.get_tag('CB')
                        self.cell.add(barcode)
                        self.enrich_reads += 1

                        UMI = record.get_tag('UB')
                    except KeyError:
                        continue
                    if barcode in self.match_barcode and gene_name in self.gene_list:
                        self.cell_valid.add(barcode)
                        self.enrich_reads_incell += 1
                        if self.add_RG:
                            record.set_tag(tag='RG', value=record.get_tag('CB'), value_type='Z')
                        writer.write(record)

            cmd = f'samtools view -b {sam_temp} -o {self.out_bam_file}; rm {sam_temp}'
            utils.run_cmd(cmd)

    def add_log(self):
        ## 有问题
        total_cell = len(self.cell)
        total_valid_cell = len(self.cell_valid)
        self.add_log_record(f"Number of Target Genes: {self.n_gene}")
        self.add_log_record(f"Number of Cells: {len(self.cell)}")
        self.add_log_record(f"Number of Valid Cells: {len(self.cell_valid)}({self.get_pct(total_valid_cell, total_cell)})")
        self.add_log_record(f"Enriched Reads: {self.enrich_reads }")
        self.add_log_record(f"Enriched Reads in Cells: {self.enrich_reads_incell }")

    def run(self):
        self.read_bam_write_filtered()
        # sort 
        utils.samtools(
            inputbam=self.out_bam_file,
            do = "sort",
            outputbam=self.out_bam_file_sorted
        )

        # index
        utils.samtools(
            inputbam=self.out_bam_file_sorted,
            do = "index"
        )

        self.add_rubbish(self.out_bam_file)
        self.clean()
        self.add_log()
        self.write_log()


def main():
    parser = argparse.ArgumentParser(description='variant calling')
    parser.add_argument('--input_bam', help='bamfile', required=True)
    parser.add_argument('--sample', help='sample name', required=True)
    parser.add_argument('--outdir', help='outdir', required=True)
    parser.add_argument("--gene_list", help='', required=True)
    parser.add_argument('--match_dir', help='', required=True)
    parser.add_argument(
            '--add_RG', help='Add tag read group: RG. RG is the same as CB(cell barcode)', action='store_true')
    args = parser.parse_args()
    runner = TargetMetrics(args)
    runner.run()

if __name__ == '__main__':
    main()
