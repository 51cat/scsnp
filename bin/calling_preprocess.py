#!/usr/bin/env python
import pymodules.utils as utils
import argparse
from pymodules.minana import MinAna

class CallingProcesser(MinAna):
    """
    ## Features
    - Perform variant calling at single cell level.

    ## Output
    - `{sample}_raw.vcf` Variants are called with bcftools default settings.
    - `{sample}_norm.vcf` Indels are left-aligned and normalized. See https://samtools.github.io/bcftools/bcftools.html#norm for more details.
    """

    def __init__(self, args):
        super().__init__(args.sample, args.outdir)
        self.input_bam = args.input_bam
        self.fasta = args.fasta

        self.splitN_bam = f'{self.outdir}/{self.sample}_splitN.bam'

    def SplitNCigarReads(self):
        cmd = (
            f'gatk '
            f'SplitNCigarReads '
            f'--do-not-fix-overhangs '
            f'-R {self.fasta} '
            f'-I {self.input_bam} '
            f'-O {self.splitN_bam} '
        )
        utils.run_cmd(cmd)

    def run(self):
        self.SplitNCigarReads()


def main():
    parser = argparse.ArgumentParser(description='extract exon from bam')
    parser.add_argument('--input_bam', help='bamfile')
    parser.add_argument('--sample', help='sample name')
    parser.add_argument('--outdir', help='outdir')
    parser.add_argument('--fasta', help='Genome sequence')
    args = parser.parse_args()
    runner = CallingProcesser(args)
    runner.run()

if __name__ == '__main__':
    main()
