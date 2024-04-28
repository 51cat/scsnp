
import bin.pymodules.utils as utils
import argparse
from pymodules.minana import MinAna

class VariantCaller(MinAna):
    """
    ## Features
    - Perform variant calling at single cell level.

    ## Output
    - `{sample}_raw.vcf` Variants are called with bcftools default settings.
    - `{sample}_norm.vcf` Indels are left-aligned and normalized. See https://samtools.github.io/bcftools/bcftools.html#norm for more details.
    """

    def __init__(self, args):
        # set
        super().__init__(args.sample, args.outdir)

        self.thread = args.thread
        self.input_bam = args.input_bam
        self.bed = args.bed_file
        self.fasta = args.fasta
        
        # out
        self.raw_bcf_file = f'{self.outdir}/{self.sample}_raw.bcf'
        self.raw_vcf_file = f'{self.outdir}/{self.sample}_raw.vcf'
        self.fixed_header_vcf = f'{self.outdir}/{self.sample}_fixed.vcf'
        self.norm_vcf_file = f'{self.outdir}/{self.sample}_norm.vcf'


    def call_variants(self):
        """
        max depth 100M
        """
        cmd = (
            f'bcftools mpileup '
            f'-f {self.fasta} '
            f'--threads {self.thread} '
            f'--annotate DP,AD -d 100000000 '
            f'-o {self.raw_bcf_file} '
            f'{self.input_bam} '
        )
        if self.bed:
            cmd += f' --regions-file {self.bed} '
        utils.run_cmd(cmd)

        cmd = (
            f'bcftools call '
            f'-mv -Ov '
            f'-o {self.raw_vcf_file} '
            f'{self.raw_bcf_file} '
        )
        utils.run_cmd(cmd)

    def bcftools_norm(self):
        cmd = (
            'bcftools norm '
            '-m- '
            f'-f {self.fasta} '
            f'{self.raw_vcf_file} '
            '| bcftools norm '
            '-d both '
            f'-o {self.norm_vcf_file} '
        )
        utils.run_cmd(cmd)

    def run(self):
        self.call_variants()
        self.bcftools_norm()


def main():
    parser = argparse.ArgumentParser(description='variant calling')
    parser.add_argument('--input_bam', help='bamfile')
    parser.add_argument('--sample', help='sample name')
    parser.add_argument('--outdir', help='outdir')
    parser.add_argument('--fasta', help='Genome sequence')
    parser.add_argument('--bed_file', help='Genome sequence')
    parser.add_argument('--thread', help='Genome sequence')
    args = parser.parse_args()
    runner = VariantCaller(args)
    runner.run()

if __name__ == '__main__':
    main()
