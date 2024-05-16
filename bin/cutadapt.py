#!/usr/bin/env python
import subprocess
import json
from pymodules.minana import MinAna
import argparse

POLY_A = '-a polyA=A{18} '

LOG_METRICS_TITLE = (
    'Total reads processed',
    'Reads with adapters',
    'Reads that were too short',
    'Reads written (passing filters)',
    'Total basepairs processed',
    'Quality-trimmed',
    'Total written (filtered)',
)


def get_cutadapt_cmd(args, input_file, output_file):
    cmd = (
        'cutadapt '
        f'{POLY_A} '
        f'-j {args.thread} '
        f'-m {args.minimum_length} '
        f'--nextseq-trim={args.nextseq_trim} '
        f'--overlap {args.overlap} '
        f'{args.cutadapt_param} '
        f'--json {args.outdir}/{args.sample}_cutadapt.json '
        f'-o {output_file} '
        f'{input_file} '
        )
    return cmd

class Cutadapt(MinAna):
    """
    ## Features
    - Trim poly A tails and user-provided adapters in R2 reads with cutadapt.

    ## Output
    - `cutadapt.log` Cutadapt output log file.
    - `{sample}_clean_2.fq.gz` R2 reads file without adapters.
    """

    def __init__(self, args):
        super().__init__(args.sample, args.outdir)

        # out files
        self.out_fq2 = f'{self.outdir}/{self.sample}_clean_2.fq'
        self.json_log = f'{self.outdir}/{self.sample}_cutadapt.json'
        self.input_fq = args.input_fq
        self.args = args


    def add_cutadapt_metrics(self):
        with open(self.json_log) as f:
            log_dict = json.load(f)

        total_reads = log_dict['read_counts']['input']
        reads_with_adapters = log_dict['read_counts']['read1_with_adapter']
        reads_too_short = log_dict['read_counts']['filtered']['too_short']
        total_base_pairs = log_dict['basepair_counts']['input']
        quality_trimmed = log_dict['basepair_counts']['quality_trimmed']
        
        self.add_log_record(f'Reads with Adapters:{reads_with_adapters}({self.get_pct(reads_with_adapters, total_reads)})')
        self.add_log_record(f'Reads too Short:{reads_too_short}({self.get_pct(reads_too_short, total_reads)})')
        self.add_log_record(f'Base Pairs Quality-Trimmed:{quality_trimmed}({self.get_pct(quality_trimmed, total_base_pairs)})')

    def run(self):
        input_file = self.input_fq
        output_file = self.out_fq2
        cmd = get_cutadapt_cmd(self.args, input_file, output_file)
        subprocess.check_call(cmd, shell=True)
        self.add_cutadapt_metrics()
        self.write_log()

def main():
    parser = argparse.ArgumentParser(description='Cutadapt')

    parser.add_argument('--sample', help='', required=True)
    parser.add_argument('--outdir', help='', required=True)
    parser.add_argument('--thread', help='', default=4)
    parser.add_argument('--input_fq', help='Required. R2 reads from step Barcode.', required=True)

    parser.add_argument(
        '--minimum_length',
        help='Discard processed reads that are shorter than LENGTH.',
        default=20
    )

    parser.add_argument(
        '--nextseq_trim',
        help="""Quality trimming of reads using two-color chemistry (NextSeq). 
Some Illumina instruments use a two-color chemistry to encode the four bases. 
This includes the NextSeq and the NovaSeq. 
In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. 
However, dark cycles also occur when sequencing “falls off” the end of the fragment.
The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.""",
        default=20,
    )
    parser.add_argument(
        '--overlap',
        help="""Since Cutadapt allows partial matches between the read and the adapter sequence,
short matches can occur by chance, leading to erroneously trimmed bases. 
For example, roughly 0.25 of all reads end with a base that is identical to the first base of the adapter. 
To reduce the number of falsely trimmed bases, the alignment algorithm requires that 
at least {overlap} bases match between adapter and read. """,
        default=10
    )
    parser.add_argument('--cutadapt_param', help='Other cutadapt parameters. For example, --cutadapt_param "-a p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" ', default="")
    
    args = parser.parse_args()
    runner = Cutadapt(args)
    runner.run()

if __name__ == '__main__':
    main()