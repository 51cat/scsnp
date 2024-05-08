import subprocess
import pandas as pd
from pymodules.minana import MinAna
import  argparse

STAR_BAM_SUFFIX = 'Aligned.out.bam'

def get_star_cmd(args, input_file, output_prefix):
    """
    output sam format to improve speed
    """
    cmd = (
        f'STAR '
        f'--runThreadN {args.thread} '
        f'--genomeDir {args.genomeDir} '
        f'--outSAMmultNmax 1 '
        f'--outFilterMultimapNmax {args.outFilterMultimapNmax} '
        f'--outSAMtype BAM Unsorted '
        f'--outFilterMatchNmin {args.outFilterMatchNmin} '
        f'{args.STAR_param} '
        f'--readFilesIn {input_file} '
        f'--outFileNamePrefix {output_prefix}_ '
    )
    return cmd


def get_star_log(logPath):
    df = pd.read_csv(logPath, sep='\t',header=None, names=['name','value'])
    log_dict = {}
    for t in df.itertuples():
        name = t.name.strip('|').strip()
        log_dict[name] = t.value

    return log_dict

class Star_mixin(MinAna):
    """
    Mixin class for STAR
    """

    def __init__(self, args):
        super().__init__(args.sample, args.outdir)

        # parse
        #self.stat_prefix = 'Reads'
        #if getattr(args, "consensus_fq", False):
        #    self.stat_prefix = 'UMIs'

        # out
        #if add_prefix:
        #    self.out_prefix += f'_{add_prefix}'
        self.args = args
        
        self.logPath = f'{self.outdir}/{self.sample}_Log.final.out'
        self.STAR_bam = f'{self.outdir}/{self.sample}_{STAR_BAM_SUFFIX}'

    def STAR(self):
        input_file = self.args.fq
        output_prefix = f"{self.outdir}/{self.sample}"
        cmd = get_star_cmd(self.args, input_file, output_prefix)
        
        if self.args.out_unmapped:
            cmd += '--outReadsUnmapped Fastx'
        if self.args.fq[-3:] == ".gz":
            cmd += '--readFilesCommand zcat'
        if self.args.STAR_param:
            cmd += (" " + self.args.STAR_param)
        self.add_log_record(f"cmd: {cmd}")
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.STAR()
        self.add_star_metrics()
        self.write_log()

    def add_star_metrics(self):
        """
        step metrics
        """
        log_dict = get_star_log(self.logPath)

        unique_reads = int(log_dict['Uniquely mapped reads number'])
        multiple_loci_reads = int(log_dict['Number of reads mapped to multiple loci'])
        too_many_loci_reads = int(log_dict['Number of reads mapped to too many loci'])
        total_reads = int(log_dict['Number of input reads'])
        multi_reads = multiple_loci_reads + too_many_loci_reads

        #self.add_metric(
        #    name='Genome',
        #    value=self.genome_name,
        #)
        
        self.add_log_record(f'Uniquely Mapped Reads: {unique_reads}({self.get_pct(unique_reads, total_reads)})')
        self.add_log_record(f'Multi-Mapped Reads: {multi_reads}({self.get_pct(multi_reads, total_reads)})')


def main():
    parser = argparse.ArgumentParser(description='Cutadapt')
    parser.add_argument('--sample', help="", required=True)
    parser.add_argument('--outdir', help="", required=True)
    parser.add_argument('--thread', help='', default=4)

    parser.add_argument(
        '--genomeDir',
        help='',
    )
    parser.add_argument(
        '--outFilterMatchNmin',
        help="""Alignment will be output only if the number of matched bases 
is higher than or equal to this value.""",
        default=50,
    )
    parser.add_argument(
        '--out_unmapped',
        help='Output unmapped reads.',
        action='store_true'
    )
    parser.add_argument('--STAR_param', help='', default="")
    parser.add_argument(
        '--outFilterMultimapNmax',
        help=(
            'How many places are allowed to match a read at most. Please note that even if this value is greater than 1, '
            'at most 1 alignment(with the highest score) will be output in the bam file.'
        ),
        default=1,
    )
    parser.add_argument(
        '--starMem',
        help='Default `30`. Maximum memory that STAR can use.',
        default=30
    )
    parser.add_argument('--fq', help="Required. R2 fastq file.", required=True)
    parser.add_argument("--consensus_fq", action='store_true',
                            help="A indicator that the input fastq has been consensused.")
    args = parser.parse_args()
    runner = Star_mixin(args)
    runner.run()

if __name__ == '__main__':
    main()