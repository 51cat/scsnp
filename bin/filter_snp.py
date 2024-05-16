#!/usr/bin/env python
import pysam
from pymodules.threshold import Threshold
from pymodules.minana import MinAna
import argparse

class SNPFilter(MinAna):
    """
    ## Features
    - Filter out `ref` and `alt` alleles that do not have enough reads to support.

    ## Output
    - `{sample}_test1_filtered.vcf` VCF file after filtering. Alleles read counts that do not have enough reads to support are set to zero. 
    Genotypes are changed accordingly.
    """
    def __init__(self, args):
        super().__init__(args.sample, args.outdir)
        
        self.vcf = args.vcf
        self.vaf = args.VAF
        self.ref_threshold_method = args.ref_threshold_method
        self.alt_threshold_method = args.alt_threshold_method
        self.ref_min_support_read = args.ref_min_support_read
        self.alt_min_support_read = args.alt_min_support_read

        self.add_log_record(
            (
                f"VAF value use:"
                f"{self.vaf}"
            )
        )    

        self.add_log_record(
            (
                f"reference allele threshold method:"
                f"{self.ref_threshold_method}"
            )
        )
        self.add_log_record(
            (
                f"alternate allele threshold method:"
                f"{self.alt_threshold_method}"
            )
        )
        # out
        self.out_vcf_file = f'{self.outdir}/{self.sample}_filtered.vcf'


    @staticmethod
    def get_count_array(record):
        """
        get allele count array

        Args
            record: pysam vcf record
            entry: str, 'ref' or 'alt'
        Return
            (ref_count_array, alt_count_array)
        """
        ref_count_array = []
        alt_count_array = []
        for sample in record.samples:
            ad = record.samples[sample]['AD']
            ref_count_array.append(ad[0])
            alt_count_array.append(ad[1])

        return ref_count_array, alt_count_array

    def get_threshold(self, count_array, threshold_method, otsu_plot_path=None):
        runner = Threshold(
            count_array,
            threshold_method=threshold_method,
            otsu_plot_path=otsu_plot_path
        )
        threshold = runner.run()
        return threshold

    @staticmethod
    def filter_array(count_array, threshold):
        return [x if x >= threshold else 0 for x in count_array]

    @staticmethod
    def get_genotype(ref_count, alt_count, VAF):
        
        if VAF >= 1:
            raise ValueError(f"VAF value must lower than 1!")

        if ref_count > 0 and alt_count > 0:
            af = alt_count / (alt_count + ref_count)
            if VAF <= af <= 1 - VAF:
                return (0,1)
            elif af < VAF:
                return (0,0)
            else:
                return (1,1)
        elif ref_count > 0 and alt_count == 0:
            return (0,0)
        elif ref_count == 0 and alt_count > 0:
            return (1,1)
        elif ref_count == 0 and alt_count == 0:
            return (None, None)

    def run(self):
        if self.ref_threshold_method == 'otsu':
            ref_otsu_plot_dir = f'{self.outdir}/ref_otsu_plots/'
            MinAna.check_mkdir(ref_otsu_plot_dir)
        
        if self.alt_threshold_method == 'otsu':
            alt_otsu_plot_dir = f'{self.outdir}/alt_otsu_plots/'
            MinAna.check_mkdir(alt_otsu_plot_dir)

        with pysam.VariantFile(self.vcf) as vcf_in:
            header = vcf_in.header
            header.add_meta('ref_threshold_method', value=self.ref_threshold_method)
            header.add_meta('alt_threshold_method', value=self.alt_threshold_method)
            header.add_meta('INFO', items=[('ID',"REF_T"), ('Number',1), ('Type','Integer'), ('Description','Reference allele count threshold')])
            header.add_meta('INFO', items=[('ID',"ALT_T"), ('Number',1), ('Type','Integer'), ('Description','Alternate allele count threshold')])
            with pysam.VariantFile(self.out_vcf_file, 'w', header=vcf_in.header) as vcf_out:
                for record in vcf_in.fetch():
                    name = f'{record.chrom}_{record.pos}'
                    ref_count_array, alt_count_array = self.get_count_array(record)
                    ref_otsu_plot_path = f'{ref_otsu_plot_dir}/{name}_ref.pdf' if self.ref_threshold_method == 'otsu' else None
                    ref_threshold = self.get_threshold(ref_count_array, self.ref_threshold_method, ref_otsu_plot_path)
                    alt_otsu_plot_path = f'{alt_otsu_plot_dir}/{name}_alt.pdf' if self.alt_threshold_method == 'otsu' else None
                    alt_threshold = self.get_threshold(alt_count_array, self.alt_threshold_method, alt_otsu_plot_path)
                    ref_threshold = max(ref_threshold, self.ref_min_support_read)
                    alt_threshold = max(alt_threshold, self.alt_min_support_read)
                    ref_filtered_count_array = self.filter_array(ref_count_array, ref_threshold)
                    alt_filtered_count_array = self.filter_array(alt_count_array, alt_threshold)

                    new_record = record.copy()
                    new_record.info.__setitem__('REF_T', ref_threshold)
                    new_record.info.__setitem__('ALT_T', alt_threshold)
                    for index, sample in enumerate(record.samples):
                        ref_count = ref_filtered_count_array[index]
                        alt_count = alt_filtered_count_array[index]
                        new_record.samples[sample]['AD'] = (ref_count, alt_count)
                        genotype = self.get_genotype(ref_count, alt_count, VAF = self.vaf)
                        new_record.samples[sample]['GT'] = genotype

                    vcf_out.write(new_record)
def main():

    parser = argparse.ArgumentParser(description='variant calling')
    parser.add_argument('--vcf', help='vcf file')
    parser.add_argument('--sample', help='sample name')
    parser.add_argument('--outdir', help='outdir')
    parser.add_argument('--ref_threshold_method', default='otsu', choices=['otsu', 'auto', 'none'], help='')
    parser.add_argument('--alt_threshold_method', default='otsu', choices=['otsu', 'auto', 'none'], help='')
    
    parser.add_argument(
        "--VAF",
        type=int,
        help='variant allele frequency (VAF) threshold', 
        default=0.2,
    )
    parser.add_argument(
        "--ref_min_support_read",
        type=int,
        help='minimum supporting read number for ref.', 
        default=2,
    )

    parser.add_argument(
        "--alt_min_support_read",
        type=int,
        help='minimum supporting read number for alt.', 
        default=2,
    )

    args = parser.parse_args()
    runner = SNPFilter(args)
    runner.run()

if __name__ == '__main__':
    main()