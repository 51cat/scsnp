#!/usr/bin/env python
import glob
import re
import sys
from collections import Counter, defaultdict
from itertools import combinations, product
from pymodules.minana import MinAna
import pysam
from xopen import xopen
import argparse
import pymodules.utils as utils
import os


BIN_DIR = os.path.dirname(os.path.abspath(__file__))
ASSETS_DIR = os.path.join(os.path.dirname(BIN_DIR), "assets")

MIN_T = 10

PATTERN_DICT = {
    'auto': None,
    'scopeV1': 'C12U8T18',
    'scopeV2.0.0': 'C8L16C8L16C8U8T18',
    'scopeV2.0.1': 'C8L16C8L16C8L1U8T18',
    'scopeV2.1.0': 'C8L16C8L16C8U12T18',
    'scopeV2.1.1': 'C8L16C8L16C8L1U12T18',
    'scopeV2.2.1': 'C8L16C8L16C8L1U12T18',
    'scopeV3.0.1': 'C9L16C9L16C9L1U12T18',
    'customized': None,
}

class Chemistry():
    """
    Auto detect chemistry from R1-read
    """

    def __init__(self, fq1, assay=None):
        '''
        'scopeV2.0.1': 'C8L16C8L16C8L1U8T18'
        'scopeV2.1.1': 'C8L16C8L16C8L1U12T18'
        'scopeV2.2.1': 'C8L16C8L16C8L1U12T18' with 4 types of linkers
        'scopeV3.0.1': 'C9L16C9L16C9L1U12T18' with 4 types of linkers
        '''
        self.fq1 = fq1
        self.assay = assay
        self.fq1_list = fq1.split(',')
        self.n_read = 10000

        self.pattern_dict_v2, * \
            _, self.linker_1_v2_set_list, self.linker_1_v2_mismatch_list = Barcode.parse_chemistry('scopeV2.1.1')
        self.pattern_dict_v2, * \
            _, self.linker_4_v2_set_list, self.linker_4_v2_mismatch_list = Barcode.parse_chemistry('scopeV2.2.1')
        self.pattern_dict_v3, *_, self.linker_v3_set_list, self.linker_v3_mismatch_list = Barcode.parse_chemistry('scopeV3.0.1')

    def check_chemistry(self):
        """check chemistry in the fq1_list"""
        chemistry_list = []
        for fastq1 in self.fq1_list:
        #    self.check_chemistry.logger.info(fastq1)
            chemistry = self.get_chemistry(fastq1)
            chemistry_list.append(chemistry)
        #if len(set(chemistry_list)) != 1:
        #    Chemistry.check_chemistry.logger.warning('multiple chemistry found!' + str(chemistry_list))
        return chemistry_list


    def seq_chemistry(self, seq):

        linker_v2 = Barcode.get_seq_str(seq, self.pattern_dict_v2["L"])
        bool_valid, _, _ = Barcode.check_seq_mismatch(
            [linker_v2], self.linker_1_v2_set_list, self.linker_1_v2_mismatch_list)
        if bool_valid:
            if seq[65:69] == "TTTT":
                return "scopeV2.0.1"
            else:
                return "scopeV2.1.1"

        linker_v3 = Barcode.get_seq_str(seq, self.pattern_dict_v3["L"])
        bool_valid, _, _ = Barcode.check_seq_mismatch(
            [linker_v3], self.linker_v3_set_list, self.linker_v3_mismatch_list)
        if bool_valid:
            return "scopeV3.0.1"

        linker_v2 = Barcode.get_seq_str(seq, self.pattern_dict_v2["L"])
        bool_valid, _, _ = Barcode.check_seq_mismatch(
            [linker_v2], self.linker_4_v2_set_list, self.linker_4_v2_mismatch_list)
        if bool_valid:
            return "scopeV2.2.1"

        return

    def get_chemistry(self, fq1):
        results = defaultdict(int)

        with pysam.FastxFile(fq1) as fh:
            for _ in range(self.n_read):
                entry = fh.__next__()
                seq = entry.sequence
                chemistry = self.seq_chemistry(seq)
                if chemistry:
                    results[chemistry] += 1
        # if it is 0, then no other linker types
        if results["scopeV2.2.1"] != 0:
            results["scopeV2.2.1"] += results["scopeV2.1.1"]
        sorted_counts = sorted(results.items(), key=lambda x: x[1], reverse=True)
        #self.get_chemistry.logger.info(sorted_counts)

        chemistry, read_counts = sorted_counts[0][0], sorted_counts[0][1]
        percent = float(read_counts) / self.n_read
        if percent < 0.5:
            print("Valid chemistry read counts percent < 0.5")
        if percent < 0.1:
            #self.get_chemistry.logger.error("Valid chemistry read counts percent < 0.1")
            raise Exception(
                'Auto chemistry detection failed! ' #+ HELP_DICT['chemistry']
            )
        #Chemistry.get_chemistry.logger.info(f'chemistry: {chemistry}')

        return chemistry


class Barcode(MinAna):
    """
    ## Features

    - Demultiplex barcodes.
    - Filter invalid R1 reads, which includes:
        - Reads without linker: the mismatch between linkers and all linkers in the whitelist is greater than 2.  
        - Reads without correct barcode: the mismatch between barcodes and all barcodes in the whitelist is greater than 1.  
        - Reads without polyT: the number of T bases in the defined polyT region is less than 10.
        - Low quality reads: low sequencing quality in barcode and UMI regions.

    ## Output

    - `01.barcode/{sample}_2.fq(.gz)` Demultiplexed R2 reads. Barcode and UMI are contained in the read name. The format of 
    the read name is `{barcode}_{UMI}_{read ID}`.
    """

    def __init__(self, args):
        super().__init__(args.sample, args.outdir)

        self.args = args
        self.assay = "snp"
        self.fq1_list = args.fq1.split(",")
        self.fq2_list = args.fq2.split(",")
        self.fq_number = len(self.fq1_list)
        if self.fq_number != len(self.fq2_list):
            raise Exception('fastq1 and fastq2 do not have same file number!')
        if args.chemistry == 'auto':
            ch = Chemistry(args.fq1, self.assay)
            self.chemistry_list = ch.check_chemistry()
        else:
            self.chemistry_list = [args.chemistry] * self.fq_number
        self.add_log_record(f"{','.join(self.chemistry_list)}")
        self.barcode_corrected_num = 0
        self.linker_corrected_num = 0
        self.total_num = 0
        self.clean_num = 0
        self.no_polyT_num = 0
        self.lowQual_num = 0
        self.no_linker_num = 0
        self.no_barcode_num = 0
        self.barcode_qual_Counter = Counter()
        self.umi_qual_Counter = Counter()
        self.pattern = args.pattern
        self.linker = args.linker
        self.whitelist = args.whitelist
        self.lowNum = args.lowNum
        self.lowQual = args.lowQual
        self.filterNoPolyT = args.filterNoPolyT
        self.allowNoLinker = args.allowNoLinker
        self.nopolyT = args.nopolyT  # true == output nopolyT reads
        self.noLinker = args.noLinker
        self.output_R1 = args.output_R1
        #self.wells = args.wells

        # out file
        self.out_fq2 = f'{self.outdir}/{self.sample}_2.fq'
        self.out_fq1 = f'{self.outdir}/{self.sample}_1.fq'
        if self.nopolyT:
            self.nopolyT_1 = f'{self.outdir}/{self.sample}_noPolyT_1.fq'
            self.nopolyT_2 = f'{self.outdir}/{self.sample}_noPolyT_2.fq'
        if self.noLinker:
            self.noLinker_1 = f'{self.outdir}/{self.sample}_noLinker_1.fq'
            self.noLinker_2 = f'{self.outdir}/{self.sample}_noLinker_2.fq'

        self.open_files()

    @staticmethod
    def get_seq_str_no_exception(seq, sub_pattern_dict):
        """get subseq with intervals in arr and concatenate"""
        return ''.join([seq[item[0]: item[1]] for item in sub_pattern_dict])

    @staticmethod
    def get_seq_str(seq, sub_pattern_dict):
        """
        Get subseq with intervals in arr and concatenate

        Args:
            seq: str
            sub_pattern_dict: [[0, 8], [24, 32], [48, 56]]

        Returns:
            str
            if sequence length is not enough, return ""

        >>> sub_pattern_dict = [[0, 8]]
        >>> seq = "A" * 7
        >>> Barcode.get_seq_str(seq, sub_pattern_dict)
        ""
        >>> seq = "A" * 8
        >>> Barcode.get_seq_str(seq, sub_pattern_dict)
        'AAAAAAAA'
        """
        seq_len = len(seq)
        ans = []
        for item in sub_pattern_dict:
            start, end = item[0], item[1]
            if end > seq_len:
                raise IndexError(f"sequence length is not enough in R1 read: {seq}")
            else:
                ans.append(seq[start:end])
        return ''.join(ans)

    @staticmethod
    def get_seq_list(seq, pattern_dict, abbr):
        """
        >>> pattern_dict = Barcode.parse_pattern("C2L3C2")
        >>> seq = "AAGGGTT"
        >>> Barcode.get_seq_list(seq, pattern_dict, "C")
        ['AA', 'TT']
        """
        
        return [seq[item[0]: item[1]] for item in pattern_dict[abbr]]

    @staticmethod
    def parse_pattern(pattern):
        """
        >>> pattern_dict = Barcode.parse_pattern("C8L16C8L16C8L1U12T18")
        >>> pattern_dict['C']
        [[0, 8], [24, 32], [48, 56]]
        >>> pattern_dict['L']
        [[8, 24], [32, 48], [56, 57]]
        """
        pattern_dict = defaultdict(list)
        p = re.compile(r'([CLUNT])(\d+)')
        tmp = p.findall(pattern)
        if not tmp:
            sys.exit()(f'Invalid pattern: {pattern}')
        start = 0
        for item in tmp:
            end = start + int(item[1])
            pattern_dict[item[0]].append([start, end])
            start = end
        return pattern_dict

    @staticmethod
    def get_abbr_len(pattern_dict, abbr):
        """
        >>> pattern_dict = Barcode.parse_pattern("C8L16C8L16C8L1U12T18")
        >>> Barcode.get_abbr_len(pattern_dict, 'C')
        24
        >>> Barcode.get_abbr_len(pattern_dict, 'L')
        33
        """
        length = 0
        for item in pattern_dict[abbr]:
            length += item[1] - item[0]

        return length
        

    @staticmethod
    def get_scope_bc(chemistry, assets_dir = ASSETS_DIR):
        """Return (linker file path, whitelist file path)"""
        try:
            linker_f = glob.glob(f'{assets_dir}/whitelist/{chemistry}/linker*')[0]
            whitelist_f = f'{assets_dir}/whitelist/{chemistry}/bclist'
        except IndexError:
            return None, None
        return linker_f, whitelist_f

    @staticmethod
    def ord2chr(q, offset=33):
        return chr(int(q) + offset)

    @staticmethod
    def qual_int(char, offset=33):
        return ord(char) - offset

    @staticmethod
    def low_qual(quals, minQ, num):
        # print(ord('/')-33)           14
        return True if len([q for q in quals if Barcode.qual_int(q) < minQ]) > num else False

    @staticmethod
    def findall_mismatch(seq, n_mismatch=1, bases='ACGTN'):
        """
        choose locations where there's going to be a mismatch using combinations
        and then construct all satisfying lists using product

        Return:
        all mismatch <= n_mismatch set. 

        >>> answer = set(["TCG", "AAG", "ACC", "ATG", "ACT", "ACN", "GCG", "ANG", "ACA", "ACG", "CCG", "AGG", "NCG"])
        >>> seq_set = Barcode.findall_mismatch("ACG")
        >>> seq_set == answer
        True
        """
        seq_set = set()
        seq_len = len(seq)
        if n_mismatch > seq_len:
            n_mismatch = seq_len
        for locs in combinations(range(seq_len), n_mismatch):
            seq_locs = [[base] for base in seq]
            for loc in locs:
                seq_locs[loc] = list(bases)
            for poss in product(*seq_locs):
                seq_set.add(''.join(poss))
        return seq_set

    @staticmethod
    def get_mismatch_dict(seq_list, n_mismatch=1):
        """
        Return:
        mismatch dict. Key: mismatch seq, value: seq in seq_list

        >>> seq_list = ["AACGTGAT", "AAACATCG"]
        >>> mismatch_dict = Barcode.get_mismatch_dict(seq_list)
        >>> mismatch_dict["AACGTGAA"] == "AACGTGAT"
        True
        """
        mismatch_dict = {}

        for seq in seq_list:
            seq = seq.strip()
            if seq == '':
                continue
            for mismatch_seq in Barcode.findall_mismatch(seq, n_mismatch):
                mismatch_dict[mismatch_seq] = seq

        return mismatch_dict

    @staticmethod
    def check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list):
        '''
        Return bool_valid, bool_corrected, corrected_seq

        >>> seq_list = ['ATA', 'AAT', 'ATA']
        >>> correct_set_list = [{'AAA'},{'AAA'},{'AAA'}]
        >>> mismatch_dict_list = [Barcode.get_mismatch_dict(['AAA'])] * 3

        >>> Barcode.check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list)
        (True, True, 'AAA_AAA_AAA')

        >>> seq_list = ['AAA', 'AAA', 'AAA']
        >>> Barcode.check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list)
        (True, False, 'AAA_AAA_AAA')
        '''
        bool_valid = True
        bool_corrected = False
        corrected_seq_list = []
        for index, seq in enumerate(seq_list):
            if seq not in correct_set_list[index]:
                if seq not in mismatch_dict_list[index]:
                    bool_valid = False
                    return bool_valid, bool_corrected, ""
                else:
                    bool_corrected = True
                    corrected_seq_list.append(mismatch_dict_list[index][seq])
            else:
                corrected_seq_list.append(seq)

        return bool_valid, bool_corrected, '_'.join(corrected_seq_list)

    @staticmethod
    def parse_whitelist_file(files: list, n_pattern: int, n_mismatch: int):
        """
        files: file paths
        n_pattern: number of sections in pattern
        n_mismatch: allowed number of mismatch bases
        Returns:
            white_set_list
            mismatch_list
        """
        n_files = len(files)
        if n_files == 1 and n_pattern > 1:
            files = [files[0]] * n_pattern
        elif n_files != n_pattern:
            sys.exit(f'number of whitelist files({n_files} files:{files}) != n_pattern({n_pattern})')
        
        white_set_list, mismatch_list = [], []
        for f in files:
            barcodes = utils.read_one_col(f)
            white_set_list.append(set(barcodes))
            barcode_mismatch_dict = Barcode.get_mismatch_dict(barcodes, n_mismatch)
            mismatch_list.append(barcode_mismatch_dict)

        return white_set_list, mismatch_list

    @staticmethod
    def parse_chemistry(chemistry):
        """
        Returns: pattern_dict, barcode_set_list, barcode_mismatch_list, linker_set_list, linker_mismatch_list
        """
        pattern = PATTERN_DICT[chemistry]
        pattern_dict = Barcode.parse_pattern(pattern)
        ######################################################################################3
        linker_file, whitelist_file = Barcode.get_scope_bc(chemistry, ASSETS_DIR)
        #whitelist_file = "/SGRNJ06/randd/USER/liuzihao/work/CeleScope-master/celescope/data/chemistry/scopeV2.2.1/bclist"
        #linker_file = "/SGRNJ06/randd/USER/liuzihao/work/CeleScope-master/celescope/data/chemistry/scopeV2.2.1/linker_4types"

        barcode_set_list, barcode_mismatch_list = Barcode.parse_whitelist_file([whitelist_file], n_pattern=len(pattern_dict['C']), n_mismatch=1)
        linker_set_list, linker_mismatch_list = Barcode.parse_whitelist_file([linker_file],n_pattern=1, n_mismatch=2)

        return pattern_dict, barcode_set_list, barcode_mismatch_list, linker_set_list, linker_mismatch_list

    @staticmethod
    def check_polyT(seq, pattern_dict, min_polyT_count=MIN_T):
        """
        Return:
            True if polyT is found
        """
        seq_polyT = Barcode.get_seq_str(seq, pattern_dict['T'])
        n_polyT_found = seq_polyT.count('T')
        if n_polyT_found >= min_polyT_count:
            return True
        return False

    def open_files(self):
        if self.output_R1:
            self.fh_fq1 = xopen(self.out_fq1, 'w')
        if not self.args.stdout:
            self.fh_fq2 = xopen(self.out_fq2, 'w')

        if self.nopolyT:
            self.fh_nopolyT_fq1 = xopen(self.nopolyT_1, 'w')
            self.fh_nopolyT_fq2 = xopen(self.nopolyT_2, 'w')

        if self.noLinker:
            self.fh_nolinker_fq1 = xopen(self.noLinker_1, 'w')
            self.fh_nolinker_fq2 = xopen(self.noLinker_2, 'w')

    def close_files(self):
        if self.output_R1:
            self.fh_fq1.close()
        if not self.args.stdout:
            self.fh_fq2.close()

        if self.nopolyT:
            self.fh_nopolyT_fq1.close()
            self.fh_nopolyT_fq2.close()
        
        if self.noLinker:
            self.fh_nolinker_fq1.close()
            self.fh_nolinker_fq2.close()

    def add_step_metrics(self):

        self.add_log_record(
            f'Raw Reads: {self.total_num}'
        )
        self.add_log_record(
            f'Valid Reads: {self.clean_num}({self.get_pct(self.clean_num, self.total_num)})'
        )

        BarcodesQ30 = sum([self.barcode_qual_Counter[k] for k in self.barcode_qual_Counter if k >= self.ord2chr(
            30)]) / float(sum(self.barcode_qual_Counter.values())) * 100
        BarcodesQ30 = round(BarcodesQ30, 2)
        BarcodesQ30_display = f'{BarcodesQ30}%'
        self.add_log_record(f"Q30 of Barcodes: {BarcodesQ30_display}")
        #self.add_metric(
        #    name='Q30 of Barcodes',
        #    value=BarcodesQ30,
        #    display=BarcodesQ30_display,
        #    help_info='percent of barcode base pairs with quality scores over Q30',
        #)

        UMIsQ30 = sum([self.umi_qual_Counter[k] for k in self.umi_qual_Counter if k >= self.ord2chr(
            30)]) / float(sum(self.umi_qual_Counter.values())) * 100
        UMIsQ30 = round(UMIsQ30, 2)
        UMIsQ30_display = f'{UMIsQ30}%'

        self.add_log_record(f"Q30 of UMIs: {UMIsQ30_display}")
        #self.add_metric(
        #    name='Q30 of UMIs',
        #    value=UMIsQ30,
        #    display=UMIsQ30_display,
        #    help_info='percent of UMI base pairs with quality scores over Q30',
        #)

        #self.add_metric(
        #    name='No PolyT Reads',
        #    value=self.no_polyT_num,
        #    total=self.total_num,
        #    show=False
        #)
        self.add_log_record(
            f'No PolyT Reads: {self.no_polyT_num}({self.get_pct(self.no_polyT_num, self.total_num)})'
        )
        self.add_log_record(
            f'Low Quality Reads: {self.lowQual_num}({self.get_pct(self.lowQual_num, self.total_num)})'
        )
        self.add_log_record(
            f'No Linker Reads: {self.no_linker_num}({self.get_pct(self.no_linker_num, self.total_num)})'
        )
        self.add_log_record(
            f'No Barcode Reads: {self.no_barcode_num}({self.get_pct(self.no_barcode_num, self.total_num)})'
        )

        self.add_log_record(
            f'Corrected Linker Reads: {self.linker_corrected_num}({self.get_pct(self.linker_corrected_num, self.total_num)})'
        )
        self.add_log_record(
            f'Corrected Barcode Reads: {self.barcode_corrected_num}({self.get_pct(self.barcode_corrected_num, self.total_num)})'
        )

        if self.clean_num == 0:
            raise Exception('no valid reads found! please check the --chemistry parameter.') #+ HELP_DICT['chemistry'])
        

    def run(self):
        """
        Extract barcode and UMI from R1. Filter reads with 
            - invalid polyT
            - low quality in barcode and UMI
            - invalid inlinker
            - invalid barcode
            
        for every sample
            get chemistry
            get linker_mismatch_dict and barcode_mismatch_dict
            for every read in read1
                filter
                write valid R2 read to file
        """

        for i in range(self.fq_number):

            chemistry = self.chemistry_list[i]
            lowNum = int(self.lowNum)
            lowQual = int(self.lowQual)
            
            if chemistry == 'scopeV1':
                lowNum = min(0, lowNum)
                lowQual = max(10, lowQual)
                #self.run.logger.info(f'scopeV1: lowNum={lowNum}, lowQual={lowQual} ')
            
            # get linker and whitelist 拿到linker和white list
            #################################################################################
            bc_pattern = PATTERN_DICT[chemistry]
            if (bc_pattern):
                linker_file, whitelist_file = self.get_scope_bc(chemistry)
                #whitelist_file = "/SGRNJ06/randd/USER/liuzihao/work/CeleScope-master/celescope/data/chemistry/scopeV2.2.1/bclist"
                #linker_file = "/SGRNJ06/randd/USER/liuzihao/work/CeleScope-master/celescope/data/chemistry/scopeV2.2.1/linker_4types"
                whitelist_files = [whitelist_file]
            else:
                bc_pattern = self.pattern
                linker_file = self.linker
                whitelist_file = self.whitelist
                whitelist_files = whitelist_file.split(',')
            if not bc_pattern:
                raise Exception("invalid bc_pattern!")

            # 解析pattern
            pattern_dict = self.parse_pattern(bc_pattern)

            bool_T = True if 'T' in pattern_dict else False
            bool_L = True if 'L' in pattern_dict else False
            bool_whitelist = (whitelist_file is not None) and whitelist_file != "None"
            C_len = sum([item[1] - item[0] for item in pattern_dict['C']])

            if bool_whitelist:
                barcode_set_list, barcode_mismatch_list = Barcode.parse_whitelist_file(whitelist_files,
                                                n_pattern=len(pattern_dict['C']), n_mismatch=1)
            if bool_L:
                linker_set_list, linker_mismatch_list = Barcode.parse_whitelist_file([linker_file],n_pattern=1, n_mismatch=2)

            with pysam.FastxFile(self.fq1_list[i], persist=False) as fq1, \
                    pysam.FastxFile(self.fq2_list[i], persist=False) as fq2:
                for entry1, entry2 in zip(fq1, fq2):
                    header1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                    header2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                    self.total_num += 1

                    # polyT filter
                    if bool_T and self.filterNoPolyT:
                        if not Barcode.check_polyT(seq1, pattern_dict):
                            self.no_polyT_num += 1
                            if self.nopolyT:
                                self.fh_nopolyT_fq1.write(
                                    '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                                self.fh_nopolyT_fq2.write(
                                    '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                            continue

                    # lowQual filter
                    C_U_quals_ascii = Barcode.get_seq_str(
                        qual1, pattern_dict['C'] + pattern_dict['U'])
                    if lowQual > 0 and Barcode.low_qual(C_U_quals_ascii, lowQual, lowNum):
                        self.lowQual_num += 1
                        continue

                    # linker filter
                    if bool_L and (not self.allowNoLinker):
                        seq_str = Barcode.get_seq_str(seq1, pattern_dict['L'])
                        bool_valid, bool_corrected, _ = Barcode.check_seq_mismatch(
                            [seq_str], linker_set_list, linker_mismatch_list)
                        if not bool_valid:
                            self.no_linker_num += 1
                            if self.noLinker:
                                self.fh_nolinker_fq1.write(
                                    f'@{header1}\n{seq1}\n{seq_str}\n{qual1}\n')
                                self.fh_nolinker_fq2.write(
                                    '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                            continue
                        elif bool_corrected:
                            self.linker_corrected_num += 1

                    # barcode filter
                    seq_list = self.get_seq_list(seq1, pattern_dict, 'C')

                    if bool_whitelist:
                        bool_valid, bool_corrected, corrected_seq = Barcode.check_seq_mismatch(
                            seq_list, barcode_set_list, barcode_mismatch_list)

                        if not bool_valid:
                            self.no_barcode_num += 1
                            continue
                        elif bool_corrected:
                            self.barcode_corrected_num += 1
                        cb = corrected_seq
                    else:
                        cb = "_".join(seq_list)

                    self.clean_num += 1
                    self.barcode_qual_Counter.update(C_U_quals_ascii[:C_len])
                    self.umi_qual_Counter.update(C_U_quals_ascii[C_len:])

                    umi = Barcode.get_seq_str(seq1, pattern_dict['U'])
                    if not umi:
                        continue

                    if self.args.stdout:
                        print(f'@{cb}:{umi}:{self.total_num}\n{seq2}\n+\n{qual2}')
                    else:
                        self.fh_fq2.write(f'@{cb}:{umi}:{self.total_num}\n{seq2}\n+\n{qual2}\n')
                    if self.output_R1:
                        self.fh_fq1.write(f'@{cb}:{umi}:{self.total_num}\n{seq1}\n+\n{qual1}\n')                   
            
                #self.run.logger.info(self.fq1_list[i] + ' finished.')

        self.close_files()
        self.add_step_metrics()
        self.write_log()



def main():
    parser = argparse.ArgumentParser(description='barcode')

    parser.add_argument(
        '--chemistry',
        help='Predefined (pattern, barcode whitelist, linker whitelist) combinations. ' ,#+ HELP_DICT['chemistry'],
        choices=list(PATTERN_DICT.keys()),
        default='auto'
    )
    
    parser.add_argument(
        '--pattern',
        help="""The pattern of R1 reads, e.g. `C8L16C8L16C8L1U12T18`. The number after the letter represents the number 
        of bases.  
- `C`: cell barcode  
- `L`: linker(common sequences)  
- `U`: UMI    
- `T`: poly T""",
    )
    
    parser.add_argument(
        '--whitelist',
        help='Cell barcode whitelist file path, one cell barcode per line.'
    )
    
    parser.add_argument(
        '--linker',
        help='Linker whitelist file path, one linker per line.'
    )
    
    parser.add_argument(
        '--lowQual',
        help='Default 0. Bases in cell barcode and UMI whose phred value are lower than \
lowQual will be regarded as low-quality bases.',
        type=int,
        default=0
    )
    
    parser.add_argument(
        '--lowNum',
        help='The maximum allowed lowQual bases in cell barcode and UMI.',
        type=int,
        default=2
    )
    
    parser.add_argument(
        '--nopolyT',
        help='Outputs R1 reads without polyT.',
        action='store_true',
    )
    
    parser.add_argument(
        '--noLinker',
        help='Outputs R1 reads without correct linker.',
        action='store_true',
    )
    parser.add_argument(
        '--filterNoPolyT',
        help="Filter reads without PolyT.",
        action='store_true'
    )
    parser.add_argument(
        '--allowNoLinker',
        help="Allow valid reads without correct linker.",
        action='store_true'
    )
    parser.add_argument(
        '--output_R1',
        help="Output valid R1 reads.",
        action='store_true'
    )

    parser.add_argument('--fq1', help='R1 fastq file. Multiple files are separated by comma.', required=True)
    parser.add_argument('--fq2', help='R2 fastq file. Multiple files are separated by comma.', required=True)
    
    parser.add_argument('--sample', help='', required=True)
    parser.add_argument('--outdir', help='', required=True)

    parser.add_argument(
        '--stdout',
        help="Write output to standard output",
        action='store_true'
)
    args = parser.parse_args()
    runner = Barcode(args)
    runner.run()

if __name__ == '__main__':
    main()
