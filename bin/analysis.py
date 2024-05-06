import pysam
import pandas as pd

import pymodules.utils as utils
import argparse
from pymodules.minana import MinAna
from pymodules.snpeffer import SNPEff


AA_DICT = {
    'Gly' : 'G',
    'Ala' : 'A',
    'Val' : 'V',
    'Leu' : 'L',
    'Ile' : 'I',
    'Phe' : 'F',
    'Trp' : 'W',
    'Tyr' : 'Y',
    'Asp' : 'D',
    'Asn' : 'N',
    'Glu' : 'E',
    'Lys' : 'K',
    'Gln' : 'Q',
    'Met' : 'M',
    'Ser' : 'S',
    'Thr' : 'T',
    'Cys' : 'C',
    'Pro' : 'P',
    'His' : 'H',
    'Arg' : 'R',
}


def parse_variant_ann(variant_ann_file):
    """
    Args:
        variant_ann_file: variant annotation file from snpEff.
    
    Returns:
        gene_list, mRNA_list, protein_list
    """
    gene_list, mRNA_list, protein_list = [], [], []

    with open(variant_ann_file) as f:
        for line in f.readlines():
            if not line.startswith("#"):
                info = line.split('\t')[7]
                anns = info.split("|")
                gene = anns[3]
                gene_list.append(gene)
            
                tmp1, tmp2 = [], []
                for ann in anns:
                    if ann.startswith("c."):
                        exon_loc = anns[anns.index(ann) - 1].split('/')[0]
                        # WARNING_TRANSCRIPT_INCOMPLETE
                        if not exon_loc:
                            continue
                        
                        exon = ann.strip("c.")
                        exon = f"exon{exon_loc}:{exon}"
                        if exon not in tmp1:
                            tmp1.append(exon)

                    if ann.startswith("p."):
                        protein = ann[2:]
                        for i in AA_DICT:
                            protein = protein.replace(i, AA_DICT[i])
                        if protein not in tmp2:
                            tmp2.append(protein)
                        
                mRNA_list.append(','.join(tmp1))
                protein_list.append(','.join(tmp2))

    return (gene_list, mRNA_list, protein_list)


def parse_vcf_to_df(vcf_file, cols=('chrom', 'pos', 'alleles'), infos=('VID', 'CID')):
    """
    Read cols and infos into pandas df
    """
    vcf = pysam.VariantFile(vcf_file)
    df = pd.DataFrame(columns=[col.capitalize() for col in cols] + infos)
    rec_dict = {}
    for rec in vcf.fetch():

        for col in cols:
            rec_dict[col.capitalize()] = getattr(rec, col)
            if col == 'alleles':
                rec_dict['Alleles'] = '-'.join(rec_dict['Alleles'])

        for info in infos:
            rec_dict[info] = rec.info[info]

        '''
        rec_dict['GT'] = [s['GT'] for s in rec.samples.values()][0]
        rec_dict['GT'] = [str(item) for item in rec_dict['GT']]
        rec_dict['GT'] = '/'.join(rec_dict['GT'])
        '''
        df_new = pd.DataFrame(rec_dict, index=[0])
        df = pd.concat([df, df_new])

    vcf.close()
    df.reset_index(drop=True, inplace=True)
    return df

def vcf_to_gt_csv(vcf_file, csv_file):
    vcf = pysam.VariantFile(vcf_file)
    
    samples = vcf.header.samples
    
    with open(csv_file, 'w') as f:
        header = ['variant'] + list(samples)
        f.write(','.join(header) + '\n')
        
        for record in vcf:
            mutation_name = f"{record.chrom}_{record.pos}"
            genotypes = []
            
            for sample in samples:
                genotype = record.samples[sample]['GT']
                g1, g2 = genotype
                
                if g1 is None:
                    genotype_str = "NA"
                else:
                    genotype_str = '/'.join([str(g1),str(g2)])
                
                genotypes.append(genotype_str)
            
            line = [mutation_name] + genotypes
            f.write(','.join(line) + '\n')


class AnalysiSNP(MinAna):
    """
    ## Features
    - Annotate variants with [snpEff](http://pcingola.github.io/SnpEff/).

    ## Output
    - `{sample}_gt.csv` Genotypes of variants of each cell. Rows are variants and columns are cells.
    - `{sample}_variant_ncell.csv` Number of cells with each genotype.
    - `{sample}_variant_table.csv` annotated with snpEff.

    """

    def __init__(self, args):
        super().__init__(args.sample, args.outdir)
        self.vcf_file = args.input_vcf
        self.gene_list = utils.read_one_col(args.gene_list)
        self.n_gene = len(self.gene_list)
        self.database = args.database
        self.variant_table = None

        # out
        self.snpeff_ann_vcf_file = None
        self.final_vcf_file = f'{self.outdir}/{self.sample}_final.vcf'
        self.plot_snp_dir = f'{self.outdir}/{self.sample}_plot_snp/'

        self.gt_file = f'{self.outdir}/{self.sample}_gt.csv'
        self.ncell_file = f'{self.outdir}/{self.sample}_variant_ncell.csv'
        self.variant_table_file = f'{self.outdir}/{self.sample}_variant_table.csv'

    def write_gt(self):
        vcf_to_gt_csv(self.final_vcf_file, self.gt_file)

    def write_ncell(self):
        """
        parse gt_file to collect each genotype cell count into ncell_file
        """
        df = pd.read_csv(self.gt_file, index_col=0)
        df_ncell = df.apply(pd.Series.value_counts, axis=1).fillna(0).astype(int)
        df_ncell.to_csv(self.ncell_file, index=True)

    def ann(self):
        runner = SNPEff(
            self.database,
            self.vcf_file,
            self.outdir
        )
        runner.run_snpEff()
        self.snpeff_ann_vcf_file = runner.variants_vcf

    def keep_in_gene(self):
        """
        Output:
            self.final_vcf_file
        """
        gene_list, _, _ = parse_variant_ann(self.snpeff_ann_vcf_file)
        with pysam.VariantFile(self.snpeff_ann_vcf_file) as vcf_in:
            with pysam.VariantFile(self.final_vcf_file, 'w', header=vcf_in.header) as vcf_out:
                for i, record in enumerate(vcf_in.fetch()):
                    if gene_list[i] in self.gene_list:
                        vcf_out.write(record)              


    def get_variant_table(self):
        """
        Returns:
            is_in_gene_list: if res[i] == True, line i is in gene_list
        """

        df_vcf = parse_vcf_to_df(self.final_vcf_file, infos=[])
        df_vcf["Gene"], df_vcf["mRNA"], df_vcf["Protein"] =  parse_variant_ann(self.final_vcf_file)
        df_ncell = pd.read_csv(self.ncell_file)
        df_vcf = pd.concat([df_vcf, df_ncell], axis=1)

        cols = ["Chrom", "Pos", "Alleles", "Gene", "0/0", "0/1", "1/1", "mRNA", "Protein"]
        cols = [col for col in cols if col in df_vcf.columns]
        df_vcf = df_vcf.loc[:, cols]
        is_in_gene_list = df_vcf.Gene.isin(self.gene_list)
        df_vcf = df_vcf[is_in_gene_list]

        self.variant_table = df_vcf
        self.variant_table.reset_index(drop=True, inplace=True)
        self.variant_table.to_csv(self.variant_table_file, index=False)


    def run(self):
        self.ann()
        self.keep_in_gene()
        self.write_gt()
        self.write_ncell()
        self.get_variant_table()



def main():
    parser = argparse.ArgumentParser(description='analysis SNP')
    parser.add_argument('--input_vcf', help='vcf file.', required=True)
    parser.add_argument('--sample', help='sample name', required=True)
    parser.add_argument('--outdir', help='outdir', required=True)
    parser.add_argument("--gene_list", help='', required=True)
    parser.add_argument("--database", help='snpEff database. Common choices are GRCh38.mane.1.0.ensembl(human) and GRCm38.99(mouse)', 
                        default='GRCh38.mane.1.0.ensembl')
    args = parser.parse_args()

    runner = AnalysiSNP(args)
    runner.run()
    

if __name__ == '__main__':
    main()