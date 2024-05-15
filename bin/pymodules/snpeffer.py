import os
import subprocess


class SNPEff:
    def __init__(self,
                 database,
                 vcf_file,
                 outdir) -> None:
        
        self.database= database
        self.vcf_file = vcf_file
        self.outdir = outdir
        self.snpeff_outdir = ""

        self.out_vcf = "variants_ann.vcf"
        self.snpeff_genes = "./snpEff_genes.txt"
        self.summary = "./snpEff_summary.html"
    
    def run_snpEff(self):
        cmd = (
            f"snpEff -Xmx8g -v {self.database} {os.path.abspath(self.vcf_file)} >  {self.out_vcf}"
        )
        subprocess.check_call(cmd, shell=True)
        os.system(f"mkdir -p {self.outdir}/snpEff_out/")
        
        for file in [self.out_vcf, self.snpeff_genes, self.summary]:
            os.system(f"mv {file} {self.outdir}/snpEff_out/")

    @property
    def variants_vcf(self):
        return os.path.abspath(f"{self.outdir}/snpEff_out/{self.out_vcf}")
    
    def make_database(self, fa, gtf):
        pass
    
def main():
    database = "GRCh38.99"
    vcf_file = "/SGRNJ06/randd/USER/liuzihao/work/scsnp/bin/bin_testdata/celescope_test_script/snp/test1/08.filter_snp/test1_filtered.vcf"

    test = SNPEff(database,vcf_file, outdir="./out")
    test.run_snpEff()
    print(test.variants_vcf)

if  __name__ == '__main__':
    main()