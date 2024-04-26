import gzip
import glob
import subprocess
import multiprocessing

def openfile(file_name, mode='rt', **kwargs):
    """open gzip or plain file"""
    if file_name.endswith('.gz'):
        file_obj = gzip.open(file_name, mode=mode, **kwargs)
    else:
        file_obj = open(file_name, mode=mode, **kwargs)
    return file_obj


def read_one_col(fn):
    """read one column file into list"""
    with openfile(fn) as f:
        return [x.strip() for x in f]


def get_matrix_file_path(match_dir, sample_name = None):
    """
    compatible with non-gzip file
    """         
                # celescope match dir: 
                # line1-3: old version
                # line4: new version
    mtx_p =     glob.glob(f"{match_dir}/*count/*_matrix_10X/barcode*") + \
                glob.glob(f"{match_dir}/*count/*_filtered_feature_bc_matrix/barcode*") + \
                glob.glob(f"{match_dir}/outs/filtered/barcode*") + \
                glob.glob(f"{match_dir}/starsolo/{sample_name}.matrix/filtered/barcode*")
        
    if len(mtx_p) > 1:
        raise FileNotFoundError(f"Multi mtx dir {mtx_p}")

    if len(mtx_p) == 0:
        raise FileNotFoundError(f"Not found mtx dir {mtx_p}")

    return  mtx_p[0]


def get_barcode_from_matrix_dir(matrix_dir):
    """
    Returns:
        match_barcode: list
        no_match_barcode: int
    """
  
    match_barcode_file = get_matrix_file_path(matrix_dir)
    match_barcode = read_one_col(match_barcode_file)
    n_match_barcode = len(match_barcode)

    return match_barcode, n_match_barcode

def samtools(inputbam, do, outputbam = None, sort_by = "pos"):
    if do == "index":
        cmd = f"samtools index {inputbam}"
    
    if do == "sort":
        cpus = min(multiprocessing.cpu_count(), 16)
        cmd = f"samtools sort {inputbam} -o {outputbam} --threads {cpus}"
        if sort_by == "name":
            cmd += " -n"
    run_cmd(cmd)

def get_bed(panel):
    pass

def get_fasta_from_index_dir(index_dir):
    pass



def run_cmd(cmd):
    subprocess.check_call(f"{cmd}", shell=True)