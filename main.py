import sys
import os
import subprocess
import argparse
from scripts.split_fasta import Fasta
from scripts.join_gffs import join_gffs

def args_from_parser():
    parser = argparse.ArgumentParser(
        description='''Script performs:\n(i) split of multisequence fasta file
        based on defined number of subsequences and overlap size\n
        (ii) TE-greedy-nester run on splitted sequence\n
        (iii) Sort, coordinates recarcultions and joining into one GFF for original unspliced sequence 
        ''')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        "-fa", "--fasta", type=str, required=True,
        help='Input FASTA file'
        )
    parser.add_argument('-subseq_nr', '--number_of_subseqs', type=int,
                        default=2)
    parser.add_argument('-overlap', '--subseq_overlap', type=int,
                        default=100_000)
    return parser.parse_args()
"""
argparse:
subseq_nr = 2
overlap = 100_00

Steps:
1. split FASTA -> fa_path = "Zmays_Chr10_2Mb_split/Zmays_Chr10_2Mb_split.fa"
2. Run Nester -> nester_output = fa_path.split("/")[0] + "/data"
3. Sort *cols.gff -> gff_sorted NOTE: could be run within join step
4. Join GFFs -> '*joined.gff'
"""
def main():
    # 1. split FASTA
    faObj = Fasta(sys.argv[1], 4)
    faObj.split_fasta()
    fa_path = faObj.split_fasta()

    # 2. run nester
    prim_path = os.getcwd()
    fa_folder, fa_name = fa_path.split("/")[0], fa_path.split("/")[1]
    fa_pwd = prim_path + "/" + fa_folder
    os.chdir(fa_pwd)
    cmd = f"nested-nester {fa_name}"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    nester_procRetCode = process.returncode
    print(f"COMMAND: {cmd}: ", nester_procRetCode)
    os.chdir(prim_path)

    # 3. sort, recalculate coordinates and join GFFs
    if nester_procRetCode == 0:  
        join_gffs(fa_path)
    
if __name__=='__main__':
    main()