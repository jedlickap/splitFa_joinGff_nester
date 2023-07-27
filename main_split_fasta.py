import sys
import os
import subprocess
from split_fasta import Fasta
from join_gffs import join_gffs

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

# 1. split FASTA
faObj = Fasta(sys.argv[1], 4)
faObj.split_fasta()
fa_path = faObj.split_fasta()

# 2. run nester
"""
prim_path = os.getcwd()
fa_folder, fa_name = fa_path.split("/")[0], fa_path.split("/")[1]
fa_pwd = prim_path + "/" + fa_folder
os.chdir(fa_pwd)
cmd = f"nested-nester {fa_name}"
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
process.wait()
print(f"COMMAND: {cmd}: ", process.returncode)
os.chdir(prim_path)
"""

# 3. sort, recalculate coordinates and join GFFs
fa_path = "Zmays_Chr10_2Mb_split/Zmays_Chr10_2Mb_split.fa"
join_gffs(fa_path)
