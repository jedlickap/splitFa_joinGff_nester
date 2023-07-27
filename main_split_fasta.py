import sys
from split_fasta import Fasta

"""
argparse:
subseq_nr = 2
overlap = 100_00

Steps:
1. split FASTA
2. Run Nester
3. Sort *cols.gff
4. Join GFFs
"""

faObj = Fasta(sys.argv[1], 4)
faObj.split_fasta()