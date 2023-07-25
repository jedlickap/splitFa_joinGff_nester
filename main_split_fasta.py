import sys
from split_fasta import Fasta

faObj = Fasta(sys.argv[1], 4)
faObj.split_fasta()