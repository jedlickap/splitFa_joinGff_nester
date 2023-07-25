import sys
from split_fasta import Fasta

faObj = Fasta(sys.argv[1], 6)
faObj.split_fasta()