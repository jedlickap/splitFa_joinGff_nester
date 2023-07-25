import sys
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def args_from_parser():
    parser = argparse.ArgumentParser(
        description='''Script performs split of multisequence fasta file
        based on defined number of subsequences and overlap size
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

def split_seq(seq_obj, chunk_size, overlap_size, seq_len):
    rl = []
    subseq_cnt = 0
    start, end = 0, chunk_size

    # 1st seq
    rl.append(SeqRecord(seq_obj.seq[start:end],seq_obj.id+"_part"+str(subseq_cnt)+f":{str(start)}|{str(end)}","",""))
    subseq_cnt += 1
    start = end - overlap_size
    end += chunk_size
    while (seq_len-end > chunk_size):
        rl.append(SeqRecord(seq_obj.seq[start:end],seq_obj.id+"_part"+str(subseq_cnt)+f":{str(start)}|{str(end)}","",""))
        subseq_cnt += 1
        start = end - overlap_size
        end += chunk_size
    rl.append(SeqRecord(seq_obj.seq[start:seq_len],seq_obj.id+"_part"+str(subseq_cnt)+f":{str(start)}|{str(end)}","",""))
    return rl

def main():
    args = args_from_parser()
    fa, subseq_nr, overlap_size = args.fasta, args.number_of_subseqs, args.subseq_overlap
    rec_list = []
    for rec in SeqIO.parse(fa, 'fasta'):
        seq_len = len(rec.seq)
        if seq_len / subseq_nr < 100_000:
            rec_list.append(rec)    # write_down sequence
        else:
            # run splitting
            chunk_size = int(seq_len / subseq_nr)
            rec_list += split_seq(rec, chunk_size, overlap_size, seq_len)
    # output directory
    fa_pref = ".".join(fa.split(".")[:-1])
    outdir = fa_pref + '_split/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # generate splitted fasta
    out_fa_path = outdir + fa_pref + "_split.fa"
    handle = open(out_fa_path,"w")
    SeqIO.write(rec_list,handle,"fasta")            

if __name__=='__main__':
    main()