import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class Fasta:
    
    def __init__(self, fasta, subseq_nr = 2, overlap = 100_000):
        self.fasta = fasta
        self.subseq_nr = subseq_nr
        self.overlap = overlap
    
    def __generate_seq_rec(self, seq_object, start, end, subseq_cnt):
        return SeqRecord(seq_object.seq[start:end],
                         seq_object.id +
                         "_part" +
                         str(subseq_cnt) +
                         f"_{str(start)}_{str(end)}",
                         "",
                         "")
    
    def __split_seq(self, seq_obj, chunk_size, overlap_size, seq_len):
        rl = []
        subseq_cnt = 0
        start, end = 0, chunk_size

        while (seq_len-end > chunk_size):   # all the chunk+overlap length subsequences
            rl.append(self.__generate_seq_rec(seq_obj, start, end, subseq_cnt))
            subseq_cnt += 1
            start = end - overlap_size
            end += chunk_size
        rl.append(self.__generate_seq_rec(seq_obj, start, seq_len, subseq_cnt)) # sequence tail 
        return rl
    
    def split_fasta(self):
        rec_list = []
        for rec in SeqIO.parse(self.fasta, 'fasta'):
            seq_len = len(rec.seq)
            if seq_len / self.subseq_nr < 100_000:
                rec_list.append(rec)    # write_down sequence
            else:
                # run splitting
                chunk_size = int(seq_len / self.subseq_nr)
                rec_list += self.__split_seq(rec, chunk_size, self.overlap, seq_len)
        # output directory
        fa_pref = ".".join(self.fasta.split(".")[:-1])
        outdir = fa_pref + '_split/'
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # generate splitted fasta
        out_fa_path = outdir + fa_pref + "_split.fa"
        handle = open(out_fa_path,"w")
        SeqIO.write(rec_list,handle,"fasta")
        # return out_fa_path