from glob import glob
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# fa = sys.argv[1]

class FastaIn:

    def __init__(self, fasta):
        self.fasta = fasta
        
    def __splitted_fas_summary(self, splitFaLens, fa):
        tot_cnt = len(splitFaLens)
        cnt_10kbp = len([fa for fa in splitFaLens if fa < 10_000])
        cnt_50kbp = len([fa for fa in splitFaLens if fa < 50_000])
        perc_10kbp, perc_50kbp, perc_over50kbp = round((cnt_10kbp/tot_cnt)*100,2), round((cnt_50kbp/tot_cnt)*100,2), round(((tot_cnt-cnt_50kbp)/tot_cnt)*100,2)
        print(f"""Fasta: \'{fa}\' was splitted into {tot_cnt} sequences.
    {cnt_10kbp} ({perc_10kbp} %) < 10 kbp
    {cnt_50kbp} ({perc_50kbp} %) < 50 kbp
    {tot_cnt-cnt_50kbp} ({perc_over50kbp} %) > 50 kbp""")
        
    def __getPattIdxList(self, fa, pattern, writeFa):
        splitFaLens = []
        rl = []
        with open(fa.replace(".fa",f"_split_{pattern}_Indexes.tab"),"w") as out_b:
            for r in SeqIO.parse(fa,"fasta"):
                idx_l = [[int(m.start(0)), int(m.end(0))] for m in re.finditer(rf'[{pattern}]+',str(r.seq))]    # find all matches of given patern with start-end indexes
                cnt = 0
                for se in idx_l:
                    out_b.write(f"{r.id}\t{cnt}\t{se[0]}\t{se[1]}\t{se[1]-se[0]}\n")    # write into table
                    if writeFa:
                        rl.append(SeqRecord(r.seq[se[0]:se[1]],r.id+'_part_'+str(cnt),"",""))   # split sequence into multiple SeqRecord objecs
                        splitFaLens.append(se[1]-se[0])
                    cnt+=1               
        if rl:
            out_fa_name = fa.replace(".fa","_splitted.fa")
            handle = open(out_fa_name,"w")
            SeqIO.write(rl, handle, "fasta")    # write down splitted fasta file
            self.__splitted_fas_summary(splitFaLens, fa)

    def split_fa_by_Ns(self):
        for k, v in {'ACGT':True,'nN':False}.items():
            self.__getPattIdxList(self.fasta, k, v)

class NesterOut():

    def __init__(self, data_path, indexes_tab):
        self.data_path = data_path
        self.indexes_tab = indexes_tab

    def __firstIdxDict(self):
        d = {}
        with open(self.indexes_tab) as acgt:
            for l in acgt:
                ll = l.rstrip().split("\t")
                d[ll[1]] = int(ll[2])
        return d

    def recalc_and_join(self):
        gff_list = glob(f"{self.data_path}/*/*cols.gff")
        out_name = self.indexes_tab.replace("_split_ACGT_Indexes.tab","_recalc_joined.gff")
        d = self.__firstIdxDict()
        with open(out_name,"w") as out:
            for gff in gff_list:
                orig_ch, splitPart = gff.split("/")[-2].split("_")[0], gff.split("/")[-2].split("_")[2]
                with open(gff) as gf:
                    for line in gf:
                        if line.startswith(orig_ch):
                            ll = line.split("\t")
                            ll[0] = orig_ch
                            ll[3], ll[4] = str(int(ll[3])+int(d[splitPart])),str(int(ll[4])+int(d[splitPart]))
                            new_line = "\t".join(ll)
                            out.write(new_line)
