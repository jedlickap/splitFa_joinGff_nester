import sys
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

fa = sys.argv[1]

def splitted_fas_summary(splitFaLens, fa):
    tot_cnt = len(splitFaLens)
    cnt_10kbp = len([fa for fa in splitFaLens if fa < 10_000])
    cnt_50kbp = len([fa for fa in splitFaLens if fa < 50_000])
    perc_10kbp, perc_50kbp, perc_over50kbp = round((cnt_10kbp/tot_cnt)*100,2), round((cnt_50kbp/tot_cnt)*100,2), round(((tot_cnt-cnt_50kbp)/tot_cnt)*100,2)
    print(f"""Fasta: \'{fa}\' was splitted into {tot_cnt} sequences.
{cnt_10kbp} ({perc_10kbp} %) < 10 kbp
{cnt_50kbp} ({perc_50kbp} %) < 50 kbp
{tot_cnt-cnt_50kbp} ({perc_over50kbp} %) > 50 kbp""")
    
def getPattIdxList(fa, pattern, writeFa):
    splitFaLens = []
    rl = []
    for r in SeqIO.parse(fa,"fasta"):
        idx_l = [[int(m.start(0)), int(m.end(0))] for m in re.finditer(rf'[{pattern}]+',str(r.seq))]    # find all matches of given patern with start-end indexes
        with open(fa.replace(".fa",f"_split_{pattern}_Indexes.tab"),"w") as out_b:
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
        splitted_fas_summary(splitFaLens, fa)

def main():
    for k, v in {'ACGT':True,'nN':False}.items():
        getPattIdxList(fa, k, v)

if __name__=='__main__':
    main()