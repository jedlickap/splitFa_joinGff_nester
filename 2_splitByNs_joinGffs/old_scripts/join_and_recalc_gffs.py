import sys
from glob import glob

nester_data_dir = sys.argv[1]
acgt_idxs = sys.argv[2]  # Zmays_Chr10_2Mb_split_ACGT_Indexes.tab

gff_list = glob("split_nesterRun/data/*/*cols.gff")
#gff_list = glob.glob(f"{nester_data_dir}/*/*cols.gff")

def firstIdxDict(acgt_idxs):
    d = {}
    with open(acgt_idxs) as acgt:
        for l in acgt:
            ll = l.rstrip().split("\t")
            d[ll[1]] = int(ll[2])
    return d

def main():
    out_name = acgt_idxs.replace("_split_ACGT_Indexes.tab","_recalc_joined.gff")
    d = firstIdxDict(acgt_idxs)
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
                        
if __name__=='__main__':
    main()