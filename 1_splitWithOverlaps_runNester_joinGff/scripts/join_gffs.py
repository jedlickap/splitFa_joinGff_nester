import os
from Bio import SeqIO

"""
Program takes nester output GFF and (i) sorts TEs within
file from 5 prime to 3 prime end in sequence;
(ii) recalculates coordinates; ane (iii) join them
into one GFF.
! Important: no TE numbering recalculation (yet) !
"""

class GFFLine:
    """
    Expect just one track line from GFF file:
    track = 'NC_041789.1	Gnomon	CDS	18512	18884	.	-	0	ID=cds-XP_028955055.1;Parent=rna-XM_029099222.1;Dbxref=GeneID:103406057,GenBank:XP_028955055.1;Name=XP_028955055.1;gbkey=CDS;gene=LOC103406057;product=uncharacterized protein LOC103406057;protein_id=XP_028955055.1'
    gfl = GFFLine_v2(track)
    """
    def __init__(self, track):
        self.__split_track(track)
        self.__attr_dict(self.attributes)
        
    def __split_track(self, track):
        track_l = track.rstrip().split("\t")
        self.seqid, self.source, self.type, self.start, self.end, self.score, self.strand, self.phase, self.attributes = track.rstrip().split("\t")
    
    def __attr_dict(self,attributes):
        if ";" in attributes and "=" in attributes:
            self.attributes_dict = {rec.split("=")[0]:rec.split("=")[1] for rec in attributes.split(";")}

def gff_paths2list(fa_heads_list, fa_path):
    # add paths of all respetive gff into list
    gff_list = []
    for p in fa_heads_list:
        split_folder = fa_path.split("/")[0]
        for gff in list(os.listdir(split_folder+'/data/'+p+'/')):
            if 'cols.gff' in gff:
                full_gff_path = split_folder+'/data/'+p+'/'+gff
                gff_list.append(full_gff_path)
    return gff_list

def sort_recalc_gff(gff):
    # divide GFF into TE specific lines
    seqid = ""
    teid = ""
    my_teid = ""
    gffObj = lambda: None
    te_dict = {}
    with open(gff) as gf:
        for line in gf:
            if '\tnested_repeat\t' in line:
                # print(line)
                gffObj = GFFLine(line)
                my_teid = f"{gffObj.seqid}|{gffObj.start}|{gffObj.end}"
                teid = gffObj.attributes_dict['ID']
                seqid = gffObj.seqid
                # print(teid,'\n',gffObj.attributes_dict)
                te_dict[my_teid] = []
                te_dict[my_teid].append(line)
            elif seqid + "\t" in line and f"Parent={teid};" in line:
                te_dict[my_teid].append(line)
    # sort dict
    te_dict_sort = dict(sorted(te_dict.items(), key = lambda x: int(x[0].split("|")[1])))
    # output sorted GFF and recalc coord
    out_name = gff.replace(".gff","_sorted.gff")
    with open(out_name, "w") as out:
        for my_teid in te_dict_sort:
            bps2add = int(my_teid.split("|")[0].split("_")[-2])
            for line in te_dict_sort[my_teid]:
                ll = line.split("\t")
                ll[0], ll[3], ll[4] = ll[0].split("_")[0], str(int(ll[3]) + bps2add), str(int(ll[4]) + bps2add)
                out.write("\t".join(ll))                
    return out_name

def sort_gffs(gff_list):
    sorted_gffs = []
    for gff in gff_list:
        sorted_gffs.append(sort_recalc_gff(gff))
    return sorted_gffs

def gff2teDict(gff):
    gffObj = lambda: None
    te_dict = {}
    with open(gff) as gf:
        for line in gf:
            if '\tnested_repeat\t' in line:
                gffObj = GFFLine(line)
                my_teid = f"{gffObj.seqid}|{gffObj.start}|{gffObj.end}"
                teid = gffObj.attributes_dict['ID']
                seqid = gffObj.seqid
                te_dict[my_teid] = []
                te_dict[my_teid].append(line)
            elif seqid + "\t" in line and f"Parent={teid};" in line:
                te_dict[my_teid].append(line)
    return te_dict

def find_intersections(aList, bList):
    out_list = []
    for teid in aList:
        if teid in bList: 
            out_list.append(teid)
    return out_list

def join_gffs(fa_path):
    """
    fa_path = "Zmays_Chr10_2Mb_split/Zmays_Chr10_2Mb_split.fa"
    fa_heads_list = [r.id for r in SeqIO.parse(fa_path,'fasta')]
    """    
    fa_heads_list = [r.id for r in SeqIO.parse(fa_path,'fasta')]
    unsorted_gffs = gff_paths2list(fa_heads_list, fa_path)
    sorted_gffs = sort_gffs(unsorted_gffs)
    
    # sort and recalculate gff coordinates
    teDictOut = {}
    for i, gff in enumerate(sorted_gffs[:-1]):
        if not teDictOut:
            aTeDict, bTeDict = gff2teDict(gff), gff2teDict(sorted_gffs[i+1])
            if find_intersections([k for k in aTeDict],[l for l in bTeDict]):
                stop_teid = find_intersections([k for k in aTeDict],[l for l in bTeDict])[0]
                aStopIndex, bStopIndex = list(aTeDict).index(stop_teid), list(bTeDict).index(stop_teid)
                aTeKeys, bTeKeys = [ak for ak in aTeDict][:aStopIndex],[bk for bk in bTeDict][bStopIndex:]
                teDictOut.update({k:aTeDict[k] for k in aTeDict if k in aTeKeys})
                teDictOut.update({k:bTeDict[k] for k in bTeDict if k in bTeKeys})
            else:
                teDictOut.update(aTeDict)
                teDictOut.update(bTeDict)
        else:
            bTeDict = gff2teDict(sorted_gffs[i+1])
            if find_intersections([k for k in teDictOut],[l for l in bTeDict]):
                stop_teid = find_intersections([k for k in teDictOut],[l for l in bTeDict])[0]
                outStopIndex, bStopIndex = list(teDictOut).index(stop_teid), list(bTeDict).index(stop_teid)
                outTeKeys, bTeKeys = [ak for ak in teDictOut][:outStopIndex],[bk for bk in bTeDict][bStopIndex:]
                teDictOut = {k:teDictOut[k] for k in teDictOut if k in outTeKeys}
                teDictOut.update({k:bTeDict[k] for k in bTeDict if k in bTeKeys})
            else:
                teDictOut.update(bTeDict)
    # write into a file            
    out_joinedGff = fa_path.split(".")[0] + "_joined.gff"
    with open(out_joinedGff,"w") as out:
        for k in teDictOut:
            for line in teDictOut[k]:
                out.write(line)
