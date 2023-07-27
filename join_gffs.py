import os
from file_read_backwards import FileReadBackwards
from Bio import SeqIO

fa_path = "Zmays_Chr10_2Mb_split/Zmays_Chr10_2Mb_split.fa"

fa_heads_list = [r.id for r in SeqIO.parse(fa_path,'fasta')]

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

def gff_paths2dict(fa_heads_list):
    # add paths of all respetive gff into dict
    path_dict = {'igv':[],'trf':[],'ltrs':[],'mite':[]}
    for p in fa_heads_list:
        split_folder = fa_path.split("/")[0]
        for gff in list(os.listdir(split_folder+'/data/'+p+'/')):
            if '.gff' in gff:
                full_gff_path = split_folder+'/data/'+p+'/'+gff
                path_dict[gff.split("_")[-1].replace(".gff","")].append(full_gff_path)
    return path_dict

# recalc coordinates for igv.gffs
repeat_key = 'igv'
print(gff_paths2dict(fa_heads_list))
# def last_rec(gff_path, repeat_type):
#     repeat_type = 'nested_repeat'
#     # find last 'nested_nester' record in the first gff
#     with FileReadBackwards(gff_path, encoding="utf-8") as gff:
#         for line in gff:
#             if f'\t{repeat_type}\t' in line:
#                 return line