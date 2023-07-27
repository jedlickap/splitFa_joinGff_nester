"""
Program takes nester output GFF and sort TEs within
file from 5 prime to 3 prime end in sequence
! Important: no TE numbering recalculation (yet) !
"""
# gff = 'Zmays_Chr10_2Mb_split/data/2Mb_part0_0_500000/2Mb_part0_0_500000_igv.gff'

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

def sort_gff(gff):
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
    # output sorted GFF
    out_name = gff.replace(".gff","_sorted.gff")
    with open(out_name, "w") as out:
        for my_teid in te_dict_sort:
            for line in te_dict_sort[my_teid]:
                out.write(line)
    return out_name

"""
['2Mb_part0_0_500000\tnested-nester\tnested_repeat\t475440\t481818\t.\t.\t.\tID=TE_BASE 0;annot=Class_I::LTR::Ty3/gypsy::chromovirus::CRM;name=TE_BASE 0;color=#a9a9a9\n',
 '2Mb_part0_0_500000\tnested-nester\trepeat_fragment\t475440\t481818\t.\t.\t.\tID=TE 0-0;Parent=TE_BASE 0;name=TE 0;color=#a9a9a9\n', 
 '2Mb_part0_0_500000\tnested-nester\tpolypeptide_conserved_region\t476055\t476832\t.\t-\t.\tID=DOMAIN 0-0-0;Parent=TE_BASE 0;annot=Class_I::LTR::Ty3/gypsy::chromovirus::CRM INT_crm;name=INT;color=#b10026\n',
 '2Mb_part0_0_500000\tnested-nester\tpolypeptide_conserved_region\t477299\t477650\t.\t-\t.\tID=DOMAIN 0-1-0;Parent=TE_BASE 0;annot=Class_I::LTR::Ty3/gypsy::chromovirus::CRM RNaseH_crm;name=RH;color=#fc4e2a\n',
 '2Mb_part0_0_500000\tnested-nester\tpolypeptide_conserved_region\t477845\t477965\t.\t-\t.\tID=DOMAIN 0-2-0;Parent=TE_BASE 0;annot=Class_I::LTR::Ty3/gypsy::chromovirus::CRM RT_crm;name=RT;color=#feb24c\n',
 '2Mb_part0_0_500000\tnested-nester\tpolypeptide_conserved_region\t478436\t478703\t.\t-\t.\tID=DOMAIN 0-3-0;Parent=TE_BASE 0;annot=Class_I::LTR::Ty3/gypsy::chromovirus::CRM AP_crm;name=PROT;color=#ffeda0\n', 
 '2Mb_part0_0_500000\tnested-nester\tpolypeptide_conserved_region\t478963\t479671\t.\t-\t.\tID=DOMAIN 0-4-0;Parent=TE_BASE 0;annot=Class_I::LTR::Ty3/gypsy::chromovirus::CRM GAG_crm;name=GAG;color=#225ea8\n', 
 '2Mb_part0_0_500000\tnested-nester\tprimer_binding_site\t481223\t481244\t.\t.\t.\tID=PBS 0;Parent=TE_BASE 0;name=ppt;color=#41ab5d\n', 
 '2Mb_part0_0_500000\tnested-nester\tRR_tract\t476008\t476022\t.\t.\t.\tID=PPT 0;Parent=TE_BASE 0;name=pbs;color=#c7e9c0\n', 
 '2Mb_part0_0_500000\tnested-nester\tlong_terminal_repeat\t481249\t481818\t.\t.\t.\tID=LTR LEFT 0-0;Parent=TE_BASE 0;name=ltr left;color=#f032e6\n',
 '2Mb_part0_0_500000\tnested-nester\tlong_terminal_repeat\t475440\t476008\t.\t.\t.\tID=LTR RIGHT 0-0;Parent=TE_BASE 0;name=ltr right;color=#f032e6\n',
 '2Mb_part0_0_500000\tnested-nester\ttarget_site_duplication\t475435\t475439\t.\t.\t.\tID=TSR RIGHT 0;Parent=TE_BASE 0;name=tsr right;color=#333333\n',
 '2Mb_part0_0_500000\tnested-nester\ttarget_site_duplication\t481819\t481823\t.\t.\t.\tID=TSR LEFT 0;Parent=TE_BASE 0;name=tsr left;color=#333333\n']
"""