#!/usr/bin/python
'''
This script is meant to bed files associated with the elements based on a grouping
such as class, order, superfamily, and family. 

By Maia whos trying O.o

USAGE: python breakup_intoBed.py TEclass data_breakup_15.txt
creates a 3 files: data_breakup_15.txt1.bed data_breakup_15.txt2.bed    data_breakup_15.txtUnknown.bed
'''
import re
import fileinput
from sys import argv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict



order_labels=["LTR","TIR","SINE","LINE","Helitron","Unknown"]
class_labels=["1","2","Unknown"]
superfamily_labels = ["LTR/Gypsy","LTR/Copia","DNA/HAT","DNA/Pogo","DNA/En-Spm","DNA/MuDR","DNA","DNA/Tc1","DNA/Harbinger","DNA/Mariner","LINE?","LINE/L1","SINE","RC/Helitron","RathE3_cons","RathE2_cons","RathE1_cons","Unassigned"]
family_labels =[] #fill in later

#   orders: LTR, TIR, SINE, LINE, Helitron, Unknown
TEorder_dict={"LTR/Gypsy": "LTR","LTR/Copia" : "LTR","DNA/HAT" : "TIR","DNA/Pogo" : "TIR","DNA/En-Spm" : "TIR","DNA/MuDR" : "TIR","DNA" : "TIR","DNA/Tc1" : "TIR","DNA/Harbinger" : "TIR","DNA/Mariner" : "TIR", "LINE?" : "LINE","LINE/L1" : "LINE","SINE" : "SINE","RC/Helitron" : "Helitron","RathE3_cons" : "Unknown","RathE2_cons" : "Unknown","RathE1_cons" : "Unknown","Unassigned" : "Unknown"}
#   classes: 1) retrotransposon, and 2) DNA transposon, and Unknown
TEclass_dict = {"LTR" : "1","Helitron" : "2","TIR" : "2","LINE" : "1","SINE" : "1","Unknown" : "Unknown"}

TE_level = argv[1]
knownTE_file = argv[2]

 #TE_level can be TEclass, TEorder, TEfamily, TEsuperfamily
#I need this to output a list for each TE_level
groups = ''
if (TE_level == 'TEclass'):
    groups = class_labels
    print('got the class!')
elif (TE_level == 'TEorder'):
    groups = order_labels
#elif(TE_level == TEfamily):
    #groups = 
elif(TE_level == 'TEsuperfamily'):
    groups = superfamily_labels
else:
    print('Invalid TE level')

out_list=['nada']

infile = open(knownTE_file,'r')
header=infile.readline().rstrip().split("\t")  #   read in a single line from the file pointer

def new_group_file(path):
	f = open(path, "w")
	f.write('\t'.join(['#chromosome','start','stop','name']))
	return f
group_files = dict(zip(groups, [new_group_file(knownTE_file + "_" + TE_level.replace("TE", "") + "_" + g.replace("/","_") + ".bed") for g in groups]))

for line in infile:

    name,orientation5prime,start,end,family,superFamily=line.rstrip().split("\t")
    chromosome="chr"+re.search("^AT[0-9+]",name).group().replace("AT","")
    
    #not all the information is in the column of that file (the order and class)
    Seg_TEorder = TEorder_dict[superFamily]
    Seg_TEclass = TEclass_dict[Seg_TEorder]

    for each in groups: # so if we're looking at classes, this would be 1,2,unknown
        if (TE_level == 'TEclass' and Seg_TEclass == each):
            group_files[each].write('\t'.join([chromosome,start,end,name])+'\n')
        elif(TE_level == 'TEorder' and Seg_TEorder == each):
            group_files[each].write('\t'.join([chromosome,start,end,name])+'\n')
        elif(TE_level == 'TEsuperfamily' and superFamily == each):
            group_files[each].write('\t'.join([chromosome,start,end,name])+'\n')
            