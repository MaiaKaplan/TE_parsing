#!/usr/bin/python
'''
This script is meant to break up the sequence by indices into testing and trainings sets. Create bed files 
ssociated with the elements based on a grouping such as class, order, superfamily, and family.

Author: Maia Kaplan

USAGE: python breakup_intoBed.py TEclass data_breakup_15.txt
creates 3 files: data_breakup_15.txt1.bed data_breakup_15.txt2.bed    data_breakup_15.txtUnknown.bed
'''
import re
import fileinput
import pandas as pd 
from sys import argv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict



order_labels=["LTR","TIR","SINE","LINE","Helitron","Unknown"]
class_labels=["1","2","Unknown"]
superfamily_labels = ["LTR/Gypsy","LTR/Copia","DNA/HAT","DNA/Pogo","DNA/En-Spm","DNA/MuDR","DNA","DNA/Tc1","DNA/Harbinger","DNA/Mariner","LINE?","LINE/L1","SINE","RC/Helitron","RathE3_cons","RathE2_cons","RathE1_cons","Unassigned"]
family_labels =[] #initialize

#   orders: LTR, TIR, SINE, LINE, Helitron, Unknown
TEorder_dict={"LTR/Gypsy": "LTR","LTR/Copia" : "LTR","DNA/HAT" : "TIR","DNA/Pogo" : "TIR","DNA/En-Spm" : "TIR","DNA/MuDR" : "TIR","DNA" : "TIR","DNA/Tc1" : "TIR","DNA/Harbinger" : "TIR","DNA/Mariner" : "TIR", "LINE?" : "LINE","LINE/L1" : "LINE","SINE" : "SINE","RC/Helitron" : "Helitron","RathE3_cons" : "Unknown","RathE2_cons" : "Unknown","RathE1_cons" : "Unknown","Unassigned" : "Unknown"}
#   classes: 1) retrotransposon, and 2) DNA transposon, and Unknown
TEclass_dict = {"LTR" : "1","Helitron" : "2","TIR" : "2","LINE" : "1","SINE" : "1","Unknown" : "Unknown"}

TE_level = argv[1]
knownTE_file = argv[2]
train_i = argv[3] # the list for that particular replicate
test_i = argv[4]

d_knownTEs = pd.read_table(knownTE_file)
TE_ids = d_knownTEs["TE_ID(s)"].fillna("") #the id column is the name of the TE element,

train_ids = TE_ids[train_i]
test_ids = TE_ids[test_i]

 #TE_level can be TEclass, TEorder, TEfamily, TEsuperfamily
#I need this to output a list for each TE_level
groups = ''
if (TE_level == 'TEclass'):
    groups = class_labels
    print 'got the class!'
elif (TE_level == 'TEorder'):
    groups = order_labels
#elif(TE_level == TEfamily):
    #groups = 
elif(TE_level == 'TEsuperfamily'):
    groups = superfamily_labels
else:
    print 'Invalid TE level'

out_list=['nada']

infile = open(knownTE_file,'r')
header=infile.readline().rstrip().split("\t")  #   read in a single line from the file pointer

for line in infile:
    name,orientation5prime,start,end,family,superFamily=line.rstrip().split("\t")
    chromosome="chr"+re.search("^AT[0-9+]",name).group().replace("AT","")
    #not all the information is in the column of that file (the order and class)
    Seg_TEorder = TEorder_dict[superFamily]
    Seg_TEclass = TEclass_dict[Seg_TEorder]

    for each in groups: # so if we're looking at classes, this would be 1,2,unknown
        TEST_name = knownTE_file + each.replace("/","_") + '.test.bed' # make a file for each of the divisions
        TRAIN_name = knownTE_file + each.replace("/","_") + '.train.bed' # for both test and train

        with open(TEST_name, "a") as TEST, open (TRAIN_name, "a") as TRAIN: 
            ## header would go here if I could figure it out. 
            #if fileinput.isfirstline():
                #OUT.write('\t'.join(['#chromosome','start','stop','name']))   
            # for each line, divide it into the group necessary and write it to the corresponding file 
            if (TE_level == 'TEclass' and Seg_TEclass == each):
                if (name == train_ids):
                    TRAIN.write('\t'.join([chromosome,start,end,name])+'\n')
                elif(name -- test_ids):
                    TEST.write('\t'.join([chromosome,start,end,name])+'\n')
                else:
                    print "Not training or testing set"

            elif(TE_level == 'TEorder' and Seg_TEorder == each):
                if (name == train_ids):
                    TRAIN.write('\t'.join([chromosome,start,end,name])+'\n')
                elif(name -- test_ids):
                    TEST.write('\t'.join([chromosome,start,end,name])+'\n')
                else:
                    print "Not training or testing set"
            elif(TE_level == 'TEsuperfamily' and superFamily == each):
                                if (name == train_ids):
                    TRAIN.write('\t'.join([chromosome,start,end,name])+'\n')
                elif(name -- test_ids):
                    TEST.write('\t'.join([chromosome,start,end,name])+'\n')
                else:
                    print "Not training or testing set'
            else:
                print 'taxonomic division not found'


