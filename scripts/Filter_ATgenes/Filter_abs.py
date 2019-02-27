
#!/usr/bin/python
'''
This script is meant to filter region annotated as introns in a focus 
species (FS) to extract those which are novel compared to related 
species (RS). 

The input is the results of liftover (orthologous regions of FS  mapped to RS), 
which are the mappedTo$RS.bed and unmappedTo$RS.bed files. 

The filter must match unmapped introns to cooresponding mapped exons.
Information must be kept about intron and exon coordinates of FS, and coordinates 
of where they map to in the RS. 

Usage: python Filter_param.py -f AL -r AT -d 80 -b 0.001 -m 0.01
Example: for i in AT CR TH SI EP; do python filter_params.py -f AL -r $i -d 80 -b 0.2 0.3 -m 0.4 0.5 > results_$i.bed; done

Fixed bug with reverse strand Nov 1st

'''
#imports:
from sys import argv
from itertools import product

import pandas as pd 
import argparse
import time

#input files:
annotation_name = 'AL_Weigel_gene_annotation'

def save_results(dir_name, mB, mM, lines, related_species):
	out_name = dir_name + '/NI_' + related_species + '.mB.' + str(mB) + '.mM' + str(mM)
	with open(out_name,'w') as f:
		f.write('\n'.join(lines))
	print 'saved to', out_name

def sp_liftovers(RS, result, mB, mM):
	file_name = annotation_name + '.' + result + 'To' + RS + '.mB' + str(mB) + '.mM' + str(mM) + '.bed'
	return file_name 

time_now = time.time() * 1000.0
def p_time():
	global time_now
	print 'time: %.2f ms' % (time.time() * 1000.0 - time_now)
	time_now = time.time() * 1000.0

def main(focus_species, related_species, max_dist, minBlock, minMatch):
	for mB, mM in product(minBlock, minMatch):
		unintron_df, mexon_df = combos_mapped(focus_species, related_species, mB, mM)
		lines = loop_through(unintron_df, mexon_df, max_dist)
		save_results('results', mB, mM, lines, related_species)


def combos_mapped(focus_species, related_species, mB, mM):
	header_starts = ['RS_chr','RS_start','RS_end', 'FS_info']
	header_info = ['AL_chr','AL_start','AL_end','AL_ann','dot','strand','dot2','geneinfo']

	mapped_filename = sp_liftovers(related_species, 'mapped', mB, mM)
	unmapped_filename = sp_liftovers(related_species, 'unmapped', mB, mM)

	mapped_file = [line.split() for line in open(mapped_filename) if line.strip()]
	unmapped_file = [line.split() for line in open(unmapped_filename) if line.strip()]

	mapped_df_extracol = pd.DataFrame(mapped_file)
	mapped_df = mapped_df_extracol.ix[:,0:3]
	mapped_df.columns = header_starts
	new_inf = pd.DataFrame(list(mapped_df.FS_info.str.split('#')))
	smaller_df = mapped_df.drop('FS_info', axis=1)
	new_inf.columns = header_info
	df_mapped = smaller_df.join(new_inf)
	mexon_df = df_mapped[df_mapped['AL_ann'].str.contains('exon')]


	unmapped_df_extrarows = pd.DataFrame(unmapped_file, columns = header_starts)
	drop_list = ['deleted','in'] #unmapped doens't have RS sequence because it doesn't map
	unmapped_df = unmapped_df_extrarows[~unmapped_df_extrarows['RS_start'].isin(drop_list)]	

	#Really bad, fix this:
	list_indexies = range(0,len(unmapped_df),1)
	unmapped_df.index = list_indexies # becasue they were offset once I removed some rows

	new_inf2 = pd.DataFrame(list(unmapped_df.FS_info.str.split('#')))
	smallerunmapped_df = unmapped_df.drop('FS_info',axis=1) 
	new_inf2.columns = header_info
	df_unampped = smallerunmapped_df.join(new_inf2)
	unintron_df = df_unampped[df_unampped['AL_ann'].str.contains('intron')]

	return (unintron_df, mexon_df)
	
def loop_through(unintron_df, mexon_df, max_dist):
	acceptors_list = []
	donors_list = []
	intron_da = {}
	# p_time()
	for index, row in unintron_df.iterrows():
		if (unintron_df.ix[index, 'strand'] == '+'):
			donor = unintron_df.ix[index, 'AL_chr'] + '_' + unintron_df.ix[index, 'AL_start']
			acceptor = unintron_df.ix[index, 'AL_chr'] + '_' + unintron_df.ix[index, 'AL_end']
		else:
			donor = unintron_df.ix[index, 'AL_chr'] + '_' + unintron_df.ix[index, 'AL_end']
			acceptor = unintron_df.ix[index, 'AL_chr'] + '_' + unintron_df.ix[index, 'AL_start']
		donors_list.append(donor)
		acceptors_list.append(acceptor)# I don't know how else to store information well. :( 
		intron_da[donor] = acceptor

	exon_5donor = {}
	exon_3acceptor = {}

	for index, row in mexon_df.iterrows():
		if (mexon_df.ix[index, 'strand'] == '+'):
			donor = mexon_df.ix[index, 'AL_chr'] + '_' + str(int(mexon_df.ix[index, 'AL_end'])+1)
			dmapped = mexon_df.ix[index, 'RS_chr'] + '_' + mexon_df.ix[index, 'RS_end']
			acceptor = mexon_df.ix[index, 'AL_chr'] + '_' + str(int(mexon_df.ix[index, 'AL_start'])-1)
			amapped = mexon_df.ix[index, 'RS_chr'] + '_' + mexon_df.ix[index, 'RS_start']
		else:
			donor = mexon_df.ix[index, 'AL_chr'] + '_' + str(int(mexon_df.ix[index, 'AL_start'])-1)
			dmapped = mexon_df.ix[index, 'RS_chr'] + '_' + mexon_df.ix[index, 'RS_start']
			acceptor = mexon_df.ix[index, 'AL_chr'] + '_' + str(int(mexon_df.ix[index, 'AL_end'])+1)
			amapped = mexon_df.ix[index, 'RS_chr'] + '_' + mexon_df.ix[index, 'RS_end']
		
		exon_5donor[donor] = dmapped, amapped #these are AT coordinates 
		exon_3acceptor[acceptor] = amapped, dmapped
	loss_donor = []
	loss_acceptor = []
	#full_loss = {}
	novel_intron = {}

	lines = []
	for idonor in intron_da.keys():  
		FS_donor = idonor # AL donor
		FS_acceptor = intron_da[idonor] #AL acceptor
		if (FS_donor in exon_5donor.keys() and FS_acceptor in exon_3acceptor.keys()): 
			AT_donorA = exon_5donor[FS_donor][0]
			AT_donorB = exon_5donor[FS_donor][1]
			AT_acceptorA = exon_3acceptor[FS_acceptor][0]
			AT_acceptorB = exon_3acceptor[FS_acceptor][1]
			# I could simplify this by taking out information and just using the coordinates 
			#	but I think there is value in keeping it for going back to check/debugging. Q
			AT_dcoordinateA = AT_donorA.split('_')[2]
			AT_dcoordinateB = AT_donorB.split('_')[2]
			AT_acoordinateA = AT_acceptorA.split('_')[2]
			AT_acoordinateB = AT_acceptorB.split('_')[2]
			
			AT_co = [int(AT_dcoordinateA), int(AT_dcoordinateB), int(AT_acoordinateA), int(AT_acoordinateB)]
			AT_co_sorted = sorted(AT_co)
			bps_btwn = abs(int(AT_co_sorted[1]) - int(AT_co_sorted[2]))
			#print bps_btwn
			#print AT_co_sorted

			if (bps_btwn < max_dist):
				#print max_dist
				#print bps_btwn
				novel_intron[idonor] = intron_da[idonor]
				scaff = FS_donor.split('_')[1]
				
				line = '\t'.join([scaff.strip(), FS_donor.split('_')[2].strip(), FS_acceptor.split('_')[2]])
				lines.append(line)
	return (lines)

#Arg parse: 
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Filter liftover results \
		to detect novel introns')
	parser.add_argument('-f', '--focus_species', type = str,
		help = 'select focus', dest = 'focus_species', required = True)
	parser.add_argument('-r', '--related_species', type = str, 
		help = 'select relative', dest = 'related_species', required = True)
	parser.add_argument('-d', '--max_dist', type = int, 
		help = 'max number of bps between mapped exons in RS, st that intron \
		was not previously present', dest = 'max_dist', required = True)
	parser.add_argument('-b', '--minBlock', type = float,
		help = 'minBlock values, can be list as 0.1 0.2', dest = 'minBlock', required = True, nargs = '+')
	parser.add_argument('-m', '--minMatch', type = float,
		help = 'minMatch values, can be list as 0.3 0.4', dest = 'minMatch', required = True, nargs = '+')


	args = parser.parse_args()
	main(args.focus_species, args.related_species, args.max_dist, args.minBlock, args.minMatch)

