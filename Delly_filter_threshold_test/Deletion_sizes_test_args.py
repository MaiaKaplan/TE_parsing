#!/usr/bin/python
'''
This script is meant to parse a bed file from a delly vcf <DEL> file that match the desired criteria. 

default thresholds are as follows:
-q > 20 : the minimin MAPQ accepted
-s > 60 : the minimum number of split reads that should support the variant
-r > 300: the number of read counts that support the variant
-l > 100: the variant, here only deletions, must be over 100 basepairs long

PASS and PRECISE are always required. 


USAGE: python Deletion_sizes_test_args.py -i Region_chr1_dellymini_DEL.vcf > Result.bed

'''
import os
import re
from sys import argv
import pandas as pd 
import argparse
import linecache
import operator
import argparse

def main(INFILE, alt, min_qual, min_sr, min_rc, length):

	data = [line.split() for line in open(INFILE) if  line.strip()]
	data = data[735:]
	data = pd.DataFrame(data[1:], columns=data[0])

	filter_w = 'PASS' # always
	precision_w = 'PRECISE' # always

	# select just the rows of the dataframe that meet the variant type (ALT) as well as the PASS 
	passing_df = data.loc[(data['ALT'] == alt) & (data['FILTER'] == filter_w)]

	for row in passing_df.itertuples():
		#print row
		INFO_col = row.INFO	
		FORMAT_col = row._10 #I know, I know
		# values of the info column that I am interested in:
		prec = INFO_col.split(';')[0]
		MAPQ_info = INFO_col.split(';')[7]
		MAPQ = MAPQ_info.split('=')[1]
		
		chrom = row[1] 
		start_value = row[2]
		end_info = INFO_col.split(';')[4]
		end_value = end_info.split('=')[1]

		#if the row has the right split read support (precise rathe than not INPRECISE)
		if  (prec== "PRECISE" and int(MAPQ) > int(min_qual)):
			sr_info = INFO_col.split(';')[11]
			sr_value = sr_info.split('=')[1]
			rc = FORMAT_col.split(':')[5]

			# next, check if the split reads and read counts are above the thresholds
			if (int(sr_value) > int(min_sr) and int(rc) >int(min_rc)) :

				# then check that the length of the deletion is long enough (not compatible to)
				if (int(end_value)-int(start_value) > length):
					print chrom, start_value, end_value, MAPQ, sr_value, rc

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='My Fancy Program')

    parser.add_argument('-i', '--input', type=str,
                        help="This is the input file!",
                        dest='input', required=True)
    # parser.add_argument('-o', '--output', type=str,
    #                     help='Output file!!!',
    #                     dest='output', required=True)
    parser.add_argument('-a', '--alt', type=str,
                        help='Type of variant',
                        dest='alt', required=False, default="<DEL>")
    parser.add_argument('-q', '--qual', type=int,
                        help='min MAPQ accepted',
                        dest='qual', required=False, default="20")
    parser.add_argument('-s', '--split_read_support', type=int,
                        help='min SR accepted',
                        dest='sr', required=False, default="60")
    parser.add_argument('-r', '--read_count_sv', type=int,
                        help='min RC accepted',
                        dest='rc', required=False, default="300")
    parser.add_argument('-l', '--length', type=int,
                        help='length of ',
                        dest='len', required=False, default="100")

    args = parser.parse_args()
    main(args.input, args.alt, args.qual, args.sr, args.rc, args.len)
