#!/usr/bin/env python3

import argparse
from itertools import combinations
from itertools import product

def pairs_intra_file(input_file, pos, field_sep):
	list_of_elem = []
	with open(input_file, 'r', errors='replace') as infile:
		for line in infile:
			list_of_elem.append( ((line.split(field_sep)[pos]).rstrip()) )
	paired_list = list(combinations(list_of_elem, 2))
	for i in paired_list:
		print( (str(i[0]) + '\t' + str(i[1])) )



def pairs_across_files(input_file1, input_file2, pos1, pos2, field_sep1, field_sep2):
	list_of_elem1 = []
	list_of_elem2 = []
	with open(input_file1, 'r', errors='replace') as infile1, open(input_file2, 'r', errors='replace') as infile2:
		for line1 in infile1:
			list_of_elem1.append( ((line1.split(field_sep1)[pos1]).rstrip()) )
		for line2 in infile2:
                        list_of_elem2.append( ((line2.split(field_sep2)[pos2]).rstrip()) )
	paired_list = list(product(list_of_elem1, list_of_elem2))
	for i in paired_list:
                print( (str(i[0]) + '\t' + str(i[1])) )


def main(args):
	if args.in2:
		pairs_across_files(args.in1, args.in2, args.field1, args.field2, args.fs1, args.fs2)
	else:
		pairs_intra_file(args.in1, args.field1, args.fs1)


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Given a file makes all pairwise combination of item found in given column (specifable), MODE: one file.\nGiven two files it makes all possible pairwise combination between items of the two files, column specifiable, MODE: two file.\nA - B is equal to B - A, only one will be reported. If two files are given: file1 -> a,b,c ; file2 -> d,e; combinations reported will be: a-d, a-e, b-d, b-e, c-d, c-e.\nThe output is printed to standard output so redirect it to save it to a file.')
	parser.add_argument('--in1', type=str, required=True, help='only mandatory flag. input file n1, if only this is given the script will work in mode one file, all lines are expected to have the requested field/column.')
	parser.add_argument('--in2', type=str, required=False, help='input file n2, if this field is also given the script will work in mode two files, intra file pairs will not be reported. all lines are expected to have the requested field/column.')
	parser.add_argument('--field1', type=int, required=False, default=0, help='tells the script what field/column to use for the generation of the pairs, default 0 first field. This applies to file --in1. Hint pass -1 for the last field and so on.')
	parser.add_argument('--field2', type=int, required=False, default=0, help='tells the script what field/column to use for the generation of the pairs, default 0 first field. This value is only used in two file mode and refers to file --in2.')
	parser.add_argument('--fs1', type=str, required=False, default='\t', help='aka field separator, tells the script what digit or substring to use as sparator of columns/fields. Default <tab> = \\t. This applies to file --in1.')
	parser.add_argument('--fs2', type=str, required=False, default='\t', help=' tells the script what digit or substring to use as sparator of columns/fields. Default <tab> = \\t. This applies to file --in2. This value is only used in two file mode.')

	args = parser.parse_args()

	main(args)

