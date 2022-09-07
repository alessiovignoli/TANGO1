#!/usr/bin/env python3


import argparse


def xml_fasta_transformer(in_xml, out_file):
	with open(in_xml, 'r') as inxml, open(out_file, 'w') as outfile:
		for line in inxml:
			if 'id="' in line:
				fasta_header = line.split('id="')[1].split('"')[0]
				outfile.write( '>' + fasta_header + '\n' )
				next_line = inxml.readline()
				if '>' in next_line and '<' in next_line:
					sequence = next_line.split('>')[1].split('<')[0]
					if sequence.isupper() and sequence.isalpha():
						outfile.write( sequence+'\n' )
	return out_file



def main(args):
	outfile_path = xml_fasta_transformer(args.xml, args.out)
	return outfile_path


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='transforming a ensembl xml file into a fasta one')
	parser.add_argument('--xml', type=str, required=True, help='input xml file with   id="someID"  and  the following line containing the sequence between  >   and   <  signs')
	parser.add_argument('--out', type=str, required=True, help='output filename or path')
	args = parser.parse_args()

	main(args)

