#!/usr/bin/env python3


import argparse


def fasta_longest(in_fasta):
	with open(in_fasta, 'r') as infasta:
		header = ''
		sequence = ''
		buffer_seq = ''							## for the case where the sequence is not on one line
		buffer_header = ''
		for line in infasta:
			if line[0] == '>':
				if not sequence:			
					## the first iteration will be empty, but in the second it will innclude the first sequence and header
					sequence = buffer_seq
					header = buffer_header
				elif len( buffer_seq ) < len( sequence ):		## using the sequence saved itself as variable for the max length found
					sequence = buffer_seq
					header = buffer_header
				buffer_seq = ''
				buffer_header = line.rstrip()
			else:
				buffer_seq += line.rstrip()			## for the sequence on more than one line
		if len( buffer_seq ) <  len( sequence ):			## last line of fasta file has not >
			print( buffer_header )
			print( buffer_seq )
		else:
			print(header)
			print(sequence)



def main(args):
	fasta_longest(args.fasta)



if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='picks the shortest sequence and corresponding header present in a fasta file and prints it to standardout')
	parser.add_argument('--fasta', type=str, required=True, help='input fasta file path')
	args = parser.parse_args()

	main(args)

