#!/usr/bin/env python3


import argparse


def aug_gff_to_fasta(in_aug_gff):
	output_filename = (in_aug_gff.split('/')[-1]).split('.')[0] + '_aug.fasta'
	with open(in_aug_gff, 'r', errors='replace') as ingff, open(output_filename, 'w') as outfile:
		transcript_count = 1							## to make unique headers for many transcripts
		buffer_sequence = ''							## to put all sequence on one line
		check = False
		for line in ingff:
			if 'name = ' in line and ')' in line:
				header = '>' + (line.split('name = ')[1]).split(')')[0]
			elif 'protein sequence = [' in line and ']' not in line:
				buffer_sequence += (line.split('protein sequence = [')[1]).rstrip()
				check = True
			elif 'protein sequence = [' in line and ']'  in line:                ## case where sequence is very short and on one line
				sequence = (line.split('protein sequence = [')[1]).split(']')[0] + '\n'
				outfile.write(sequence)
			elif check and ']' not in line:
				buffer_sequence += (line.split(' ')[-1]).rstrip()                 ## sequences are preceded by an # and a space
			elif check and ']' in line:
				buffer_sequence += (line.split(' ')[-1]).rstrip()[:-1] + '\n'     ## removing square brachets at end
				outfile.write( header +  '.' + str(transcript_count) + '\n' )
				outfile.write(buffer_sequence)
				check = False
				buffer_sequence = ''
				transcript_count += 1
			
def main(args):
	aug_gff_to_fasta(args.gff)



if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='takes as output the augustus gff standard format output, it then extraxts all the sequences it finds. It looks for two things: for the header the pattern   <name = >  and <)> as end   major and minor signs not included, after such pattern should be the name of the sequence, an additional  .1 or .2 ecc.. will be added in case or more transcripts;\n for the sequence it looks for the pattern     <protein sequence = [>  as start and   <]> as end.')
	parser.add_argument('--gff', type=str, required=True, help='input augustus gff')
	args = parser.parse_args()

	main(args)

