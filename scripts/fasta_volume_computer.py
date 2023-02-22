#!/usr/bin/env python3

import argparse



def volume_computer(infasta, outname):
	
	# Volume scale used, from   StephenJ. PERKINS 
	# Protein volumes and hydration effects
	# The calculations of partial specific volumes, neutron scattering matchpoints
	# and 280 nm absorption coefficientsfor proteins and glycoproteins from amino acid sequences
	# Eur. J. Biochem. 157, 169-180 (1986)

	# the scale is n x 10totheminus3 nm3
	
	dict_volumes = {'A':[87.8, 2.3], 'C':[105.4, 5.0], 'D':[115.4, 2.2], 'E':[140.9, 5.3], 'F':[189.7, 7.4], 'G':[59.9, 2.2], 'H':[156.3, 6.1], 'I':[166.1, 3.4], 'K':[172.7, 5.9], 'L':[168.0, 4.3], 'M':[165.2, 1.8], 'N':[120.1, 4.1], 'P':[123.3, 1.8], 'Q':[145.1, 5.1], 'R':[188.2, 9.6], 'S':[91.7, 1.8], 'T':[118.3, 2.3], 'W':[227.9, 3.8], 'Y':[191.2, 8.0], 'V':[138.8, 3.6] }

	"""
	print(len(dict_volumes))
	for i in dict_volumes:
		print(i, dict_volumes[i])
	"""
	
	with open(infasta, 'r') as infa, open(outname, 'w') as out_file:
		tot_aa = 0							## used as counter
		tot_volume = 0.0						## used as buffer
		tot_err = 0.0							## comupting mean error
		tmp_id = infa.readline().rstrip().split(' ')[0][1:] 			## just for warning message and to not divide by zero on first header

		# Header of outfile
		out_file.write('seq ID\tavg Volume\tstd err\tlen seq\n')

		for line in infa:
			if line[0] == '>':
				# Compute averages
				avg_volume = round((tot_volume / tot_aa), 1)
				avg_err = round((tot_err / tot_aa), 1)
	
				# write to outfile
				out_file.write( (tmp_id + '\t' + str(avg_volume) + '\t' + str(avg_err) + '\t' + str(tot_aa) + '\n') )
	
				# Reset buffer variables
				tot_aa = 0
				tot_volume = 0.0
				tot_err = 0.0
				tmp_id = line.rstrip().split(' ')[0][1:]	
				
			else:
				for letter in line.rstrip():
					if letter.upper() in dict_volumes:
						tot_aa += 1
						tot_volume += dict_volumes[letter.upper()][0]
						tot_err += dict_volumes[letter.upper()][1]

					# if the digit is not in the canonical alphabet send warning and skip from computation
					else:
						print('WARNING digit not recognized: ', letter, '  in fasta : ', tmp_id)
		

		# Write the info of the last sequence
		avg_volume = round((tot_volume / tot_aa), 1)
		avg_err = round((tot_err / tot_aa), 1)
		out_file.write( (tmp_id + '\t' + str(avg_volume) + '\t' + str(avg_err) + '\t' + str(tot_aa) + '\n') )
	
	return out_file



def main(args):
	#print(args.infa, args.out)
	
	out_filename = volume_computer(args.infa, args.out)
	


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='This script computes the average volume of a peptide. Knowing the volume each residue occupies it is just summing and averaging for the total number of residues in the sequence. A fasta file with many sequence can be provided and the average volume is computed for each of them.   \n###  WARNING  ####\n non canonical amminoacid and any other letter not in the aa alphabet will not be considered for the computation of the final Volume.')
	parser.add_argument('-f', '--infa', type=str, required=True, metavar="FILE", help='path to the input fasta file')
	parser.add_argument('-o', '--out', type=str, required=True, metavar="FILE", help='output filename/path, the file where the info will be written')
	
	args = parser.parse_args()
	main(args)
