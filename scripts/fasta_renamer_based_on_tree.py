#!/usr/bin/env python3

# pp = posterior probability     always
# tm = trans membrane            always




import sys



def  match_and_rewrite(in_tree, in_fasta, out_fasta, delimiter):
	taable_of_correspondancefile = "correspondance_fasta_tree.tab"
	with open(in_tree,'r') as intree:
		newick_ids = (intree.readline().replace("(", "").replace(")", "")[:-2]).split(",")
		for elem in newick_ids:
			searchkeys = [(delimiter + elem), (delimiter + elem + delimiter)]
			with open(in_fasta, 'r') as infasta, open(out_fasta, 'a') as outfasta, open(taable_of_correspondancefile, 'a') as outtab:
				for fastaline in infasta:
					if '>' in fastaline and (searchkeys[0] in fastaline or searchkeys[1] in fastaline):
						fastaid = ((fastaline.split(' ')[0]).replace('tr|', '').replace('sp|', '')).split('|')[0][1:]
						outfasta.write(( '>' + elem + '.' + fastaid + '\n'))
						sequence = infasta.readline()
						outfasta.write( sequence )
						outtab.write(( elem + ' ' + fastaid + '\n' ))
						break
				




if __name__ == "__main__":
	try:
		input_treefile = sys.argv[1]
		input_fastafile = sys.argv[2]
		output_fastafile = sys.argv[3]
		input_delimiter = str(sys.argv[4])
	except Exception:
		print('Program usage: text.py <a newick format tree, with the id that have to have a one to one correspondance with a field in the header of a fasta sequence in the fasta file> <a fasta file \n####  WARNING  ####\n sequences on one line,\n the fasta file that as said has to have in each header a correspondance with the elements found in the newick tree> <the output file path of the renamed fasta> <an optionalstring argument that is used as delimiter, use " or \' to encapsulate such string it can also be an empty string, the delimiter is directly concatenated with the id found in the tree and searched in the headers, two different concatenation are done  delimiter+id and delimiter+id+delimiter, both are valid searches>\n\n the script will also produce in the dir in which the script is lounched a fiell called     correspondance_fasta_tree.tab    that stores on each line the tree id space fasta sequence id, that is "found" by the script', file=sys.stderr)
		raise SystemExit
	else:
		match_and_rewrite(input_treefile, input_fastafile, output_fastafile, input_delimiter)
