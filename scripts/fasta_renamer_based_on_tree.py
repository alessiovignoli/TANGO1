#!/usr/bin/env python3

# pp = posterior probability     always
# tm = trans membrane            always




import argparse



def  match_and_rewrite(in_tree, in_fasta, out_fasta, delimiter=""):
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
				



def main(args):
	
	if args.delimiter:
		match_and_rewrite(args.tree, args.fasta, args.output, args.delimiter)
	else:
		match_and_rewrite(args.tree, args.fasta, args.output)




if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='This script rename the headers of a fasta file and reorder the sequences based on the input tree file. The renaming convention is that the new id will be the first argument of the header split on spaces, with "tr|", "sp|" when present removed and then split again if possible on "|" taking first field. SO as an example from a header like this :\n>tr|A0A7N5JUE6|A0A7N5JUE6_AILME MIA SH3 domain ER export factor 3 OS=Ailuropoda melanoleuca OX=9646 GN=MIA3 PE=4 SV=1\nthe new header will be :\n>9646.A0A7N5JUE6\nwhere 9646 is the id foun in the newick tree associated with this sequnce.\nIn fact th emaqndatory thing for this script to work is that the keywords present in the newick are also present in the header of the sequences. If an id match is not found it will not be written to output.\n Lastly the scripts outputs a file in the directory where the script is lounched with the matching found between tree and fasta called -> correspondance_fasta_tree.tab ')
	parser.add_argument('-t', '--tree', type=str, required=True,  metavar="FILE",  help='path to the tree file in newick format')
	parser.add_argument('-f', '--fasta', type=str, required=True,  metavar="FILE",  help='path to fasta file with sequences on one line')
	parser.add_argument('-o', '--output', type=str, required=True,  metavar="FILE",  help='path to the output file')
	parser.add_argument('-d', '--delimiter', type=str, required=False, metavar="STR", help='optionalstring argument that is used as delimiter, use " or \' to encapsulate such string it can also be an empty string, the delimiter is directly concatenated with the id found in the tree and searched in the headers, two different concatenation are done  delimiter+id and delimiter+id+delimiter, both are valid searches')
	
	args = parser.parse_args()
	main(args)

