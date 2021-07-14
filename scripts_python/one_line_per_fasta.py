#!/usr/bin/env python3


import sys

def one_line_seq(infile, outfile):
    with open(infile, 'r') as infasta, open(outfile, 'w') as outfasta:
        #print(infasta, infasta)
        first_iteration_counter = 0
        for line in infasta:
            if line[0] == '>' and first_iteration_counter == 0:
                outfasta.write(line)
                first_iteration_counter = 1
            elif line[0] == '>' and first_iteration_counter == 1:
                header = "\n" + line
                outfasta.write(header)
            else:
                stripped_line = line.rstrip()
                outfasta.write(stripped_line)
        outfasta.write("\n")




if __name__ == "__main__":
    try:
        input_fasta = sys.argv[1]
        output_fasta = sys.argv[2]
    except:
        print('Program usage: text.py <input multifatsa or single fasta file that have the sequencespread on more than one line> <output filename of the one line sequenceper fasta file>', file=sys.stderr)
        raise SystemExit
    else:
        one_line_seq(input_fasta, output_fasta)
