#!/usr/bin/env python3

import sys

def converter(infilename, outfilename):
    with open(infilename, 'r') as infile, open(outfilename, 'w') as outfile:
        infile.readline()
        pdb_id = ''
        resolution = ''
        for line in infile:
            if line[0] == '"':
                tempo = line.rstrip()
                splits = tempo.split('",')
                pdb_id = splits[0][1:]
                resolution = splits[1][1:]
                len_seq = len(splits[2][1:])
                if len_seq >= 50:
                    header = '>'  + pdb_id + '_' + splits[3][1] + '-' + resolution + '\n'
                    #print(header)
                    outfile.write(header)
                    len_seq = len(splits[2][1:])
                    #print(splits[2][1:len_seq])
                    seq = splits[2][1:len_seq] + '\n'
                    outfile.write(seq)
                else:
                    continue
            else:
                tmp = line.rstrip()
                splits1 = tmp.split(',"')
                len_seq1 = len(splits1[1]) -1
                if len_seq1 >= 50:
                    header1 = '>' + pdb_id + '_' + splits1[2][0] + '-' + resolution + '\n'
                    #print(header1)
                    outfile.write(header1)
                    len_seq1 = len(splits1[1]) -1
                    #print(splits1[1][:len_seq1])
                    seq1 = splits1[1][:len_seq1] + '\n'
                    outfile.write(seq1)
                else:
                    continue

if __name__ == "__main__":
    try:
        csv_file_path = sys.argv[1]
        output_filepath = sys.argv[2]
    except:
        print('Program usage: text.py <custom pdb tabular report with only pdb id resolution sequence and chain fields> <the filepath to the output multi-fasta file>', file=sys.stderr)
        raise SystemExit
    else:
        converter(csv_file_path, output_filepath)
