#!/usr/bin/env python3


import sys

def colist_prep(inrf, ff):
    #print('bubba')
    with open(inrf, 'r') as ranges_file, open(ff, 'r') as fasta_file:
        #print(ranges_file, fasta_file)
        list_fasta = fasta_file.readlines()
        for line in ranges_file:
            uniprot_id = line.split(' ')[0]
            start = int(line.split(' ')[1])
            end = int((line.split(' ')[2]).rstrip())
            header = ''
            for list_elem in list_fasta:
                if uniprot_id in list_elem:
                    header = (list_elem.rstrip())[1:]
            #print(uniprot_id, start, end, header)
            for i in range((start-5),(end+6)):
                #print(i)
                if i < start:
                    print(header, i, '2')
                elif i > end:
                    print(header, i, '2')
                else:
                    print(header, i, '0')

if __name__ == "__main__":
    try:
        input_filename = sys.argv[1]
        fasta_filename = sys.argv[2]
    except:
        print('Program usage: text.py <uniprot_id TM_range .txt file> <fasta_file.txt>', file=sys.stderr)
        raise SystemExit
    else:
        colist_prep(input_filename, fasta_filename)
