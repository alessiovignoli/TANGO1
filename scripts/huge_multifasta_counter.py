#!/usr/bin/env python3

import sys

def line_counter(in_multifastapath):
    #print(in_multifastapath, n_of_outputfiles)
    n_lines = 0
    n_ids = 0
    with open(in_multifastapath, 'r') as infasta:
        for line in infasta:
            n_lines += 1
            if line[0] == '>':
                n_ids += 1
    print('number of lines = ', n_lines, '      number of ids = ', n_ids, '     file = ', in_multifastapath, end='')
    #print('number of ids = ', n_ids)



if __name__ == "__main__":
    try:
        in_multifastapath = sys.argv[1]
        #n_of_outputfiles = sys.argv[2]
    except:
        print('Programm usage: <Huge input file, could be enything, but if it is a multifasta (fasta in general) one also computes the number of ids>', file=sys.stderr)
        raise SystemExit
    else:
        line_counter(in_multifastapath)
