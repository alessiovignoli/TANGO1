#!/usr/bin/env python3

import sys

def splitter(in_multifastapath, n_of_outputfiles, n_of_lines):
    #print(in_multifastapath, n_of_outputfiles)
    #seq_countinuity_counter = 0
    line_per_split = (n_of_lines // n_of_outputfiles) + 1
    split_on = line_per_split
    line_number = 0
    suffix = 1
    actual_filename = in_multifastapath
    if '/' in in_multifastapath:
        n_dirs_inpath = len(in_multifastapath.split('/')) - 1
        actual_filename = in_multifastapath.split('/')[n_dirs_inpath]
    n_suffix_inpath = len(actual_filename.split('.')) - 2
    prefix = actual_filename.split('.')[n_suffix_inpath]
    split_filename = prefix + '-splt' + str(suffix) + '.fasta'
    #print(line_per_split)
    with open(in_multifastapath, 'r') as infasta:
        for line in infasta:
            line_number += 1
            if line_number >= split_on:
                print('line number:', line_number, 'line_per_split :', line_per_split )
                if line[0] == '>':
                    #print('I m in the important if :', line)
                    split_on = line_per_split + line_number
                    suffix += 1
                    split_filename = prefix + '-splt' + str(suffix) + '.fasta'
                    split_file = open(split_filename, 'a+')
                    split_file.write(line)
                    split_file.close()
                else:
                    split_file = open(split_filename, 'a+')
                    stripped_line = line.rstrip()
                    split_file.write(stripped_line)
                    split_file.close()
            elif line[0] == '>':
                split_file = open(split_filename, 'a+')
                split_file.write('\n')
                split_file.write(line)
                split_file.close()
            else:
                split_file = open(split_filename, 'a+')
                stripped_line = line.rstrip()
                split_file.write(stripped_line)
                split_file.close()



if __name__ == "__main__":
    try:
        in_multifastapath = sys.argv[1]
        n_of_outputfiles = int(sys.argv[2])
        n_of_lines = int(sys.argv[3])
    except:
        print('Programm usage: <Huge input file, could be enything, one mustt lso computes the number of ids present in such file (using the script huge_multifasta_counter.py)> <number of times the file must be split, is the number of the output files that youwant to get at the end> < number of lines contained in the input file, as specified this has to be computed before the lounch of this script>', file=sys.stderr)
        raise SystemExit
    else:
        splitter(in_multifastapath, n_of_outputfiles, n_of_lines)

