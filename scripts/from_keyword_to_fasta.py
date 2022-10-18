#!/usr/bin/env python3

import sys
import gzip

def fasta_searcher(keywords_file, file_to_search, output_filepath):
    list_of_headers = []
    #print(output_filepath)
    with open(output_filepath, 'w') as outfile, open(keywords_file, 'r') as infile: 
        #print(infile, outfile)
        for in_line in infile:
            #print(in_line, end='')
            search_key = in_line.rstrip()
            list_of_headers.append(search_key)
        searchfile = None
        if file_to_search.endswith('.gz'):
            searchfile = gzip.open(file_to_search, 'r')
        else:
            searchfile = open(file_to_search, 'r')
        #print(searchfile)
        linebuffer = ''
        hit_counter = 0
        for searchline in searchfile:
            #print(type(searchline))
            decoded_line = searchline
            if isinstance(searchline, bytes):
                decoded_line = searchline.decode()
            if decoded_line[0] == '>':
                if hit_counter == 0:
                    for keyword in list_of_headers:
                        if keyword in decoded_line:
                            outfile.write(decoded_line)
                            hit_counter = 1
                            break
                else:
                    hit_counter = 0
                    seq_line = linebuffer + '\n'
                    outfile.write(seq_line)
                    linebuffer = ''
                    for keyword in list_of_headers:
                        if keyword in decoded_line:
                            outfile.write(decoded_line)
                            hit_counter = 1
                            break
            elif hit_counter == 1:
                stripped_seq_line = decoded_line.rstrip()
                linebuffer += stripped_seq_line
            else:
                continue
        if linebuffer != '':
            towrite = linebuffer + '\n'
            outfile.write(towrite)
        searchfile.close()


if __name__ == "__main__":
    try:
        ids_file_path = sys.argv[1]
        filepath_to_search = sys.argv[2]
        output_file = sys.argv[3]
    except:
        print('Program usage: text.py <one id per line file> <a fasta file path to the file that is supposed to contain the input ids, it can be gzipped> <The path to the file on which the outputfile has to be written>', file=sys.stderr)
        raise SystemExit
    else:
        fasta_searcher(ids_file_path, filepath_to_search, output_file)

