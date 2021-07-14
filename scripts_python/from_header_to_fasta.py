#!/usr/bin/env python3

import sys

def fasta_retriever(headers_file, list_files_to_search, output_filepath):
    list_of_headers = []
    print(output_filepath)
    with open(output_filepath, 'w') as outfile, open(headers_file, 'r') as infile:
        #print(infile, outfile)
        len_list = len(list_files_to_search)
        list_lines_search = []
        for i in range(0,len_list):
            #print(i)
            file_to_open = list_files_to_search[i]
            #print(file_to_open)
            with open(file_to_open, 'r') as searchfile:
                #print(searchfile)
                linebuffer = ''
                for line in searchfile:
                    if line[0] == '>':
                        linebuffer = line
                        #print(linebuffer, end='')
                    else:
                        linebuffer += line
                        list_lines_search.append(linebuffer)
        for in_line in infile:
            #print(in_line, end='')
            search_key = in_line.rstrip()
            for elem in list_lines_search:
                if search_key in elem:
                    new_header = elem.replace(' ', '_', 1) 
                    outfile.write(new_header)
                    print(elem, end='')
                    break



if __name__ == "__main__":
    try:
        headers_file_path = sys.argv[1]
        one_strins_filenames = sys.argv[2]
        output_path = sys.argv[3]
        list_of_filepaths = (one_strins_filenames.rstrip()).split(',')
    except:
        print('Program usage: text.py <one header per line file> <a bunch of fasta file paths (better if absolutes) separated by a comma like /home/file1.fasta,/home/file2.fasta and so on> <The path (absolute) to the file on which the outputfile has to be written with last slash>', file=sys.stderr)
        raise SystemExit
    else:
        #print(headers_file_path)
        #print(list_of_filepaths)
        #print(output_path)
        #print(list_of_filepaths, end= '')
        fasta_retriever(headers_file_path, list_of_filepaths, output_path)

