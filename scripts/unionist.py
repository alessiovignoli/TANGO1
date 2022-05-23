#!/usr/bin/env python3

# this script takes in input two arguments: 
# a bunch of file path (better if absolutes) separated by a comma like /home/file1.fasta,/home/file2.fasta and so on
# The path (absolute) to the directory on which the outputfile has to be written (the script automatically generates the name
# does the union based on the one-line sequence associated to each header
# it can do the union for n files and it also filters out the sequences that do present one or more <X> charachters
# the output name is generated automatically, simply by putting <u> <the number of sequences (n)> 
# and all the filenames used for the union without extention separated by <->

import sys

def unionist(list_of_files):
    #print(output_dir)
    output_filename = 'u' + str(len(list_of_files))
    dict_of_all_files ={}
    for filepath in list_of_files:
        path_splitted = filepath.split('/')
        filename = (path_splitted[(len(path_splitted)-1)]).split('.')[0]
        #print(filename)
        output_filename += '-' + filename
        with open(filepath, 'r') as f1:
            #print(f1)
            heder_buffer = ''
            for line in f1:
                if line[0] == '>':
                    #print(line)
                    heder_buffer = line.rstrip()
                else:
                    if 'X' in line:
                        continue
                    else:
                        stripped_line = line.rstrip()
                        dict_of_all_files[stripped_line] = heder_buffer
    #print(len(dict_of_all_files))
    output_file = output_filename + '.fasta'
    #print(output_file)
    with open(output_file, 'w') as out_file:
        for seq in dict_of_all_files:
            if ' | ' in dict_of_all_files[seq]:
                elems_of_header = dict_of_all_files[seq].split(' | ')
                if elems_of_header[2] == '1:1 ortholog':
                    new_header = elems_of_header[0]+'-'+elems_of_header[1]+'-'+elems_of_header[3]+'\n'
                    out_file.write(new_header)
                    out_file.write(seq)
                    out_file.write('\n')
                    #print(new_header)
                    #print(seq)
                else:
                    new_header = elems_of_header[0]+'-'+elems_of_header[1]+'-'+elems_of_header[2]+'\n'
                    out_file.write(new_header)
                    out_file.write(seq)
                    out_file.write('\n')
                    #print(new_header)
                    #print(seq)
            else:
                print(dict_of_all_files[seq])
                print(seq)

if __name__ == "__main__":
    try:
        one_strins_filenames = sys.argv[1]
        #output_dir = sys.argv[2]
        list_of_filepaths = (one_strins_filenames.rstrip()).split(',')
    except:
        print('Program usage: text.py  <a bunch of file path (better if absolutes) separated by a comma like /home/file1.fasta,/home/file2.fasta and so on> <The path (absolute) to the directory on which the outputfile has to be written with last slash (yes home/dir/) (the script automatically generates the name)>', file=sys.stderr)
        raise SystemExit
    else:
        #print(list_of_filepaths, end= '')
        unionist(list_of_filepaths)#, output_dir)
