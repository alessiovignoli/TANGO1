#!/usr/bin/env python3

# tm = trans membrane            always

import sys

def hydro_colist_prep(trimmflnm, renameflnm, outfilepath,  scale, maskrange=None):
    #print(trimmflnm, renameflnm, outfile, maskrange)
    #print(outfilepath)
    majority_hydro_dict = {'A':6, 'C':7, 'D':1, 'E':0, 'F':8, 'G':5, 'H':3, 'I':9, 'K':1, 'L':7, 'M':6, 'N':2, 'P':3, 'Q':2, 'R':0, 'S':4, 'T':4, 'V':8, 'W':9, 'Y':5, 'X':'', 'B':'10', 'U':'10', 'Z':'10'}
    kyte_doolittle_hydro_dict = {'A':6, 'C':7, 'D':2, 'E':2, 'F':8, 'G':6, 'H':3, 'I':9, 'K':0, 'L':8, 'M':7, 'N':1, 'P':3, 'Q':1, 'R':0, 'S':5, 'T':5, 'V':9, 'W':4, 'Y':4, 'X':'10', 'B':'10', 'U':'10', 'Z':'10'}
    GES_hydro_dict = {'A':6, 'C':7, 'D':0, 'E':1, 'F':9, 'G':5, 'H':3, 'I':8, 'K':1, 'L':8, 'M':9, 'N':2, 'P':4, 'Q':2, 'R':0, 'S':4, 'T':5, 'V':7, 'W':6, 'Y':3, 'X':'10', 'B':'10', 'U':'10', 'Z':'10'}
    UHS_hydro_dict = {'A':6, 'C':4, 'D':0, 'E':1, 'F':9, 'G':7, 'H':3, 'I':9, 'K':0, 'L':8, 'M':7, 'N':2, 'P':2, 'Q':3, 'R':1, 'S':4, 'T':5, 'V':8, 'W':5, 'Y':6, 'X':'10', 'B':'10', 'U':'10', 'Z':'10'}
    to_be_used_dict = None
    if scale == "kyte":
        to_be_used_dict = kyte_doolittle_hydro_dict
    elif scale == "GES":
        to_be_used_dict = GES_hydro_dict
    elif scale == 'UHS':
        to_be_used_dict = UHS_hydro_dict
    elif scale == 'majority':
        to_be_used_dict = majority_hydro_dict
    else:
        print('please specify a valid argument for the selection of the Hydrophobicity scale, valid are "kyte", "GES", "UHS" and "majority" it has been given:', scale, file=sys.stderr)
        raise SystemExit
    with open(renameflnm, 'r') as rename_file, open(outfilepath, 'w') as outfile:
        #print(fasta_file, rename_file)
        for line in rename_file:
            old_header = (line.split()[0]).rstrip()
            new_header = (line.split()[1]).rstrip()
            i = 0
            #print(old_header, new_header)
            #print(i)
            with open(trimmflnm, 'r') as fasta_file:
                for fline in fasta_file:
                    #print(fline, end='')
                    if old_header in fline:
                        i += 1
                        #print(fline, old_header, end='')
                    elif i == 1:
                        seq_line = fline.rstrip()
                        if maskrange == None:
                            l_ext = 0
                            r_ext = len(seq_line)
                            #print('>>', line + fline, end='')
                            #print(old_header, '  ', r_ext, end='')
                            #print(old_header + '\n' + seq_line[:r_ext])
                            for aa_pos in range(l_ext, r_ext):
                                aa = seq_line[aa_pos]
                                to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + str(to_be_used_dict[aa]) + '\n'
                                outfile.write(to_write_line)
                        else:
                            l_ext = int(maskrange.split(',')[0]) - 1
                            r_ext = int(maskrange.split(',')[1]) 
                            #print(l_ext, r_ext)
                            for aa_pos in range(l_ext, r_ext):
                                aa = seq_line[aa_pos]
                                to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + str(to_be_used_dict[aa]) + '\n'
                        i = 0



if __name__ == "__main__":
    mask_range = None
    out_filepath = None
    try:
        trimmed_fasta_filepath = sys.argv[1]
        rename_filename = sys.argv[2]
        out_filepath = sys.argv[3]
        mask_range = sys.argv[4]
    except:
        if mask_range is None and out_filepath is not None:
            hydro_colist_prep(trimmed_fasta_filepath, rename_filename, out_filepath)
        else:
            print('Program usage: text.py <a trimmed/normal multifasta file> <two column file, space separeted, having the final ids for each header of the multi-fasta file in the second row> <optional field that tells the script from which position on the sequence has to start the colouring and where to finish, specified extemity are coloured>', file=sys.stderr)
            raise SystemExit
    else:
        hydro_colist_prep(trimmed_fasta_filepath, rename_filename, out_filepath, mask_range)
