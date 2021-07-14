#!/usr/bin/env python3

# tm = trans membrane            always

import sys

def hydro_colist_prep(trimmflnm, renameflnm, outfilepath,  maskrange=None):
    #print(trimmflnm, renameflnm, outfile, maskrange)
    #print(outfilepath)
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
                                if aa == 'E' or aa == 'R':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '0' + '\n' 
                                    outfile.write(to_write_line)
                                elif aa == 'D' or aa == 'K':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '1' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'N' or aa == 'Q':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '2' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'H' or aa == 'P':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '3' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'S' or aa == 'T':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '4' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'Y' or aa == 'G':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '5' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'M' or aa == 'A':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '6' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'C' or aa == 'L':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '7' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'F' or aa == 'V':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '8' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'I' or aa == 'W':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '9' + '\n'
                                    outfile.write(to_write_line)
                        else:
                            l_ext = int(maskrange.split(',')[0]) - 1
                            r_ext = int(maskrange.split(',')[1]) 
                            #print(l_ext, r_ext)
                            for aa_pos in range(l_ext, r_ext):
                                aa = seq_line[aa_pos]
                                if aa == 'E' or aa == 'R':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '0' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'D' or aa == 'K':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '1' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'N' or aa == 'Q':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '2' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'H' or aa == 'P':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '3' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'S' or aa == 'T':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '4' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'Y' or aa == 'G':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '5' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'M' or aa == 'A':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '6' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'C' or aa == 'L':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '7' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'F' or aa == 'V':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '8' + '\n'
                                    outfile.write(to_write_line)
                                elif aa == 'I' or aa == 'W':
                                    to_write_line = new_header + ' ' + str((aa_pos+1)) + ' ' + '9' + '\n'
                                    outfile.write(to_write_line)

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
