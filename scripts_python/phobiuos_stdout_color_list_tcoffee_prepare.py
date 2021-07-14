#!/usr/bin/env python3

# pp = posterior probability     always
# tm = trans membrane            always

import sys
import os

def phobstout_colist_prep(phobstdout, output_f, renameflnm, plp_dir=False, trimm_val=False, signalpept_val=False, switch_col=False):
    #print(phobstdout, renameflnm, plp_dir, trimm_val)
    cut_site = 0
    if plp_dir != False:
                list_of_plps = os.listdir(plp_dir)
                with open(phobstdout, 'r') as phob_out, open(output_f, 'w') as out_file:
                    counter = 0    # used as a check for the presence of a plp file with keyword inside it
                    plp_keyword = ''
                    for phob_line in phob_out:
                        if 'ID' in phob_line:
                            cut_site = 0
                            plp_keyword = phob_line.split(' ')[3].rstrip()
                            #print(plp_keyword)
                            counter = 0
                            for plp_filename in list_of_plps:
                                if plp_keyword in plp_filename and plp_filename.endswith(".plp"):
                                    counter = 1
                                    #print(plp_filename)
                                    plp_filepath = plp_dir + plp_filename
                                    with open(plp_filepath, 'r') as plp_file:
                                        for res_line in plp_file:
                                            if res_line[0] == '#':
                                                continue
                                            else:
                                                posterior_prob = float(res_line.split()[4])
                                                res_index = int(res_line.split()[0])
                                                #print(posterior_prob)
                                                if posterior_prob >= 0.90 and res_index >= signalpept_val:   # the threshold for signal peptides
                                                    left_cut = int(trimm_val[0])
                                                    cut_site = max((res_index - left_cut), 0)
                                    break
                            if counter == 0:
                                print('This sequence has not been found in the specified plp directory: ', plp_keyword, file=sys.stderr)
                        elif 'TRANSMEM' in phob_line and counter == 1:
                            #print(cut_site)
                            #print(phob_line, end='')
                            pred_tm_start = int(phob_line[14:21].strip())
                            pred_tm_end = int(phob_line[21:28].strip()) + 1
                            with open(renameflnm, 'r') as rename_file:
                                for renem_line in rename_file:
                                    if plp_keyword in renem_line:
                                        #print(renem_line, end='')
                                        new_name = renem_line.split(' ')[1].rstrip()
                                        for i in range (pred_tm_start, pred_tm_end):
                                            j = i - cut_site
                                            length_trimm = None
                                            try:
                                                length_trimm = int(trimm_val[0]) + int(trimm_val[1])
                                            except:
                                                length_trimm = 50000
                                            if j > 0 and j < length_trimm:
                                                if 'SPECIAL' in phob_line and switch_col != False:
                                                    #print(new_name, j, 0)
                                                    line_to_write = new_name + ' ' + str(j) + ' 0\n'
                                                    out_file.write(line_to_write)
                                                else:
                                                    #print(new_name, j, 1)
                                                    line_to_write = new_name + ' ' + str(j) + ' 1\n'
                                                    out_file.write(line_to_write)
                                        break


    else:
        with open(phobstdout, 'r') as phob_out, open(output_f, 'w') as out_file:
            plp_keyword = ''
            for phob_line in phob_out:
                if 'ID' in phob_line:
                    plp_keyword = phob_line.split(' ')[3].rstrip()
                elif 'TRANSMEM' in phob_line:
                    tm_start = int(phob_line[14:21].strip())
                    tm_end = int(phob_line[21:28].strip()) + 1
                    with open(renameflnm, 'r') as rename_file:
                        for renem_line in rename_file:
                            if plp_keyword in renem_line:
                                new_name = renem_line.split(' ')[1].rstrip()
                                for m in range (tm_start, tm_end):
                                    if 'SPECIAL' in phob_line and switch_col != False:
                                        #print(new_name, m, 0)
                                        line_to_write = new_name + ' ' + str(m) + ' 0\n'
                                        out_file.write(line_to_write)
                                    else:
                                        #print(new_name, m, 1)
                                        line_to_write = new_name + ' ' + str(m) + ' 1\n'
                                        out_file.write(line_to_write)
                                break



if __name__ == "__main__":
    rename_filename = None
    input_plp_dir = None
    trimm_value = None
    signalpept_value = None
    special_tm = None
    try:
        phob_stout_filename = sys.argv[1]
        output_file = sys.argv[2]
        rename_filename = sys.argv[3]
        input_plp_dir = sys.argv[4]
        trimm_value = sys.argv[5].split(',')
        signalpept_value = int(sys.argv[6])
        special_tm = sys.argv[7]
    except:
        if trimm_value is not None and input_plp_dir is not None and rename_filename is not None and signalpept_value is not None:
            phobstout_colist_prep(phob_stout_filename, output_file, rename_filename, input_plp_dir, trimm_value, signalpept_value)
        elif trimm_value is None and input_plp_dir is None and rename_filename is not None:
            phobstout_colist_prep(phob_stout_filename, output_file, rename_filename)
        else:
            print('Program usage: text.py <a phobius default (-long flag) standard output redirected into a file, if the trimm value option is given, the ID lines must contain a name that is a keyword to identify a .plp file of that specific sequence> <two column file, space separeted, having the final ids in the second row. the sequence name found after ID has to be a substring of the old name found in the rename file and the new name the actual name we want t-coffee to recognize> <if the trimm value option is given, The directory in which search all .plp files of the sequences present in the phobius standard output input file for this script> <trimm_value used to cut the sequences, this enebles all previous options, the sequences have been cutted from another script trimm_multifasta.py, go check it to see how the cut is performed, trimm_value has to be the same number used for that script> < A signal peptide value used when trimming good pp values will be ignored if their position on the sequence is lower than the value specified it must be an integer> <a boolean parameter used to switch to a different (two colouring) coloration mostly used when there aree two different types of transmembranes>', file=sys.stderr)
            raise SystemExit
    else:
        phobstout_colist_prep(phob_stout_filename, output_file, rename_filename, input_plp_dir, trimm_value, signalpept_value, special_tm)
