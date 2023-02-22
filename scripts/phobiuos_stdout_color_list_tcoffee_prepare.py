#!/usr/bin/env python3


##
#
#
#                       This script is a mess and should be re-written from scratch some day
#
##






# pp = posterior probability     always
# tm = trans membrane            always

import argparse
import sys
import os
sys.path.append(os.path.abspath("/home/alessio/Desktop/erasmus-internship/scripts"))
from phobius_short_prediction_field_retriever import phobius_short_pred_field_selecter


def phobstout_colist_prep(phobstdout, renameflnm, output_f, trimm_val, plp_dir, signalpept_val, switch_col):
    #print(phobstdout, output_f, renameflnm, plp_dir, trimm_val)
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
                                    plp_filepath = os.path.join(plp_dir, plp_filename)
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




def phobstout_short_colist_prep(phobstdout, renameflnm, output_f, trimm_val, plp_dir, signalpept_val, switch_col):
    #print(phobstdout, output_f, renameflnm, plp_dir, trimm_val, signalpept_val, switch_col)
    if trimm_val:
        list_of_plps = os.listdir(plp_dir)
        with open(phobstdout, 'r') as phob_out, open(output_f, 'w') as out_file:
            #print('this is a mess')
            
            for line in phob_out:
                seq_id = None
                list_boundaries_normal = None
                list_boundaries_special = []                                                            #to not throw error afterwards
                if switch_col:
                    seq_id, list_boundaries_normal = phobius_short_pred_field_selecter(line, 'n')
                    tmp, list_boundaries_special = phobius_short_pred_field_selecter(line, 's')
                else:
                    seq_id, list_boundaries_normal = phobius_short_pred_field_selecter(line, '-')
                cut_site = 0
                counter = 0    # used as a check for the presence of a plp file with keyword inside it
                for plp_filename in list_of_plps:
                    if seq_id in plp_filename and plp_filename.endswith(".plp"):
                        counter = 1
                        plp_filepath = os.path.join(plp_dir, plp_filename)
                        #print(plp_filepath)
                        with open(plp_filepath, 'r') as plp_file:
                            for res_line in plp_file:
                                if res_line[0] == '#':
                                    continue
                                else:
                                    posterior_prob = float(res_line.split()[4])
                                    res_index = int(res_line.split()[0])
                                    if posterior_prob >= 0.90 and res_index >= signalpept_val:   # the threshold for signal peptides
                                        cut_site = max((res_index - int(trimm_val[0])), 0)
                        break
                if counter == 0:
                    print('This sequence has not been found in the specified plp directory: ', seq_id, file=sys.stderr)
                with open(renameflnm, 'r') as rename_file:
                    for renem_line in rename_file:
                        if seq_id in renem_line:
                            new_name = renem_line.split(' ')[1].rstrip()
                            #print(new_name, list_boundaries_normal, list_boundaries_special, cut_site)
                            for tm in list_boundaries_normal:
                                for tm_res in range(tm[0], (tm[1]+1)):
                                    out_file.write( new_name + ' ' + str(tm_res - cut_site) + ' 0\n' )
                            for sp_tm in list_boundaries_special:
                                for sp_tm_res in range(sp_tm[0], (sp_tm[1]+1)):
                                    out_file.write( new_name + ' ' + str(sp_tm_res - cut_site) + ' 1\n' )
                            break
        
    else:
        with open(phobstdout, 'r') as phob_out, open(output_f, 'w') as out_file:
            for line in phob_out:
                seq_id = None
                list_boundaries_normal = None
                list_boundaries_special = []								#to not throw error afterwards
                if switch_col:
                    seq_id, list_boundaries_normal = phobius_short_pred_field_selecter(line, 'n')
                    tmp, list_boundaries_special = phobius_short_pred_field_selecter(line, 's')
                else:
                    seq_id, list_boundaries_normal = phobius_short_pred_field_selecter(line, '-')
                with open(renameflnm, 'r') as rename_file:
                    for renem_line in rename_file:
                        if seq_id in renem_line:
                            new_name = renem_line.split(' ')[1].rstrip()
                            #print(new_name, list_boundaries_normal, list_boundaries_special)
                            for tm in list_boundaries_normal:
                                for tm_res in range(tm[0], (tm[1]+1)):
                                    out_file.write( new_name + ' ' + str(tm_res) + ' 0\n' )
                            for sp_tm in list_boundaries_special:
                                for sp_tm_res in range(sp_tm[0], (sp_tm[1]+1)):
                                    out_file.write( new_name + ' ' + str(sp_tm_res) + ' 1\n' )
                            break



def main(args):
    #print(args.phob_stout, args.rename_file, args.output_file, args.phobius_mode, args.trimm_value, args.plp_dir, args.signalpept, args.two_colors)

    trimm_value = args.trimm_value
    if args.trimm_value:
        try:
            left_trimm_value = int(args.trimm_value.split(',', 1)[0])
            right_trimm_value = int(args.trimm_value.split(',', 1)[1])
        except:
            parser.error('--trimm_value has to be a string comprised by two ingers comma separated like 50,10   It was given  :  ', args.trimm_value)
        else:
            trimm_value = args.trimm_value.split(',')
            if not args.plp_dir:
                parser.error('--trimm_value requires --plp_dir to be specified')
        

    if args.phobius_mode == 'short':
        phobstout_short_colist_prep(args.phob_stout, args.rename_file, args.output_file, trimm_value, args.plp_dir, args.signalpept, args.two_colors)
    else:
        phobstout_colist_prep(args.phob_stout, args.rename_file, args.output_file, trimm_value, args.plp_dir, args.signalpept, args.two_colors)




if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='This script colour an alignment given the prediction of phobius. Now it has many fncionalities mainly due to the possible ways phobius can be run. It can work on a trimmed alignment, so it needs to go and recover that info. It also necessitates a rename file, the one used by tcoffee to rename sequences. The output file is a collection of aminoa acids positions along the sequence that have been predicted by phobius to be Transmembrane, the couors reflect the confidence of the model given by posterior probability scores.')
    parser.add_argument('--phob_stout', type=str, required=True, help='Mandatory flag, a phobius standard output redirected into a file, if the --trimm_value option is given, the ID lines must contain a name that is a keyword to identify a .plp file of that specific sequence. From the plp file the trimming positions will be recovered (recomputed)')
    parser.add_argument('--rename_file', type=str, required=True, help='Mandatory flag, two column file, space separeted, having the final sequence name in the second row. The first row is th ename of the sequence we want to change and it is what tcoffe will try to look for, (probably by string replacment like sed)')
    parser.add_argument('--output_file', type=str, required=True, help='Mandatory flag, the output filename. This file will be a three column one, space separated, where on the first column there is the new name of the sequence. On the second column the position of the aa that has to be coloured (this value is automatically adjusted in case of trimming). The third column is the numeric value 1-9 corresponding to the color intensity, being 9 over 90%% of posterior probability (dark blue).')
    parser.add_argument('--phobius_mode', default='short', choices=['short', 'long'], required=False, help='optional flag, default short, a string stating the type of phobius output, either long or short, referring to the two types of output of phobius')
    parser.add_argument('--trimm_value', type=str, required=False, help='trimm_value used to cut the sequences, the sequences have been cutted from another script trimm_multifasta.py, go check it to see how the cut is performed, trimm_value has to be the same number used for that script')
    parser.add_argument('--plp_dir', type=str, required=False, help='if the --trimm_value option is given, The directory in which search all .plp files of the sequences present in the --phob_stout file.')
    parser.add_argument('--signalpept', type=int, default=0, required=False, help='A signal peptide value used when trimming, default 0, good pp values will be ignored if their position on the sequence is lower than the value specified it must be an integer. Usefull when signal peptides are mistaken for TM by phobius.') 
    parser.add_argument('--two_colors', type=str, required=False, help='a boolean parameter used to switch to a different (two colouring) coloration mostly used when there aree two different types of transmembranes')

    args = parser.parse_args()
    main(args) 


