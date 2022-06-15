#!/usr/bin/env python3

# pp = posterior probability     always
# tm = trans membrane            always

import sys

def pp_colist_prep(ppflnm, renameflnm, trimm_val=False, signpept_val=False, switch_col=False, threshold=0.9):
    #print('bubba')
    #print(ppflnm, renameflnm, trimm_val, signpept_val, threshold)
    cut_site = None
    pp_file = open(ppflnm, 'r')
    left_trimm_val = 0
    right_trimm_val = 5000000               # just to set a very high number since after it could be sumed to an integer
    if trimm_val != False:
        left_right_list = trimm_val.split(',')
        left_trimm_val = int(left_right_list[0])
        if len(left_right_list) == 2:
            right_trimm_val = int(left_right_list[1])
        for res_line in pp_file:
            if res_line[0] == '#':
                continue
            else:
                posterior_prob = float(res_line.split()[4])
                res_index = int(res_line.split()[0])
                #print(posterior_prob)
                if posterior_prob >= threshold and res_index >= int(signpept_val):
                    cut_site = max((res_index - left_trimm_val), 0)
    else:
        cut_site = 0
    pp_file.close()
    #print(right_cut_site)
    with open(ppflnm, 'r') as pp_file, open(renameflnm, 'r') as rename_file:
        #print(pp_file, rename_file)
        #list_renames = rename_file.readlines()
        tmp = (pp_file.readline().split())[1]
        old_header = ''
        if ':' in tmp:
            old_header = tmp.split(':')[0] + '_' + tmp.split(':')[1]
        else:
            old_header = tmp
        #priint(old_header)
        new_header = ''
        for line in rename_file:
            if old_header in line:
                #print(line)
                new_header = (line.split()[1]).rstrip()
        #print(new_header)
        classes_pp_file = pp_file.readline()
        for line in pp_file:
            print(line)
            signal_pept_pp = float(line.split()[5])
            #print(signal_pept_pp)
            tm_pp = 0.0
            if switch_col is  False or switch_col == 'false' or switch_col == 'False':
                tm_pp = float(line.split()[4])
            else:
                tm_pp = float(line.split()[6])
            #print(tm_pp)
            aa_number = str(int(line.split()[0]) - cut_site)
            #print(aa_number)
            if signal_pept_pp != 0.0 or int(aa_number) < 0 or int(aa_number) >= (right_trimm_val + left_trimm_val):
                continue
            #print(line, end='')
            elif tm_pp > 0.9:
                #print(line, end='')
                print(new_header, aa_number, '9')
            elif tm_pp <= 0.1 and tm_pp >= 0.01:         # here the lower threshold is used
                print(new_header, aa_number, '0')
            elif tm_pp > 0.1 and tm_pp <= 0.2:
                print(new_header, aa_number, '1')
            elif tm_pp > 0.2 and tm_pp <= 0.3:
                print(new_header, aa_number, '2')
            elif tm_pp > 0.3 and tm_pp <= 0.4:
                print(new_header, aa_number, '3')
            elif tm_pp > 0.4 and tm_pp <= 0.5:
                print(new_header, aa_number, '4')
            elif tm_pp > 0.5 and tm_pp <= 0.6:
                print(new_header, aa_number, '5')
            elif tm_pp > 0.6 and tm_pp <= 0.7:
                print(new_header, aa_number, '6')
            elif tm_pp > 0.7 and tm_pp <= 0.8:
                print(new_header, aa_number, '7')
            elif tm_pp > 0.8 and tm_pp <= 0.9:
                print(new_header, aa_number, '8')


if __name__ == "__main__":
    posterior_prob_filename = None
    rename_filename = None
    trimm_value = None
    signal_pept_cutoff = None
    special_tm_switch = None 
    try:
        posterior_prob_filename = sys.argv[1]
        rename_filename = sys.argv[2]
        trimm_value = (sys.argv[3])
        signal_pept_cutoff = sys.argv[4]
        special_tm_switch = sys.argv[5]
        possible_threshold = float(sys.argv[6])
    except Exception:
        if posterior_prob_filename is not None and rename_filename is not None:
            if trimm_value is None:
                pp_colist_prep(posterior_prob_filename, rename_filename)
            elif trimm_value is not None and signal_pept_cutoff is None:
                pp_colist_prep(posterior_prob_filename, rename_filename, False, False, trimm_value)
            elif trimm_value is not None and special_tm_switch is None:
                pp_colist_prep(posterior_prob_filename, rename_filename, trimm_value, signal_pept_cutoff)
            else:
                pp_colist_prep(posterior_prob_filename, rename_filename, trimm_value, signal_pept_cutoff, special_tm_switch)
        else:
            print('Program usage: text.py <a postirior probability file, phobius -plp option output file> <two column file, space separeted, having the final ids in the second row> < <trimm_value used to cut the sequences, the sequences have been cutted from another script trimm_multifasta.py, go check it to see how the cut is performed, trimm_value has to be the same number used for that script> <a boolean parameter used to switch to a different set of aa to color, mostly used when there are two different types of transmembranes> <optional field to set the threshold, give the number as 0.70 no comma, by default set to 0.9, this is used to consider where the sequence has found the reference for the cut, the last aa wuth plp value bigger than threshold>', file=sys.stderr)
            raise SystemExit
    else:
        pp_colist_prep(posterior_prob_filename, rename_filename, trimm_value, signal_pept_cutoff, special_tm_switch, possible_threshold)
