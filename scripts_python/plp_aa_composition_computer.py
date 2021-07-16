#!/usr/bin/env python3

# pp = posterior probability     always
# tm = trans membrane            always

# take a look at phobius_short_prediction_field_retriever.py  it is called to create the most important inputs for this script
# take a look at the print error for the call of the function

# this scripts creates an average plp value of the specified column <field_keyword> and aa range <list_of_ranges>
# of an plp file, in this case a range must be specified manually

# this script also deals automatically with multi-plp files given that every protein prediction is separated by something like this
#       # sp|Q8BI84|TGO1_MOUSE
#       #pos    aa      i       o       M       S       m       l
# and that there is a related txt phobius short output prediction file containing all the proteins present in the multiplp in the same order


import sys
import os
sys.path.append(os.path.abspath("/home/alessio/Desktop/erasmus-internship/scripts"))
from phobius_short_prediction_field_retriever import phobius_short_pred_field_selecter



def plp_aacomp(single_plp_file, field_keyword, range_of_aa=False):
    print("this script has been thought and preepared but has not been written, as   jun 28 2021", file=sys.stderr)
    check_keyword = {"i":2, "o":3, "-":4, "n":4, "c":5, "s":6, "l":7}
    if field_keyword not in check_keyword:
        print('please give a keyword for selecting the field \nit can be passed from command line with --KEYWORD "keyword" in nextflow or with third argument in python', file=sys.stderr)
        print("the allowed keywords are: [c, i, o, -, n, s, l] \n c = signal peptide, i = inside membrane(cytoplasm), o = outside membrane, - = helix (in phobius originalmodel) \n (only in phobius-M7or later) => -n- = normal-helix  -s- = special-helix and -l- = loop-inramembrane", file=sys.stderr)
        raise SystemExit

def multi_plp_aacomp(in_pred_txt, multi_plp_file, field_keyword, number_max_iter=False):
    #print(in_pred_txt, multi_plp_file, field_keyword, number_max_iter)
    check_keyword = {"i":2, "o":3, "-":4, "n":4, "c":5, "s":6, "l":7}
    label_strings = {"s":0, "n":1, "l":2, "i":3, "o":4, "-":1, "c":5}
    if field_keyword not in check_keyword:
        print('please give a keyword for selecting the field \nit can be passed from command line --KEYWORD in nextflow or with third argument in python', file=sys.stderr)
        print("the allowed keywords are: [c, i, o, -, n, s, l] \n c = signal peptide, i = inside membrane(cytoplasm), o = outside membrane, - = helix (in phobius originalmodel) \n (only in phobius-M7or later) => -n- = normal-helix  -s- = special-helix and -l- = loop-inramembrane", file=sys.stderr)
        raise SystemExit
    aa_dict = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0, 'X':0, 'B':0, 'U':0, 'Z':0}
    with open(in_pred_txt, 'r') as intxt, open(multi_plp_file, 'r') as inplp:
        #print('protein_id length_segment A C D E F G H I K L M N P Q R S T V W Y X B U Z class_label')
        first_iter_counter = True
        next_header = ''
        iter_num = 0
        for txtline in intxt:
            if number_max_iter != False and iter_num == number_max_iter:
                break
            else:
                seq_id, list_boundaries = phobius_short_pred_field_selecter(txtline, field_keyword)
                #print(seq_id, type(seq_id))
                list_boundaries.reverse()
                #print(list_boundaries, type(list_boundaries))
                header = ''
                columns = ''
                if first_iter_counter:
                    header = inplp.readline()
                    columns = inplp.readline()
                    first_iter_counter = False
                else:
                    header = next_header
                    columns = inplp.readline()
                #print('header:',header, end='')
                #print('columns :',columns, end='')
                aa_num = 0
                if seq_id in header:
                    for plpline in inplp:
                        #print(plpline, end='')
                        if '#' in plpline:
                            next_header = plpline
                            #print('next_header:',next_header, end=''
                            if list_boundaries == []:
                                break
                            else:
                                len_seg1 = aa_num - int(list_boundaries[-1][0]) + 1
                                print(seq_id, list_boundaries[0], len_seg1, end=' ')
                                for elem1 in aa_dict:
                                    print((int(aa_dict[elem1])/len_seg1), end=' ')
                                print(label_strings[field_keyword])
                                list_boundaries.pop()
                                iter_num += 1
                                aa_dict = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0, 'X':0, 'B':0, 'U':0, 'Z':0}
                                break
                        else:
                            aa_num = int(plpline.split()[0])
                            #print(aa_num, list_boundaries)
                            if list_boundaries == []:
                                continue
                            else:
                                left_extr = int(list_boundaries[-1][0])
                                right_extr = None
                                if list_boundaries[-1][1] == 'inf':
                                    right_extr = 9999999999
                                else:
                                    right_extr = int(list_boundaries[-1][1])
                                if aa_num >= left_extr and aa_num <= right_extr:
                                    #print(plpline, end='')
                                    requested_field = (plpline.split())[1]
                                    aa_dict[requested_field] += 1
                                    #print(seq_id, aa_num, requested_field)
                                    if aa_num == right_extr:
                                        len_seg = right_extr - left_extr + 1
                                        print(seq_id, list_boundaries[-1], len_seg, end=' ')
                                        for elem in aa_dict:
                                            print((int(aa_dict[elem])/len_seg), end=' ')
                                        print(label_strings[field_keyword])
                                        list_boundaries.pop()
                                        iter_num += 1
                                        aa_dict = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0, 'X':0, 'B':0, 'U':0, 'Z':0}
                        if number_max_iter != False and iter_num == number_max_iter:
                            break
                else:
                    print('the sequence ID obtained : ', seq_id, 'is not present in the plp file examined or check the orders of the protein of the input files, they must have the same protein sequence order', file=sys.stderr)
                    raise SystemExit


if __name__ == "__main__":
    try:
        switch = sys.argv[1]
    except:
        print('Program usage: text.py < a mandatory string argument that can have as arguments either "single" or "multi" , this flag changes the behaviour of the script from using  the function plp_aacomp with single or function multi_plp_aacomp with multi> <the field keyword that hastobe used to retrieve the extremes of the phobius prediction forsuch label>', file=sys.stderr)
        raise SystemExit
    else:
        if switch == "single":
            try:
                plp_file = sys.argv[2]
                field_keyword = sys.argv[3]
                range_of_aa = sys.argv[4]
            except:
                if field_keyword != None:
                    plp_aacomp(plp_file, field_keyword)
                else:
                    print('Program usage: text.py < a mandatory string argument that can have as arguments either "single" or "multi" , this flag changes the behaviour of the script from using  the function plp_aacomp with single or function multi_plp_aacomp with multi, in this case is single> < the multi plp file to used> <the field keyword that has to be used to retrieve the column of the plp of the phobius prediction for such label> < an optional field to declare the range of aa to be selected for the average, otherwise is going to be executed on the whole column>', file=sys.stderr)
                    raise SystemExit
            else:
                plp_aacomp(plp_file, field_keyword, range_of_aa)
        elif switch == "multi":
            try:
                plp_file = sys.argv[2]
                field_keyword = sys.argv[3]
                phobius_short_stdout_txt_file = sys.argv[4]
                num_max_iterations = int(sys.argv[5])
            except:
                if phobius_short_stdout_txt_file != None:
                    multi_plp_aacomp(phobius_short_stdout_txt_file, plp_file, field_keyword)
                else:
                    print('Program usage: text.py < a mandatory string argument that can have as arguments either "single" or "multi" , this flag changes the behaviour of the script from using  the function plp_aacomp with single or function multi_plp_aacomp with multii, in this case is multi> < the plp file to be used> <the field keyword that has to be used to retrieve the column of the plp of the phobius prediction for such label> < an optional field for the maxumum number of iterations, iterations are intended as number of ranges found by the phobius_short_pred_field_selecter funtion, so if ranges for tm helices are selected an iteration is a single tm helix not a protein>', file=sys.stderr)
                    raise SystemExit
            else:
                multi_plp_aacomp(phobius_short_stdout_txt_file, plp_file, field_keyword, num_max_iterations)
        else:
            print('Program usage: text.py < a mandatory string argument that can have as arguments either "single" or "multi" , this flag changes the behaviour of the script from using  the function plp_aacomp with single or function multi_plp_aacomp with multi> The argument passed is not allowed, it has been passed: ', switch, file=sys.stderr)
            raise SystemExit

