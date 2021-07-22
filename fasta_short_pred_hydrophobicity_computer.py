#!/usr/bin/env python3

# pp = posterior probability     always
# tm = trans membrane            always

# take a look at phobius_short_prediction_field_retriever.py  it is called to create the most important inputs for this script
# take a look at the print error for the call of the function




import sys
import os
sys.path.append(os.path.abspath("/home/alessio/Desktop/erasmus-internship/scripts"))
from phobius_short_prediction_field_retriever import phobius_short_pred_field_selecter



def hydro_segment_computer(in_txt, in_fasta, field_keyword, scale, number_max_iter=False):
    #print(in_txt, in_fasta, field_keyword, number_max_iter)
    check_keyword = {"i":2, "o":3, "-":4, "n":4, "c":5, "s":6, "l":7}
    label_strings = {"s":0, "n":1, "l":2, "i":3, "o":4, "-":1, "c":5}
    if field_keyword not in check_keyword:
        print('please give a keyword for selecting the field \nit can be passed from command line --KEYWORD in nextflow or with third argument in python', file=sys.stderr)
        print("the allowed keywords are: [c, i, o, -, n, s, l] \n c = signal peptide, i = inside membrane(cytoplasm), o = outside membrane, - = helix (in phobius originalmodel) \n (only in phobius-M7or later) => -n- = normal-helix  -s- = special-helix and -l- = loop-inramembrane", file=sys.stderr)
        raise SystemExit
    kyte_doolittle_hydro_dict = {'A':1.8, 'C':2.5, 'D':-3.5, 'E':-3.5, 'F':2.8, 'G':-0.4, 'H':-3.2, 'I':4.5, 'K':-3.9, 'L':3.8, 'M':1.9, 'N':-3.5, 'P':-1.6, 'Q':-3.5, 'R':-4.5, 'S':-0.8, 'T':-0.7, 'V':4.2, 'W':-0.9, 'Y':-1.3, 'X':0.0, 'B':0.0, 'U':0.0, 'Z':0.0}
    GES_hydro_dict = {'A':1.6, 'C':2.0, 'D':-9.2, 'E':-8.2, 'F':3.7, 'G':1.0, 'H':-3.0, 'I':3.1, 'K':-8.8, 'L':2.8, 'M':3.4, 'N':-4.8, 'P':-0.2, 'Q':-4.1, 'R':12.3, 'S':0.6, 'T':1.2, 'V':2.6, 'W':1.9, 'Y':0.7, 'X':0.0, 'B':0.0, 'U':0.0, 'Z':0.0}
    to_be_used_dict = None
    if scale == "kyte":
        to_be_used_dict = kyte_doolittle_hydro_dict
    elif scale == "GES":
        to_be_used_dict = GES_hydro_dict
    else:
        print('please specify a valid argument for the selection of the Hydrophobicity scale, valid are "kyte" and "GES" it has been given:', scale, file=sys.stderr)
        raise SystemExit
    with open(in_txt, 'r') as intxt, open(in_fasta, 'r') as infasta:
        iter_num = 0
        for txtline in intxt:
            if number_max_iter != False and iter_num == number_max_iter:
                break
            else:
                seq_id, list_boundaries = phobius_short_pred_field_selecter(txtline, field_keyword)
                #print(seq_id, type(seq_id))
                #print(list_boundaries, type(list_boundaries))
                header = infasta.readline()
                if seq_id in header:
                    seq = infasta.readline().rstrip()
                    #print('bubba', seq_id, 'is in ', header)
                    for segment in list_boundaries:
                        summatory = 0
                        #print(segment)
                        left_extr = int(segment[0]) - 1
                        right_extr = None
                        if segment[1] == 'inf':
                            right_extr = len(seq)
                        else:
                            right_extr = int(segment[1])
                        for aa in seq[left_extr:right_extr]:
                            #print(aa, end='')
                            summatory += to_be_used_dict[aa]
                        #print((right_extr - left_extr))
                        average_hydro_segment = summatory / (right_extr - left_extr)
                        print(seq_id, segment, average_hydro_segment, label_strings[field_keyword])
                        iter_num += 1
                        if number_max_iter != False and iter_num == number_max_iter:
                            break
                else:
                    seq = infasta.readline()
                    print('the sequence ID obtained : ', seq_id, 'is not present in the fasta file examined or check the orders of the protein of the input files, they must have the same protein sequence order', file=sys.stderr)
                




if __name__ == "__main__":
    try:
        input_txt = sys.argv[1]
        input_fasta = sys.argv[2]
        field_keyword = sys.argv[3]
        hydr_scale_id = sys.argv[4]
        num_max_iterations = int(sys.argv[5])
    except:
        if field_keyword != None:
            hydro_segment_computer(input_txt, input_fasta, field_keyword, hydr_scale_id)
        else:
            print('Program usage: text.py < a txt file containing the redirection of phobius short prediction, where every line is the prediction for one protein> <the fasta or multifasta file that generated the txt file, it has to have the sequences on the same order of the prediction in the txt> <the field keyword that hastobe used to retrieve the extremes of the phobius prediction forsuch label> < the type of scale to use for the computation of hydrophob supported Kyte-doolittle and GES , given as arguments  kyte  and  GES> <the number max number of elemnts to be considered (predicted) not proteins>', file=sys.stderr)
            raise SystemExit
    else:
        hydro_segment_computer(input_txt, input_fasta, field_keyword, hydr_scale_id, num_max_iterations)
