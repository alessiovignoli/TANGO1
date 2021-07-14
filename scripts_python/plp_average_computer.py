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



def plp_averager(single_plp_file, field_keyword, range_of_aa=False):
    print("this script has been thought and preepared but has not been written, as   may 26 2021", file=sys.stderr)
    check_keyword = {"i":2, "o":3, "-":4, "n":4, "c":5, "s":6, "l":7}
    if field_keyword not in check_keyword:
        print('please give a keyword for selecting the field \nit can be passed from command line with --KEYWORD "keyword" in nextflow or with third argument in python', file=sys.stderr)
        print("the allowed keywords are: [c, i, o, -, n, s, l] \n c = signal peptide, i = inside membrane(cytoplasm), o = outside membrane, - = helix (in phobius originalmodel) \n (only in phobius-M7or later) => -n- = normal-helix  -s- = special-helix and -l- = loop-inramembrane", file=sys.stderr)
        raise SystemExit




def multi_plp_averager(in_pred_txt, multi_plp_file, field_keyword, number_max_iter=False):
    #print(in_pred_txt, multi_plp_file, field_keyword, number_max_iter)
    check_keyword = {"i":2, "o":3, "-":4, "n":4, "c":5, "s":6, "l":7}
    if field_keyword not in check_keyword:
        print('please give a keyword for selecting the field \nit can be passed from command line --KEYWORD in nextflow or with third argument in python', file=sys.stderr)
        print("the allowed keywords are: [c, i, o, -, n, s, l] \n c = signal peptide, i = inside membrane(cytoplasm), o = outside membrane, - = helix (in phobius originalmodel) \n (only in phobius-M7or later) => -n- = normal-helix  -s- = special-helix and -l- = loop-inramembrane", file=sys.stderr)
        raise SystemExit
    with open(in_pred_txt, 'r') as intxt, open(multi_plp_file, 'r') as inplp:
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
                summatory = 0.0
                #summatory1 = 0.0        # comment it
                max_plp = 0.0
                if seq_id in header:
                    for plpline in inplp:
                        #print(plpline, end='')
                        if '#' in plpline:
                            next_header = plpline
                            #print('next_header:',next_header, end='')
                            break
                        else:
                            aa_num = int(plpline.split()[0])
                            #print(aa_num, list_boundaries)
                            if list_boundaries == []:
                                continue
                            else:
                                if aa_num >= float(list_boundaries[-1][0]) and aa_num <= float(list_boundaries[-1][1]):
                                    #print(plpline, end='')
                                    requested_field = float((plpline.split())[check_keyword[field_keyword]]) #+ float((plpline.split())[4])
                                    #requested_field1 = float((plpline.split())[4])   # comment it
                                    #print(aa_num, requested_field)
                                    summatory += requested_field
                                    #summatory1 += requested_field1     # comment it
                                    #print(summatory)
                                    max_plp = max(max_plp, requested_field)
                                    if aa_num == list_boundaries[-1][1]:
                                        average_plp = summatory / (list_boundaries[-1][1] - list_boundaries[-1][0] + 1)
                                        #average_plp1 = summatory1 / (list_boundaries[-1][1] - list_boundaries[-1][0] + 1)   # comment it
                                        print(average_plp, seq_id, list_boundaries[-1], max_plp)       # uncomment it
                                        #print(average_plp, average_plp1, seq_id, list_boundaries[-1], max_plp)        # comment it
                                        list_boundaries.pop()
                                        summatory = 0.0
                                        #summatory1 = 0.0
                                        iter_num += 1
                        if number_max_iter != False and iter_num == number_max_iter:
                            break
                else:
                    print('the sequence ID obtained : ', seq_id, 'is not present in the plp file examined or check the orders of the protein of the input files, they must have the same protein sequence order', file=sys.stderr)
                    raise SystemExit



if __name__ == "__main__":
    try:
        switch = sys.argv[1]
    except:
        print('Program usage: text.py < a mandatory string argument that can have as arguments either "single" or "multi" , this flag changes the behaviour of the script from using  the function plp_averager with single or function multi_plp_averager with multi> <the field keyword that hastobe used to retrieve the extremes of the phobius prediction forsuch label>', file=sys.stderr)
        raise SystemExit
    else:
        if switch == "single":
            try:
                plp_file = sys.argv[2]
                field_keyword = sys.argv[3]
                range_of_aa = sys.argv[4]
            except:
                if field_keyword != None:
                    plp_averager(plp_file, field_keyword)
                else:
                    print('Program usage: text.py < a mandatory string argument that can have as arguments either "single" or "multi" , this flag changes the behaviour of the script from using  the function plp_averager with single or function multi_plp_averager with multi, in this case is single> < the multi plp file to used> <the field keyword that has to be used to retrieve the column of the plp of the phobius prediction for such label> < an optional field to declare the range of aa to be selected for the average, otherwise is going to be executed on the whole column>', file=sys.stderr)
                    raise SystemExit
            else:
                plp_averager(plp_file, field_keyword, range_of_aa)
        elif switch == "multi":
            try:
                plp_file = sys.argv[2]
                field_keyword = sys.argv[3]
                phobius_short_stdout_txt_file = sys.argv[4]
                num_max_iterations = int(sys.argv[5])
            except:
                if phobius_short_stdout_txt_file != None:
                    multi_plp_averager(phobius_short_stdout_txt_file, plp_file, field_keyword)
                else:
                    print('Program usage: text.py < a mandatory string argument that can have as arguments either "single" or "multi" , this flag changes the behaviour of the script from using  the function plp_averager with single or function multi_plp_averager with multii, in this case is multi> < the plp file to be used> <the field keyword that has to be used to retrieve the column of the plp of the phobius prediction for such label> < an optional field for the maxumum number of iterations, iterations are intended as number of ranges found by the phobius_short_pred_field_selecter funtion, so if ranges for tm helices are selected an iteration is a single tm helix not a protein>', file=sys.stderr)
                    raise SystemExit
            else:
                multi_plp_averager(phobius_short_stdout_txt_file, plp_file, field_keyword, num_max_iterations)
        else:
            print('Program usage: text.py < a mandatory string argument that can have as arguments either "single" or "multi" , this flag changes the behaviour of the script from using  the function plp_averager with single or function multi_plp_averager with multi> The argument passed is not allowed, it has been passed: ', switch, file=sys.stderr)
            raise SystemExit
