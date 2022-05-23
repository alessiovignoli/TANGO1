#!/usr/bin/env python3

# pp = posterior probability     always
# tm = trans membrane            always


# take a look at plp_average_computer.py beacuase one of the input of the script is actually the output
# of such script

# this scripts computed the sstandard deviation plp value of the specified column <field_keyword> and aa range <list_of_ranges>
# of an plp file, in this case a range must be specified manually

# this script also deals automatically with multi-plp files given that every protein prediction is separated by something like this
#       # sp|Q8BI84|TGO1_MOUSE
#       #pos    aa      i       o       M       S       m       l
# and that there is a related id in the output file of the plp_average_computer.py containing all the proteins present in the multiplp in the same order
# the average file should look something like this
#       0.850000 sp|Q8BI84|TGO1_MOUSE [10, 30] ecc.. ecc..
# where the important field are the first four with first the mean the second the id present in the plp and multi plp (in same order)
# and the third and fourth field constitute the range of aa that has produced such average value


import sys
from math import sqrt
#import os
#sys.path.append(os.path.abspath("/home/alessio/Desktop/erasmus-internship/scripts"))
#from phobius_short_prediction_field_retriever import phobius_short_pred_field_selecter



def plp_stdeviationer(single_plp_file, field_keyword, in_average):
    print("this script has been thought and preepared but has not been written, as   jun 03 2021", file=sys.stderr)
    check_keyword = {"i":2, "o":3, "-":4, "n":4, "c":5, "s":6, "l":7}
    if field_keyword not in check_keyword:
        print('please give a keyword for selecting the field \nit can be passed from command line with --KEYWORD "keyword" in nextflow or with third argument in python', file=sys.stderr)
        print("the allowed keywords are: [c, i, o, -, n, s, l] \n c = signal peptide, i = inside membrane(cytoplasm), o = outside membrane, - = helix (in phobius originalmodel) \n (only in phobius-M7or later) => -n- = normal-helix  -s- = special-helix and -l- = loop-inramembrane", file=sys.stderr)
        raise SystemExit




def multi_plp_stdeviationer(in_average, multi_plp_file, field_keyword, number_max_iter=False):
    #print(in_average, multi_plp_file, field_keyword, number_max_iter)
    check_keyword = {"i":2, "o":3, "-":4, "n":4, "c":5, "s":6, "l":7}
    if field_keyword not in check_keyword:
        print('please give a keyword for selecting the field \nit can be passed from command line --KEYWORD in nextflow or with third argument in python', file=sys.stderr)
        print("the allowed keywords are: [c, i, o, -, n, s, l] \n c = signal peptide, i = inside membrane(cytoplasm), o = outside membrane, - = helix (in phobius originalmodel) \n (only in phobius-M7or later) => -n- = normal-helix  -s- = special-helix and -l- = loop-inramembrane", file=sys.stderr)
        raise SystemExit
    with open(in_average, 'r') as inaverage, open(multi_plp_file, 'r') as inplp:
        first_iter_counter = True
        header = ''
        columns = ''
        iter_num = 0
        #print("#average      standard deviation      id      range of segment    max")
        for txtline in inaverage:
            if number_max_iter != False and iter_num == number_max_iter:
                break
            else:
                seq_id = txtline.strip().split()[1]
                #print(seq_id, type(seq_id))
                left_extr = int((txtline.strip().split()[2])[1:-1])
                right_extr = int((txtline.strip().split()[3])[:-1])
                #print(left_extr, type(left_extr), end = ' ')
                #print(right_extr)
                #print(txtline, end='')
                #header = ''
                #columns = ''                
                if first_iter_counter:
                    header = inplp.readline()
                    columns = inplp.readline()
                    first_iter_counter = False
                ok_go = False
                if seq_id in header:
                    #print("if seq_id in header section", seq_id, header, end='')
                    ok_go = True
                else:
                    for plpline in inplp:
                        if plpline[0] == '#':
                            if seq_id in plpline:
                                ok_go = True
                                header = plpline
                                columns = inplp.readline()
                                #print("if seq_id in plpline", seq_id, plpline, end='')
                                break
                            else:
                                print('the sequence ID obtained : ', seq_id, 'is not present in the plp file examined at line :', plpline, ' check the orders of the protein of the input files, they must have the same protein sequence order', file=sys.stderr)
                summatory = 0.0
                if ok_go:
                    #print("if ok_go section", seq_id, header, end='')
                    for plpline in inplp:
                        aa_num = int(plpline.split()[0])
                        if aa_num >= left_extr and aa_num <= right_extr:
                            #print(aa_num, end=' ')
                            requested_field = float((plpline.split())[check_keyword[field_keyword]])
                            #print(aa_num, requested_field)
                            summatory += (( requested_field - float(txtline.strip().split()[0]) )**2)
                            #print(summatory)
                            if aa_num == right_extr:
                                stddeviation = sqrt ( summatory / ( right_extr - left_extr + 1 ))
                                print(txtline.split()[0], stddeviation, txtline.split()[1], txtline.split()[2], txtline.split()[3], txtline.split()[4])
                                break



if __name__ == "__main__":
    try:
        switch = sys.argv[1]
    except:
        print('Program usage: text.py < a mandatory string argument that can have as arguments either "single" or "multi" , this flag changes the behaviour of the script from using  the function plp_stdeviationer with single or function multi_plp_averager with multi> <the field keyword that has to be used to retrieve the extremes of the phobius prediction forsuch label>', file=sys.stderr)
        raise SystemExit
    else:
        if switch == "single":
            try:
                plp_file = sys.argv[2]
                field_keyword = sys.argv[3]
                average_output_file = sys.argv[4]
            except:
                print('Program usage: text.py < a mandatory string argument that can have as arguments either "single" or "multi" , this flag changes the behaviour of the script from using  the function plp_stdeviationer with single or function multi_plp_stdeviationer with multi, in this case is single> <the field keyword that has to be used to retrieve the column of the plp of the phobius prediction for such label> < the plp_average_computer.py redirection output file, simply use the > sign to as file>', file=sys.stderr)
                raise SystemExit
            else:
                plp_stdeviationer(plp_file, field_keyword, average_output_file)
        elif switch == "multi":
            try:
                plp_file = sys.argv[2]
                field_keyword = sys.argv[3]
                average_output_file = sys.argv[4]
                num_max_iterations = int(sys.argv[5])
            except:
                if average_output_file != None:
                    multi_plp_stdeviationer(average_output_file, plp_file, field_keyword)
                else:
                    print('Program usage: text.py < a mandatory string argument that can have as arguments either "single" or "multi" , this flag changes the behaviour of the script from using  the function plp_stdeviationer with single or function multi_plp_stdeviationer with multii, in this case is multi> < the multi plp file that needs to be computed on > <the field keyword that has to be used to retrieve the column of the plp of the phobius prediction for such label> < the plp_average_computer.py redirection file, simply use the > to redirect stdout to a file > < an optional field for the maxumum number of iterations, iterations are intended as number of ranges found by the phobius_short_pred_field_selecter funtion, so if ranges for tm helices are selected an iteration is a single tm helix not a protein>', file=sys.stderr)
                    raise SystemExit
            else:
                multi_plp_stdeviationer(average_output_file, plp_file, field_keyword, num_max_iterations)
        else:
            print('Program usage: text.py < a mandatory string argument that can have as arguments either "single" or "multi" , this flag changes the behaviour of the script from using  the function plp_averager with single or function multi_plp_averager with multi> The argument passed is not allowed, it has been passed: ', switch, file=sys.stderr)
            raise SystemExit
