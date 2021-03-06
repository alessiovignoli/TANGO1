#!/usr/bin/env python3

# pp = posterior probability     always
# tm = trans membrane            always

# take a look at phobius_short_prediction_field_retriever.py  it is called to create the most important inputs for this script
# take a look at the print error for the call of the function



import sys
import sys
import os
sys.path.append(os.path.abspath("/home/alessio/Desktop/erasmus-internship/scripts"))
from phobius_short_prediction_field_retriever import phobius_short_pred_field_selecter



def  adjacent_annot_finder(in_domain_info, in_shortpred, input_field, prefix):
    list_of_range = []
    range_uncertainty_value = 5                     ## change this to deal with how much is used to consider adjacency
                                                    ## the left adj is considered to have the right margin at most range_uncertainty_value
                                                    ## form the right end ogf the region considered in short pred the right is the opposite
    with open(in_shortpred, 'r') as in_pred:
        for line in in_pred:
            seq_id, list_boundaries = phobius_short_pred_field_selecter(line, input_field)
            for single_boundary in list_boundaries:                 ## this changes a lot but here ram memory can be decreased
                                                                    ## code will be more complicated though not a for here but late
                single_boundary.append(seq_id.split('_')[-1])       ## done so that Uniref50_A0A6L2KD28 gives only the real id part
                list_of_range.append(single_boundary)
    #print(list_of_range)
    output_filenames = [(prefix + '_left.adj'), (prefix + '_center.adj'), (prefix + '_right.adj')]
    with open(in_domain_info, 'r') as in_domain, open(output_filenames[0], 'w') as left_adj, open(output_filenames[1], 'w') as center_adj, open(output_filenames[2], 'w') as right_adj:
        for domain_line in in_domain:
            for boundary in list_of_range:
                if boundary[2] == domain_line.split()[0].split('_')[-1]:           ## exact string match of prot ids
                    left_margin = boundary[0]
                    right_margin = boundary[1]
                    #print(boundary[2])
                    tmp_left_adjacent = 0 
                    buffer1 = [0, None]
                    tmp_right_adjacent = float('inf')
                    buffer2 = [0, None]
                    for annotation in domain_line.split('[')[1:]:
                        for i, left_extr in enumerate(annotation.split(';')[2::2]):
                            #print('i:', i, ' left extremity:', left_extr)           ## i has to be multiplied by two to represent the index
                            right_extr = None
                            try:
                                right_extr = int((annotation.split(';')[2:][(i*2+1)]).split(']')[0]) # it can end with the closed square brachet 
                                expedient = int(left_extr)		# it is rare but a comma might be inside the annotation name in this case
                                                                        # the script would crash and the best thing is go with the following code that's skips the annotation
                            except Exception:
                                #print(domain_line)			# this whole except section id made to catch the error given by a
                                #print('annotation:', annotation)	# open square brachet inside the annotation name
                                #print('boundary', boundary)		# solution is not to compute anything in that scenario 
                                #print('i:', i, ' left extremity:', left_extr)
                                #raise SystemExit
                                break
                            else:
                                pass
                            if int(left_extr) <= left_margin and right_extr >= right_margin:
                                #print('0: ', left_extr, right_extr,  annotation.split(';')[1], left_margin, right_margin)
                                center_adj.write((left_extr + ' ' + str(right_extr) + ' ' + str(left_margin) + ' ' + str(right_margin) + ' [;' + annotation.split(';')[1] + '\n' ))
                            elif right_extr < (right_margin - range_uncertainty_value) and right_extr >= tmp_left_adjacent:
                                tmp_left_adjacent = right_extr
                                buffer1 = [left_extr, annotation.split(';')[1]]
                            elif int(left_extr) > (left_margin + range_uncertainty_value) and int(left_extr) <= tmp_right_adjacent:
                                tmp_right_adjacent = int(left_extr)
                                buffer2 = [str(right_extr), annotation.split(';')[1]]
                    #print('-1: ', tmp_left_adjacent, left_margin, right_margin)
                    #print('+1: ', tmp_right_adjacent, left_margin, right_margin)
                    #print(buffer1, buffer2)
                    if buffer1[1] is not None:
                        left_adj.write((buffer1[0] + ' ' + str(tmp_left_adjacent) + ' ' + str(left_margin) + ' ' + str(right_margin) +  ' [;' + buffer1[1] + '\n' ))
                    if buffer2[1] is not None:
                        right_adj.write((str(tmp_right_adjacent) + ' ' + buffer2[0]  + ' ' + str(left_margin) + ' ' + str(right_margin) + ' [;' + buffer2[1] + '\n' ))
                else:
                    #print('shortpred id:', boundary[2], '  domain_info id:', domain_line.split()[0])
                    continue
                    


if __name__ == "__main__":
    try:
        input_domain_info_file = sys.argv[1]
        input_phobius_shortpred = sys.argv[2]
        input_field = sys.argv[3]
        prefix_output_name = sys.argv[4]
    except Exception:
        print('Program usage: text.py <a domain_info file that has on each line a protein annotation in a special format,\nlike:\nUPI00068AA065 [IPR036514;sgnh hydrolase superfamily;471;634] [IPR002656;acyltransferase 3 domain;5;321] [IPR043968;sgnh domain;400;626]\nwhere the first field is the protein id and then there are as many lists as there have been found annotations for this entry in uniprot, there is the keyword associted with each feature a word descriptive feature and the extremities of the feature, there can be more of them\nto generate this file the script   uniprot_rest_query.py    look into that for more details> <the second file is the short output of phobius, where on each line is present the prediction and position of the latter on the sequence, with this info the neighbouring annotations are found,\n####  WARNING  #####\n the ids are suppurted only in uniref and uniprot format, in the sense that they are split on _ > <a mandatory field letter for knowing which feature predicted to look at in short pred if you do not know what to put take alook at      phobius_short_prediction_field_retriever.py    eror messages> < a prefix for the three output files that are going to be written each one representing a type of adjacent annotation, mainly closest to left, incorporating, closest to the right they can not be present in a pgiven pro>', file=sys.stderr)
        raise SystemExit
    else:
        adjacent_annot_finder(input_domain_info_file, input_phobius_shortpred, input_field, prefix_output_name)
