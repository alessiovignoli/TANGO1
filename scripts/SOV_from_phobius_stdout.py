#!/usr/bin/env python3

#tm = Transmembrane always.

import sys

def sov_computer(reference_filepath, prediction_filepath, ref_keyword, pred_keyword, special_switch=False):
    with open(reference_filepath, 'r') as ref_file: #, open(prediction_filepath, 'r') as pred_file:
        #print(ref_file, pred_file)
        seq_count = 0
        ref_tm_range_list = []
        pred_tm_range_list = []
        SOV_seq = 0
        special_switch_count = 0
        for ref_line in ref_file:
            if 'ID' in ref_line:
                #print('ref_line', ref_line, end='')
                seq_count += 1
                seq_check = 0
                ref_seqid = (ref_line.split()[1]).rstrip()
                #print('ref_seqid', ref_seqid)
                ref_tm_range_list = []
                pred_tm_range_list = []
                found_seqid_conter = False
                with open(prediction_filepath, 'r') as pred_file:
                    #pred_id_line = pred_file.readline()
                    #print(pred_id_line, end='')
                    while seq_check == 0:
                        pred_line = (pred_file.readline()).rstrip()
                        #print('pred_line', pred_line, seq_check)
                        if pred_line == '':
                            print('A problem occured while parsing:\n', prediction_filepath, '\ntwo things could have occurred:\nfirst the file has less ID-lines or sequence than the reference file (first specified) \nsecond the file does not have // sequence-terminating lines in it', file=sys.stderr)
                            raise SystemExit
                        elif '//' in pred_line and found_seqid_conter == True:
                            seq_check = 1
                        elif ref_seqid in pred_line:
                            found_seqid_conter = True
                            #print(ref_seqid, pred_line, found_seqid_conter)
                        elif pred_keyword in pred_line and found_seqid_conter == True:
                            #print('pred_line', pred_line)
                            pred_tm_start = int(pred_line[14:21].strip())
                            pred_tm_end = int(pred_line[21:28].strip())
                            pred_tm_range_list.append((pred_tm_start, pred_tm_end))
                            if special_switch != False:
                                special_switch_count += 1
                #print(pred_tm_range_list)
            elif ref_keyword in ref_line:
                #print('ref_keywod', ref_line, end='')
                ref_tm_start = int(ref_line[14:21].strip())
                ref_tm_end = int(ref_line[21:28].strip())
                ref_tm_range_list.append((ref_tm_start, ref_tm_end))
            elif '//' in ref_line:
                #print('ref_line_//', ref_line, end='')
                #print(ref_tm_range_list)
                #print(pred_tm_range_list)
                if len(ref_tm_range_list) == 0 and len(pred_tm_range_list) != 0:
                    SOV_seq += 0
                    #print('ref_tm_range_list is empty')
                elif len(ref_tm_range_list) != 0 and len(pred_tm_range_list) == 0:
                    SOV_seq += 0
                    #print('pred_tm_range_list is empty')
                elif len(ref_tm_range_list) == 0 and len(pred_tm_range_list) == 0:
                    SOV_seq += 100
                    #print(pred_id_line, end='')
                    #print(SOV_seq)
                    #print('ref_tm_range_list is empty and pred_tm_range_list is empty')
                else:
                    summation = 0
                    N_tm_res = 0
                    overlap_conter = 0
                    for ref_range in ref_tm_range_list:
                        overlap_conter = 0
                        ref_range_start = ref_range[0]
                        ref_range_end = ref_range[1]
                        len_ref_range = ref_range_end - ref_range_start + 1
                        #N_tm_res += len_ref_range
                        #print('ref_ranges : ', ref_range_start, ref_range_end, N_tm_res)
                        for pred_range in pred_tm_range_list:
                            pred_range_start = pred_range[0]
                            pred_range_end = pred_range[1]
                            len_pred_range = pred_range_end - pred_range_start + 1
                            #print('pred ranges : ', pred_range_start, pred_range_end, len_pred_range)
                            if ref_range_start <= pred_range_start and pred_range_start <= ref_range_end <= pred_range_end:
                                #print('<if> ref Tm :', ref_range_start, ref_range_end, '   pred Tm :', pred_range_start, pred_range_end)
                                maxov = max(ref_range_end, pred_range_end) - min(ref_range_start, pred_range_start) + 1
                                minov = min(ref_range_end, pred_range_end) - max(ref_range_start, pred_range_start) + 1
                                delta = min((maxov - minov), minov, (len_ref_range/2), (len_pred_range/2))
                                summation += ((minov + delta)/maxov *len_ref_range)
                                #print((maxov - minov), maxov, minov, (len_ref_range/2), (len_pred_range/2))
                                #print('sum =', summation)
                                N_tm_res += len_ref_range
                                overlap_conter += 1
                            elif pred_range_start <= ref_range_start <= pred_range_end and pred_range_end <= ref_range_end:
                                #print('<First elif> ref Tm :', ref_range_start, ref_range_end, '   pred Tm :', pred_range_start, pred_range_end)
                                maxov = max(ref_range_end, pred_range_end) - min(ref_range_start, pred_range_start) + 1
                                minov = min(ref_range_end, pred_range_end) - max(ref_range_start, pred_range_start) + 1
                                delta = min((maxov - minov), minov, (len_ref_range/2), (len_pred_range/2))
                                summation += ((minov + delta)/maxov *len_ref_range)
                                #print((maxov - minov), maxov, minov, (len_ref_range/2), (len_pred_range/2))
                                #print('sum =', summation)
                                N_tm_res += len_ref_range
                                overlap_conter += 1
                            elif ref_range_start < pred_range_start and pred_range_end < ref_range_end:
                                #print('<Second elif> ref Tm :', ref_range_start, ref_range_end, '   pred Tm :', pred_range_start, pred_range_end)
                                maxov = max(ref_range_end, pred_range_end) - min(ref_range_start, pred_range_start) + 1
                                minov = min(ref_range_end, pred_range_end) - max(ref_range_start, pred_range_start) + 1
                                delta = min((maxov - minov), minov, (len_ref_range/2), (len_pred_range/2))
                                summation += ((minov + delta)/maxov *len_ref_range)
                                #print((maxov - minov), maxov, minov, (len_ref_range/2), (len_pred_range/2))
                                #print('sum =', summation)
                                N_tm_res += len_ref_range
                                overlap_conter += 1
                            elif pred_range_start < ref_range_start and ref_range_end < pred_range_end:
                                #print('<Third elif> ref Tm :', ref_range_start, ref_range_end, '   pred Tm :', pred_range_start, pred_range_end)
                                maxov = max(ref_range_end, pred_range_end) - min(ref_range_start, pred_range_start) + 1
                                minov = min(ref_range_end, pred_range_end) - max(ref_range_start, pred_range_start) + 1
                                delta = min((maxov - minov), minov, (len_ref_range/2), (len_pred_range/2))
                                summation += ((minov + delta)/maxov *len_ref_range)
                                #print((maxov - minov), maxov, minov, (len_ref_range/2), (len_pred_range/2))
                                #print('sum =', summation)
                                N_tm_res += len_ref_range
                                overlap_conter += 1
                            else:
                                #if special_switch == False:
                                #N_tm_res += len_ref_range
                                #print('no overlap')
                                continue
                        if overlap_conter == 0 and special_switch == False:
                            N_tm_res += len_ref_range
                            #print('no overlap')
                    if N_tm_res == 0:
                        SOV_seq += 0
                    else:
                        SOV_seq += (100*(1/N_tm_res)*(summation))
                        #print('SOV per seq : ', SOV_seq)
                        #print('N_tm_res : ', N_tm_res)
                #print(SOV_seq)
            else:
                continue
                #print(ref_line, end='')
        #check_on_ids = pred_file.readline()
        #if 'ID' in check_on_ids:
        #    print(prediction_filepath, '\n  has at least one sequence more than:\n', reference_filepath, '\nbe sure to compare reference and new model prediction files containing same amount of sequences')
        #    raise SystemExit
        #else:
        if special_switch == False:
            SOV = SOV_seq/seq_count
            #print(SOV_seq)
            print('number of sequences = ', seq_count)
            print('SOV = ', SOV, '\n')
        else:
            if special_switch_count == 0:
                SOV = 0
                print('special count = ', special_switch_count)
                print('SOV = ', SOV, '\n')
            else:
                if special_switch == 'true':
                    SOV = SOV_seq/special_switch_count
                    #print(SOV_seq)
                    print('special count = ', special_switch_count)
                    print('SOV = ', SOV, '\n')
                else:
                    SOV = SOV_seq/int(special_switch)
                    #print(SOV_seq)
                    print('special count = ', special_switch_count)
                    print('number of sequences = ', special_switch)
                    print('SOV = ', SOV, '\n')



if __name__ == "__main__":
    special_case_switch = None
    structure_pred_keyword = None
    try:
        reference_filename = sys.argv[1]
        prediction_filename = sys.argv[2]
        structure_ref_keyword = sys.argv[3]
        structure_pred_keyword = sys.argv[4]
        special_case_switch = sys.argv[5]
    except:
        if special_case_switch is None and structure_pred_keyword is not None:
            sov_computer(reference_filename, prediction_filename, structure_ref_keyword, structure_pred_keyword)
        else:
            print('Program usage: text.py <reference model phobius default stdoutput re-diredted as it is in a file (it can be made from many fasta sequence concatenated outputs)> <same thing as before but this time the phobius stdout should be from another model prediction (since the SOV is computed between a reference an a newmodel prediction)> <The third field has to be a keyword used from phobius to point out a specific region in the reference file (most of the times is DOMAIN or TRANSMEM)> <The fourth field is again a keyword used from phobius to point out a specific region in the predicted file> <the fifth field is optional but is used to switch to a sligtly different handle of the divider of the SOV_seq, this is a consequence of having the original phobius model predicting some special helices in some sequemnces and not in others, an annotation step is required (transfer_annotation) in order to have an estimate of where the special helix starts and ends in a reference sequence that actually do not have it predicted>', file=sys.stderr)
            raise SystemExit
    else:
        sov_computer(reference_filename, prediction_filename, structure_ref_keyword, structure_pred_keyword, special_case_switch)
