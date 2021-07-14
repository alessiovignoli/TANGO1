#!/usr/bin/env python3

import sys

def counter_of_many_things(headers_file, phobius_short):
    with open(headers_file, 'r') as h_infile, open(phobius_short, 'r') as phob_infile:
        #print(h_infile, phob_infile)
        tm_found = 0
        tp = 0
        fn = 0
        fp = 0
        tn = 0
        special_counter = 0
        for header in h_infile:
            #with open(phobius_short, 'r') as phob_infile:
                for pred_line in phob_infile:
                    pred_id = pred_line.split(' ')[0]
                    #print(pred_id)
                    #print(header, end='')
                    if pred_id in header:
                        tm_found += 1
                        if '-n-' in pred_line:
                            #print(pred_id, header, pred_line, end='')
                            tp += 1
                        else:
                            #print(pred_id, header, pred_line, end='')
                            fn += 1
                        if '-s-'in pred_line:
                            #print(pred_id, header, pred_line, end='')
                            special_counter += 1
                        if tm_found != (tp + fn):
                            print(pred_id, header, pred_line, end='')
                            raise SystemExit
                        break
                    else:
                        if '-n-' in pred_line:
                            #print(pred_id, header, pred_line, end='')
                            fp += 1
                        else:
                            #print(pred_id, header, pred_line, end='')
                            tn += 1
                        if '-s-'in pred_line:
                            #print(pred_id, header, pred_line, end='')
                            special_counter += 1
        print('TM headers found in the queried file (second one) : ', tm_found)
        print('True posittive : ', tp)
        print('False negative : ', fn)
        print('False posittive : ', fp)
        print('True negative : ', tn)
        print('special helices found in the real TM ids : ', special_counter)

if __name__ == "__main__":
    try:
        headers_file_path = sys.argv[1]
        phobius_short_output_infile = sys.argv[2]
    except:
        print('Program usage: text.py <one header per line sorted -u file> <only second line of phobius_M7.pl -short output for many proteins, sorted -u>', file=sys.stderr)
        raise SystemExit
    else:
        counter_of_many_things(headers_file_path, phobius_short_output_infile)
