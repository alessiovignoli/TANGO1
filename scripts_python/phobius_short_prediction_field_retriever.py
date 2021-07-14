#!/usr/bin/env python3


# This script is not thought to be launche on a file but rather on a phobius short output line
# this script is written to be used inside another script block that has already opemed the input file

import sys

def phobius_short_pred_field_selecter(short_pred_line, keyword):
    #check_keyword = {"i":2, "o":3, "-":4, "n":4, "c":5, "s":6, "l":7}
    check_keyword = {"i":"i", "o":"o", "-":"-", "n":"-n-", "c":"c", "s":"-s-", "l":"-l-"}
    #print(keyword)
    #print(check_keyword)
    if keyword == False or keyword not in check_keyword:
        print('please give a keyword for selecting the field \nit can be passed from command line with --KEYWORD "keyword"', file=sys.stderr)
        print("the allowed keywords are: [c, i, o, -, n, s, l] \n c = signal peptide, i = inside membrane(cytoplasm), o = outside membrane, - = helix (in phobius originalmodel) \n (only in phobius-M7or later) => -n- = normal-helix  -s- = special-helix and -l- = loop-inramembrane", file=sys.stderr)
        raise SystemExit
    seqid = short_pred_line.strip().split(' ')[0]
    pre_section = short_pred_line.strip().split(' ')[-1]
    list_extr = []
    #print(seqid, pre_section)
    if keyword not in pre_section:
        print(seqid, " does not have the requested field present in it's prediction,   requested fieldkeyword =", keyword,  file=sys.stderr)
        return seqid, list_extr
        #raise SystemExit
    if keyword == '-' and '-n-' in pre_section:
        print(seqid, " is required to be split using '-' as keyword but '-n-' has been found in the prediction, '-' is used as spatiator in phobius.pl for helices, why are you doing it on a phobius-M7.pl short prediction?", file=sys.stderr)
        raise SystemExit
    fields = pre_section.split(check_keyword[keyword])
    if keyword == '-' and 'n' in fields[0]:
        fields = fields[1:]
    #print(fields)
    #list_extr = []
    for i in range (0, (len(fields)-1)):
        left_extr = 1
        right_extr = 'inf'
        extr_couple = []
        if '' == fields[i]:
            #print(seqid, left_extr, end=' ')
            extr_couple.append(left_extr)
            if '' == fields[(i+1)]:
                #print(right_extr)
                extr_couple.append(right_extr)
                #print(extr_couple)
                list_extr.append(extr_couple)
                break
            for char in range(0, len(fields[i+1])):
                try:
                    temp = int(fields[i+1][char])
                    #print("try section tmp =", tmp)
                except:
                    #print("except section", fields[i+1][char], type(fields[i+1][char]))
                    right_extr = int(fields[i+1][:char])
                    #print(right_extr)
                    extr_couple.append(right_extr)
                    #print(extr_couple)
                    list_extr.append(extr_couple)
                    break
        else:
            for digi in range((len(fields[i])-1), -1, -1):
                #print(fields[i][digi], type(fields[i][digi]))
                try:
                    tmp = int(fields[i][digi])
                    #print("try section tmp =", tmp)
                except:
                    #print("except section", fields[i][digi], type(fields[i][digi]))
                    left_extr = int(fields[i][(digi+1):])
                    #print(left_extr, end=' ')
                    extr_couple.append(left_extr)
                    #print(extr_couple)
                    break
            if '' == fields[(i+1)]:
                #print(right_extr)
                extr_couple.append(right_extr)
                #print(extr_couple)
                list_extr.append(extr_couple)
                break
            for char in range(0, len(fields[i+1])):
                #print(fields[i+1][char], type(fields[i+1][char]))
                try:
                    temp = int(fields[i+1][char])
                    #print("try section tmp =", tmp)
                except:
                    #print("except section", fields[i+1][char], type(fields[i+1][char]))
                    right_extr = int(fields[i+1][:char])
                    #print(right_extr)
                    extr_couple.append(right_extr)
                    #print(extr_couple)
                    list_extr.append(extr_couple)
                    break
    #print(list_extr)
    return seqid, list_extr


if __name__ == "__main__":
    try:
        phobius_short_pred_line = sys.argv[1]
        field_keyword = sys.argv[2]
    except:
        print('Program usage: text.py <a phobius short output line> <the field keyword that hastobe used to retrieve the extremes of the phobius prediction forsuch label>', file=sys.stderr)
        raise SystemExit
    else:
        phobius_short_pred_field_selecter(phobius_short_pred_line, field_keyword)
