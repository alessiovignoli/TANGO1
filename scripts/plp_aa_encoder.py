#!/usr/bin/env python3

import sys

def hard_encoder(in_plp, label, max_iter=False):
    #print(in_plp, max_iter)
    encoder_dict = {'A':'000000000000000000000001', 'C':'000000000000000000000010', 'D':'000000000000000000000100', 'E':'000000000000000000001000', 'F':'000000000000000000010000', 'G':'000000000000000000100000', 'H':'000000000000000001000000', 'I':'000000000000000010000000', 'K':'000000000000000100000000',  'L':'000000000000001000000000', 'M':'000000000000010000000000', 'N':'000000000000100000000000', 'P':'000000000001000000000000', 'Q':'000000000010000000000000', 'R':'000000000100000000000000', 'S':'000000001000000000000000', 'T':'000000010000000000000000', 'V':'000000100000000000000000', 'W':'000001000000000000000000', 'Y':'000010000000000000000000', 'X':'000100000000000000000000', 'B':'001000000000000000000000', 'U':'010000000000000000000000', 'Z':'100000000000000000000000'}
    #for elem in encoder_dict:
    #    print(elem, ' ', encoder_dict[elem], ' ', len(encoder_dict[elem]))
    check = 1000000000
    if max_iter != False:
        check = max_iter
    with open(in_plp, 'r') as inplp:
        #print(inplp)
        print('protein_id position_of_aa_in_sequence amminoacid_code plp_cytopplasmic plp_non_cytoplasmic plp_normal_membrane plp_signal_peptide plp_special_helix plp_loop_intramembrane class_label')
        num_seq = -1
        for line in inplp:
            if num_seq == check:
                break
            if line[0] == '#':
                seq_id = line.split()[1]
                num_seq += 1
                column = inplp.readline()
                #print(seq_id)
            else:
                #print(line, end='')
                aa_type = line.split()[1]
                print(seq_id, line.split()[0], encoder_dict[aa_type], line.split()[2], line.split()[3], line.split()[4], line.split()[5], line.split()[6], line.split()[7], label_str)
                


if __name__ == "__main__":
    try:
        input_plp = sys.argv[1]
        label_str = sys.argv[2]
        max_iter = int(sys.argv[3])
    except:
        if label_str != None:
            hard_encoder(input_plp, label_str)
        else:
            print('Program usage: text.py <an input plp file either single ormulti plp> < a max iteration number, this is the number ofsequence toconsider in the case of a multi plp file>', file=sys.stderr)
            raise SystemExit
    else:
        hard_encoder(input_plp, label_str, max_iter)
