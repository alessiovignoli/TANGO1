#!/usr/bin/env python3

import os
import sys

def trimmer(in_multifasta, in_plp_dir, out_trimmedfilename, cutoff_vals, signpept_cutoff, threshold=0.90):
    #print(in_multifasta, in_plp_dir, out_trimmedfilename, cutoff_val)
    plp_lstfile = os.listdir(in_plp_dir)
    #print(plp_lstfile)
    with open(in_multifasta, 'r') as in_fasta, open(out_trimmedfilename, 'w') as out_fasta:
        left_cut_site = False
        right_cut_site = False
        header_buffer = ''
        for line in in_fasta:
            if line[0] == '>':
                left_cut_site = False
                right_cut_site = False
                #print(line, end='')
                for filename in plp_lstfile:
                    if filename.endswith(".plp"):
                        seq_id = (filename.split('.')[0]).split(':')[0]
                        #print(seq_id)
                        if seq_id in line:
                            #print(line, end='')
                            #print(seq_id)
                            plp_filepath = in_plp_dir + '/' + filename.rstrip()
                            with open(plp_filepath, 'r') as plp:
                                #print(plp)
                                for res_line in plp:
                                    #print(res_line)
                                    if res_line[0] == '#':
                                        continue
                                    else:
                                        #print(res_line, end='')
                                        posterior_prob = float(res_line.split()[4])
                                        res_index = int(res_line.split()[0])
                                        #print(posterior_prob)
                                        if posterior_prob >= threshold and res_index >= signpept_cutoff:
                                            left_and_right_cuts = cutoff_vals.split(',')
                                            left_cut_val = int(left_and_right_cuts[0])
                                            left_cut_site = max((res_index - left_cut_val), 0)
                                            if len(left_and_right_cuts) == 2:
                                                right_cut_val = int(left_and_right_cuts[1])
                                                right_cut_site = max((res_index + right_cut_val), res_index)
                                            #print(cut_site)
                                #print(cut_site)
                            break
                header_buffer = line
            else:
                aa_line = line.rstrip()
                #print(left_cut_site, right_cut_site, header_buffer)
                if left_cut_site == False:
                    print('the sequence:', header_buffer, 'has not a vaiable ancor for the trimm, it may due for lack of phobius prdiction of TM, or no residue higher than the threshold and the signal peptide cut off values. Or final hypothesis the plp file corresponding to the sequencewas not found.')
                    #out_fasta.write(header_buffer)
                    #out_fasta.write(line)
                elif right_cut_site == False:
                    new_seq = aa_line[left_cut_site:] + '\n'
                    out_fasta.write(header_buffer)
                    out_fasta.write(new_seq)
                else:
                    new_seq = aa_line[left_cut_site:right_cut_site] + '\n'
                    out_fasta.write(header_buffer)
                    out_fasta.write(new_seq)

if __name__ == "__main__":
    signal_pept_cutoff = None
    try:
        in_multifastapath = sys.argv[1]
        in_plp_dirpath = sys.argv[2]
        out_trimmedpath = sys.argv[3]
        cut_off_values = sys.argv[4]
        signal_pept_cutoff = int(sys.argv[5])
        threshold_for_cut_reference = float(sys.argv[6])
    except Exception:
        if signal_pept_cutoff is not None:
            trimmer(in_multifastapath, in_plp_dirpath, out_trimmedpath, cut_off_values, signal_pept_cutoff)
        else:
            print('Program usage: text.py <multifasta file that we want to trim the sequences> <plp folder, that means posterior probablity, so the folder where all plp files outputed from phobius are stored> <path to the output filename> <cut off values needed for the trim, must be 2 integers divided by a comma (left cut off value,right cut off value) right cut off value can be missing, all residues before (last residue with posterior prob >= 0.9 - left cut off value) are not present in the output, this is done for each sequence in the input, also all residues after (last residue with posterior prob >= 0.9 + right cut off value) are not present if right cut off value is specified> <signal peptide threshold, the residue index below which posterior prob values are not considered, to avoid including signal peptides>\n###WARNING###\nif no residue is found to respect the conditions for being the position to anchor the trimm the resulting fasta will not contain such ID and sequence.', file=sys.stderr)
            raise SystemExit
    else:
        trimmer(in_multifastapath, in_plp_dirpath, out_trimmedpath, cut_off_values, signal_pept_cutoff, threshold_for_cut_reference)
