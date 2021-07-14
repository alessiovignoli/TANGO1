#!/usr/bin/env python3

# pp = posterior probability     always
# tm = trans membrane            always

import sys
import os

def phobstout_colist_prep(aln_flnm, phobstdout, plp_dir=False, trimm_val=False, signalpept_val=False):
    #print(aln_flnm, phobstdout, plp_dir, trimm_val, signalpept_val)
    dict_seqids_trimpos = {}
    one_line_all_seqids_aln = ''
    reference_seqid = ''                    # remember it has to be the last id in the alignment
    ref_tm_start = None
    ref_tm_len = None
    tmp_aln_file = open(aln_flnm, 'r')
    first_tmp_aln_file = tmp_aln_file.readline()
    second_tmp_aln_file = tmp_aln_file.readline()
    for tmp_line in tmp_aln_file:
        if tmp_line[0].isspace():
            break
        else:
            one_line_all_seqids_aln += (tmp_line.split(' ')[0])
            reference_seqid = tmp_line.split(' ')[0]
    tmp_aln_file.close()
    #print(reference_seqid)
    with open(phobstdout, 'r') as phob_out:
        counter = 0    # used as a check for the presence of a plp file with keyword inside it
        plp_keyword = ''
        cut_site = 0
        for phob_line in phob_out:
            if 'ID' in phob_line:
                cut_site = 0
                counter = 0
                plp_keyword = phob_line.split(' ')[3].rstrip()
                #print(plp_keyword)
            elif 'TRANSMEM' in phob_line:
                counter += 1
                if plp_keyword in reference_seqid and counter == 1:
                    ref_tm_start = int(phob_line[14:21].strip())
                    ref_tm_end = int(phob_line[21:28].strip()) + 1
                    ref_tm_len = ref_tm_end - ref_tm_start - 1
            elif ('//' in phob_line and counter < 2) or (plp_keyword in reference_seqid):
                if plp_dir != False:
                    list_of_plps = os.listdir(plp_dir)
                    for plp_filename in list_of_plps:
                        if plp_keyword in plp_filename and plp_filename.endswith(".plp") and plp_keyword in one_line_all_seqids_aln:
                            #print(plp_filename)
                            plp_filepath = plp_dir + plp_filename
                            with open(plp_filepath, 'r') as plp_file:
                                for res_line in plp_file:
                                    if res_line[0] == '#':
                                        continue
                                    else:
                                        posterior_prob = float(res_line.split()[4])
                                        res_index = int(res_line.split()[0])
                                        #print(posterior_prob)
                                        if posterior_prob >= 0.90 and res_index >= signalpept_val:   # the threshold for signal peptides
                                            left_cut_val = int(trimm_val.split(',')[0])
                                            cut_site = max((res_index - left_cut_val + 1), 0)
                            #print(plp_keyword, cut_site)
                            dict_seqids_trimpos[plp_keyword] = cut_site
                            break
                else:
                    if plp_keyword in one_line_all_seqids_aln:
                        #print(plp_keyword, cut_site)
                        dict_seqids_trimpos[plp_keyword] = cut_site
    #print(one_line_all_seqids_aln)
    #print(reference_seqid)
    #print(dict_seqids_trimpos)
    #print(reference_seqid, ref_tm_start, ref_tm_len)    
    with open(aln_flnm, 'r') as alnfile:
        firstline = alnfile.readline()
        secondline = alnfile.readline()
        aln_matrix = []
        seq_id_l = []
        mcounter = 0
        end_counter = 1
        #print(end_counter, (len(dict_seqids_trimpos)))
        for alnline in alnfile:
            if end_counter != (len(dict_seqids_trimpos)):
                #print(end_counter, (len(dict_seqids_trimpos)))
                seqid = (alnline.split(' ')[0])
                if alnline == '\n' and mcounter == 0:
                    mcounter = 1
                    aa_pos = 0
                    #print(end_counter, (len(dict_seqids_trimpos) - 1))
                    while aa_pos < len(aln_matrix[0]) and end_counter <= (len(dict_seqids_trimpos) - 1):
                        #print(aa_pos)
                        seq_n = 0
                        while seq_n < len(aln_matrix) and end_counter <= (len(dict_seqids_trimpos) - 1):
                            if dict_seqids_trimpos[reference_seqid] == ref_tm_start:
                                print('ID  ', (seq_id_l[seq_n]))
                                print('FT   TRANSMEM  ', (dict_seqids_trimpos[(seq_id_l[seq_n])]), ' ', (dict_seqids_trimpos[(seq_id_l[seq_n])] + ref_tm_len))
                                print('//')
                                end_counter += 1
                                #print(end_counter)
                            elif aln_matrix[seq_n][aa_pos] != '-':
                                #print(aln_matrix[seq_n][aa_pos])
                                dict_seqids_trimpos[(seq_id_l[seq_n])] += 1
                            seq_n += 1
                        aa_pos += 1
                    aln_matrix = []
                else:
                    if seqid in dict_seqids_trimpos or seqid == reference_seqid:
                        mcounter = 0
                        #print(seqid)
                        aln_aa_row = alnline.split()[1]
                        #print(aln_aa_row)
                        aln_matrix.append(aln_aa_row)
                        seq_id_l.append(seqid)
            else:
                return

if __name__ == "__main__":
    phob_stout_filename = None
    input_plp_dir = None
    trimm_value = None
    signalpept_value = None
    try:
        aln_file = sys.argv[1]
        phob_stout_filename = sys.argv[2]
        input_plp_dir = sys.argv[3]
        trimm_value = sys.argv[4]
        signalpept_value = int(sys.argv[5])
    except:
        if trimm_value is None and signalpept_value is None and input_plp_dir is None and phob_stout_filename is not None:
            phobstout_colist_prep(aln_file, phob_stout_filename)
        else:
            print('<alignment generated by t-coffe using (some of) the sequences found in the phob_stout_filename, however the sequence identifiers must be the same exact string of thoose used in the phob_stout_filename ID line (exact string matching used in thescript)> <a phobius default (-long flag) standard output redirected into a file, if the trimm value is specified, the ID lines must contain a name that is a keyword to identify a .plp file of that specific sequence> <if the trimm value option is given, The directory in which search all .plp files of the sequences present in the aln_file input file for this script> <trimm_value used to cut the sequences, this enebles all previous options, the sequences have been cutted from another script trimm_multifasta.py, go check it to see how the cut is performed, trimm_value has to be the same number used for that script> <A signal peptide value used when trimming good pp values will be ignored if their position on the sequence is lower than the value specified it must be an integer>', file=sys.stderr)
            raise SystemExit
    else:
        phobstout_colist_prep(aln_file, phob_stout_filename, input_plp_dir, trimm_value, signalpept_value)
