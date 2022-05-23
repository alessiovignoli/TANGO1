#!/usr/bin/env python3

import sys
import requests

def input_fasta_prepare(in_fasta):
    list_of_uniprot_http = []
    list_of_uniparc_http = []
    with open (in_fasta, 'r') as in_fa:
        for line in in_fa:
            if line[0] == '>':
                tmp_identifier = (line.split()[0])
                identifier = ''
                if '_' in tmp_identifier:
                    identifier = tmp_identifier.split('_')[1]
                else:
                    identifier = tmp_identifier[1:]
                if identifier[0:3] == 'UPI':
                    http_to_request = 'https://www.uniprot.org/uniparc/' + identifier + '.xml'
                    #print(http_to_request)
                    list_of_uniparc_http.append(http_to_request)
                else:
                    http_to_request = 'https://www.uniprot.org/uniprot/' + identifier + '.gff'
                    #print(http_to_request)
                    list_of_uniprot_http.append(http_to_request)
    #print(list_of_http)
    return list_of_uniprot_http, list_of_uniparc_http


def uniparc_rest_dealer(uniparc_https, out_file):
    with open(out_file, 'a') as outfile:
        for url in uniparc_https:
            resp = requests.get(url)
	    #print(type(resp.content))
            tmp_str = str(resp.content)
            if 'signatureSequenceMatch ' in tmp_str:
                tmp_list = tmp_str.split('signatureSequenceMatch')
                #print(tmp_list[1:])
                domains_of_id = (url.split('/')[-1]).split('.')[0]          ## protein id creation
                for substring in tmp_list[1:-1]:
                    #print(substring)
                    if 'ipr id' in substring:
                        domains_of_id +=  ' [' + substring.split('"')[5] + ';' + substring.split('"')[7].lower() + ';'  # made for unifing all digits in decriptive names of domains and so on to lower case ids are left as they are 
                        for elem in substring.split('"')[9::2]:
                            #print(elem, end=' ')
                            domains_of_id += elem
                            domains_of_id += ';'
                        domains_of_id = domains_of_id[:-1] + ']' 
                    elif '>\\n<' ==  substring:
                        #print('elif', substring)
                        continue
                    else:
                        print(substring)
                        raise SystemExit
                outfile.write(domains_of_id + '\n')
            else:
                print(url, file=sys.stderr)
                #print(resp.content, file=sys.stderr)
                #raise SystemExit


def uniprot_rest_dealer(uniprot_https, out_file):
    with open(out_file, 'w') as outfile:
        for url in uniprot_https:
            response = requests.get(url)
            #print(response.content)
            tmp_string = str(response.content)
            if '##gff-version' in tmp_string:
                #print(tmp_string)
                line_tmp_dict = {}
                protein_id = url.split('/')[-1].split('.')[0]
                for gff_line in tmp_string.split('\\n')[2:-1]:
                    if len(gff_line.split('\\t')) != 10:
                        print('this line does not have 10 fields', len(url, gff_line, file=sys.stderr))
                        #print(gff_line.split('\\t')[2], '|',  gff_line.split('\\t')[8])
                    #print(len(gff_line.split('\\t')))
                    key_for_dict = (gff_line.split('\\t')[2]).lower()		# set to lowercase digits for string match simplicity later on
                    feature_id = ''							# for thoose that have an ID it is going to be one or more space separated ids
                    descriptive_field = gff_line.split('\\t')[8]
                    if descriptive_field[:5] == 'Note=':
                        key_for_dict += (' ' + (descriptive_field.split(';')[0]).split('=')[1]).lower()		# set to lowercase digits for string match simplicity later on
                        if '|' in descriptive_field:
                            feature_id += (descriptive_field.split('|')[-1]).split(';')[0]
                        #print(key_for_dict, feature_id)
                    try:
                        line_tmp_dict[key_for_dict] += (';' + gff_line.split('\\t')[3] + ';' + gff_line.split('\\t')[4] + ';' + feature_id)
                    except:
                        line_tmp_dict[key_for_dict] = (gff_line.split('\\t')[3] + ';' + gff_line.split('\\t')[4] + ';' + feature_id )
                    else:
                        continue
                #print(line_tmp_dict)
                outfile.write(protein_id)
                for dict_key in line_tmp_dict:
                    arguments_dict = line_tmp_dict[dict_key].split(';')
                    ids_of_feature = arguments_dict[2]
                    extremities = arguments_dict[0] + ';' + arguments_dict[1]
                    for n in range(3, len(arguments_dict), 3):
                        extremities += (';' + arguments_dict[n] + ';' + arguments_dict[(n+1)])
                        if arguments_dict[(n+2)] in ids_of_feature:
                            continue
                        else:
                            ids_of_feature += (' ' + arguments_dict[(n+2)])
                    outfile.write((' [' + ids_of_feature + ';' + dict_key + ';' + extremities + ']'))
                outfile.write('\n')
            else:
                print(url, file=sys.stderr)

if __name__ == "__main__":
    try:
        input_file = sys.argv[1]
        output_file_name =  sys.argv[2]
    except:
        print('Program usage: text.py <a canonical fasta file that has the identifier of uniprot or uniparc right after the > sign> <the output file:\n that has on ech line the information retrieved by uniprot about the features of the entry \n in the case of uniprot ids the file considered to obtain the informations is the .gff format, in the case of uniparc ids is the xml file the one inquired\n the informationa are stored in a very specific way first field is the proteind identifier then there is a space a [ and the feaature ID (inerpro, PFAM ecc..) then a ; and then the name or decription of the feature in lower cases then a ; and then the number of the stert of the first occurrence of that feature then ; the end of the the first occurrence then ; and so on then a ] then a space then a [ and the same structure described for the feature\n like this\nUPI0006CF02F4 [IPR003343;Bacterial Ig-like, group 2;1044;1109] [IPR003343;Bacterial Ig-like, group 2;435;508;1043;1119] [IPR008964;Invasin/intimin cell-adhesion fragments;429;510]>', file=sys.stderr)
        raise SystemExit
    else:
        list_uniprot, list_uniparc = input_fasta_prepare(input_file)
        #print(list_uniprot)
        #print(list_uniparc)
        uniprot_rest_dealer(list_uniprot, output_file_name)
        uniparc_rest_dealer(list_uniparc, output_file_name)		
