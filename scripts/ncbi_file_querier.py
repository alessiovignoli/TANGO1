#!/usr/bin/env python3

import sys

def ncbi_categories_db_formatter (in_ncbi):
    list_ncbi_db = []
    with open(in_ncbi, 'r') as in_db:
        for db_line in in_db:
            kingdom_key = db_line.split()[0]
            species_id_string = db_line.rstrip().split()[2] #+ ' '
            check = False
            number_sublist = 0
            for num_sublist in range(0, len(list_ncbi_db)):
                #print("sublist :", list_ncbi_db[num_sublist], "  index", num_sublist)
                if kingdom_key in list_ncbi_db[num_sublist]:
                    check = True
                    number_sublist = num_sublist
            if check:
                list_ncbi_db[number_sublist][1].append(species_id_string)
            else:
                list_ncbi_db.append([kingdom_key, [species_id_string]])
        #print(list_ncbi_db[0][0], "\t", list_ncbi_db[1][0], "\t", list_ncbi_db[2][0], "\t", list_ncbi_db[3][0], "\t", list_ncbi_db[4][0], "\t", list_ncbi_db[5][0])
        #print(len(list_ncbi_db[0][1]), len(list_ncbi_db[1][1]), len(list_ncbi_db[2][1]), len(list_ncbi_db[3][1]), len(list_ncbi_db[4][1]), len(list_ncbi_db[5][1]))
    return list_ncbi_db


def ncbi_nodes_db_formatter (in_ncbi_nodes):
    mapping_dict_kindom_id = ['B', 'E', 'E', 'V', 'E', 'E', 'E', 'O', 'U',  'V', 'E', 'env']
    list_nodes_db = []
    with open(in_ncbi_nodes, 'r') as db:
        for line in db:
            division_id_key = int(line.split('|')[4].strip())
            kingdom_key = mapping_dict_kindom_id[division_id_key]
            #print(kingdom_key, end=', ')
            species_id = line.split()[0]
            #print(species_id, end=', ')
            check = False
            number_sublist = 0
            for num_sublist in range(0, len(list_nodes_db)):
                if kingdom_key in list_nodes_db[num_sublist]:
                    check = True
                    number_sublist = num_sublist
            if check:
                list_nodes_db[number_sublist][1].append(species_id)
            else:
                list_nodes_db.append([kingdom_key, [species_id]])
    #print(list_nodes_db[0][0], "\t", list_nodes_db[1][0], "\t", list_nodes_db[2][0], "\t", list_nodes_db[3][0], "\t", list_nodes_db[4][0], "\t", list_nodes_db[5][0])
    #print(len(list_nodes_db[0][1]), len(list_nodes_db[1][1]), len(list_nodes_db[2][1]), len(list_nodes_db[3][1]), len(list_nodes_db[4][1]), len(list_nodes_db[5][1]))
    #print(list_nodes_db)
    return list_nodes_db


def fasta_ncbi_querier (in_fasta, list_ncbi, not_found, swith=False):
    list_of_accepted_taxid_keywords = ['TaxID=', 'Taxid=', 'taxid=', 'taxID=', 'TaxID =', 'OX=', 'ox=', 'OX =']
    list_kingdom_nums = ['B', 0, 'A', 0, 'E', 0, 'V', 0, 'U', 0, 'O', 0]
    with open(in_fasta, 'r') as fasta, open(not_found, 'w') as outfile:
        for fasta_line in fasta:
            if fasta_line[0]=='>':
                presence_taxid_keyword = False
                #print(fasta_line, end='')
                for taxid_keyword in list_of_accepted_taxid_keywords:
                    if taxid_keyword in fasta_line:
                        #print("taxid keyword:", taxid_keyword)
                        tmp = fasta_line.split(taxid_keyword)[1].lstrip()
                        species_id = ''
                        presence_taxid_keyword = True
                        for digit in tmp:
                            if digit.isnumeric():
                                species_id += digit
                            else:
                                break
                        #print("species id :", species_id)
                        presence_of_taxid_in_db = False
                        for n_sublist in range(0, len(list_ncbi)):
                            if species_id in list_ncbi[n_sublist][1]:
                                #print(species_id, list_ncbi[n_sublist][0])
                                presence_of_taxid_in_db = True
                                if list_ncbi[n_sublist][0] == 'env':                ## aybe put this in the above def
                                    outfile.write(species_id)                        ## maybe put this in the above def
                                for letter_index in range(0, len(list_kingdom_nums)):
                                    if list_ncbi[n_sublist][0] == list_kingdom_nums[letter_index]:
                                        list_kingdom_nums[(letter_index+1)] += 1
                        if presence_of_taxid_in_db==False:                              ## GET BACK HERE FOR THE PRINT OF NOT FUND
                            #print(species_id)
                            if swith!=False:
                                outfile.write(species_id + '\n')
                            else:
                                outfile.write(fasta_line)
                        break
                #print(presence_taxid_keyword)
                if presence_taxid_keyword==False:
                    print("the taxid keyworfd used in this line is not present in the possible ones that the script looks for, if you want add it at line 112 in the list called <list_of_accepted_taxid_keywords>\n", "the line inquired is:\n", fasta_line, "\nthe possible taxids searched by this file are:\n", list_of_accepted_taxid_keywords, file=sys.stderr)
                    raise SystemExit
            else:
                continue
    print(in_fasta)
    print(list_kingdom_nums)
    total_taxid_found = 0
    for i in range(1, len(list_kingdom_nums), 2):
        total_taxid_found += list_kingdom_nums[i]
    print('total taxid found =  ', total_taxid_found)


if __name__ == "__main__":
    not_found_file_name = None
    try:
        infasta = sys.argv[1]
        in_database = sys.argv[2]
        not_found_file_name = sys.argv[3]
        switch_db = sys.argv[4]
    except:
        if not_found_file_name != None:
            formatted_categories = ncbi_categories_db_formatter(in_database)
            fasta_ncbi_querier(infasta, formatted_categories, not_found_file_name)
        else:
            print('Program usage: text.py <a fasta file that have in the header (the line called header has to start with >) and again in the header the tax id associated to such id/sequence, the format is actually the one of Uniprot with OX=taxid_num or the Uniref format with TaxID=tax_num and variation of them lice all lower case ecc. to see them take a look at the actual code if sections.> < a ncbi downloaded file with three columd tab separsted, with first column is the kingdom one letter identifier, the second is astrand taxid, and the third the species id (both fields are checked for a match) ### ATTENTION  ### or another databese of ncbi called nodes.dmp with a specific architecture, look at kingdom_freq_computer.nf help section or at https://ftp.ncbi.nih.gov/pub/taxonomy/ readme file> <third mandatory field that defines the name of the output file of mot found taxids> <fourth optional field for the switch to the use of nodes.dmp databes command line would be like:\n inputfasta nodes.dmp outputfilename true>\n  WARNING \n the script spits out the ids that are not fond, with replicas, if it does not find 849 every time it encounters it along the inputfastalike it will print it, if you want to avoid this behaviour just add a third field with false', file=sys.stderr)
            raise SystemExit
    else:
        formatted_nodes = ncbi_nodes_db_formatter(in_database)
        fasta_ncbi_querier(infasta, formatted_nodes, not_found_file_name, switch_db)

