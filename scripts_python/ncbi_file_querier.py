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


def fasta_ncbi_querier (in_fasta, list_ncbi, not_found_verbose=False)
    list_of_accepted_taxid_keywords = ['TaxID=', 'Taxid=', 'taxid=', 'taxID=', 'TaxID =', 'OX=', 'ox=', 'OX =']
    with open(in_fasta, 'r') as fasta:
        for fasta_line in fasta:
            if fasta_line[0]=='>':
                presence_taxid_keyword = False
                print(fasta_line, end='')
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
                        print("species id :", species_id)
                        for n_sublist in range(0, len(list_ncbi)):
                            if species_id in list_ncbi[n_sublist][1]:
                                print(species_id, list_ncbi_db[n_sublist][0])
                        break
                #print(presence_taxid_keyword)
                if presence_taxid_keyword==False:
                    print("the taxid keyworfd used in this line is not present in the possible ones that the script looks for, if you want add it at line 112 in the list called <list_of_accepted_taxid_keywords>\n", "the line inquired is:\n", fasta_line, "\nthe possible taxids searched by this file are:\n", list_of_accepted_taxid_keywords, file=sys.stderr)
                    raise SystemExit
            else:
                continue


if __name__ == "__main__":
    in_database  = None
    try:
        infasta = sys.argv[1]
        in_database = sys.argv[2]
        silencer = sys.argv[3]
    except:
        if in_database != None:
            formatted_categories = ncbi_categories_db_formatter(infasta, in_database)
            fasta_ncbi_querier(infasta, formatted_categories)
        else:
            print('Program usage: text.py <a fasta file that have in the header (the line called header has to start with >) and again in the header the tax id associated to such id/sequence, the format is actually the one of Uniprot with OX=taxid_num or the Uniref format with TaxID=tax_num and variation of them lice all lower case ecc. to see them take a look at the actual code if sections.> < a ncbi downloaded file with three columd tab separsted, with first column is the kingdom one letter identifier, the second is astrand taxid, and the third the species id (both fields are checked for a match) >\n  WARNING \n the script spits out the ids that are not fond, with replicas, if it does not find 849 every time it encounters it along the inputfastalike it will print it, if you want to avoid this behaviour just add a third field with false', file=sys.stderr)
            raise SystemExit
    else:
        formatted_categories = ncbi_categories_db_formatter(infasta, in_database)
        fasta_ncbi_querier(infasta, formatted_categories, silencer)

