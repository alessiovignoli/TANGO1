#!/usr/bin/env python3

import sys

def fasta_header_tax_id_ncbi_query (in_fasta, in_ncbi, not_found_verbose=False):
    print(in_fasta, in_ncbi, not_found_verbose)

"""    count_dict = {}
    with open(in_ids, 'r') as inheader, open(in_fasta, 'r') as infasta:
        tot_count= 0
        species_id = 0
        tmp = None
        for line in inheader:
            seq_id = line.split()[0]
            for line2 in infasta:
                tmp = seq_id
                if line2[0] == '>' and seq_id in line2:
                    species_id = None
                    if 'TaxID=' in line2:
                        species_id = line2.split('TaxID=')[1].split(' ')[0]
                    elif 'Taxid=' in line2:
                        species_id = line2.split('Taxid=')[1].split(' ')[0]
                    elif 'taxid=' in line2:
                        species_id = line2.split('taxid=')[1].split(' ')[0]
                    elif 'taxID=' in line2:
                        species_id = line2.split('taxID=')[1].split(' ')[0]
                    elif 'TaxID =' in line2:
                        species_id = line2.split('TaxID =')[1].split(' ')[0]
                    elif 'OX=' in line2:
                        species_id = line2.split('OX=')[1].split(' ')[0]
                    elif 'ox=' in line2:    
                        species_id = line2.split('ox=')[1].split(' ')[0]
                    elif 'OX =' in line2:
                        species_id = line2.split('OX =')[1].split(' ')[0]
                    else:
                        print('there is no tax id field identifier tried in this script in the header, please check if there is TaxID= or OX= in the headers, this is the problematic line', line2,  file=sys.stderr)
                        raise SystemExit
                    print(seq_id, species_id)
                    break
        if species_id == 0:
            print('this is id has not been found in the fasta file :', tmp, file=sys.stderr)
"""

if __name__ == "__main__":
    in_database  = None
    try:
        infasta = sys.argv[1]
	    in_database = sys.argv[2]
	    silencer = sys.argv[3]
    except:
        if in_database != None:
            fasta_header_tax_id_ncbi_query(infasta, in_database)
        else:
            print('Program usage: text.py <a fasta file that have in the header (the line called header has to start with >) and again in the header the tax id associated to such id/sequence, the format is actually the one of Uniprot with OX=taxid_num or the Uniref format with TaxID=tax_num and variation of them lice all lower case ecc. to see them take a look at the actual code if sections.> < a ncbi downloaded file with three columd tab separsted, with first column is the kingdom one letter identifier, the second is astrand taxid, and the third the species id (both fields are checked for a match) >\n  WARNING \n the script spits out the ids that are not fond, with replicas, if it does not find 849 every time it encounters it along the inputfastalike it will print it, if you want to avoid this behaviour just add a third field with false', file=sys.stderr)
            raise SystemExit
    else:
        fasta_header_tax_id_ncbi_query(infasta, in_database, silencer)
