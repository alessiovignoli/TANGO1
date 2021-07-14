#!/usr/bin/env python3

# This script check whether a sequence is in common between the multifasta output of OMA and the 
# output multifasta of Metaphor. Better if the filenames are given in this order.
# it takes in input two fasta files, one from OMA, one from Metaphor and prints to screen the identical sequences 
# by string exact matching, with their respecvtive header.
# it also reports if there are sequences from same species but different in exact string matching.
# It uses the set and intersect operations so is not stable, sequences are gonna come out randomized.

import sys

def intersectioner(file1, file2):
    #print('bubba')
    setf1 = set()
    setf2 = set()
    dictf1 = {}
    dictf2 = {}
    #print(type(set1))
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        #print(f1, f2)
        #print(m)
        list_linesf1 = f1.readlines()
        list_linesf2 = f2.readlines()
        n = 0
        seq_buffer = 'bubba'
        for m in range(0,(len(list_linesf1))):
            #print(line, end='')
            if list_linesf1[m][0] == '>':
                #print('bubba')
                setf1.add(seq_buffer)
                dictf1[seq_buffer] = list_linesf1[n].rstrip()
                n = m
                seq_buffer = ''
            elif m == (len(list_linesf1) -1):
                #print('bubba')
                #print(list_linesf1[m], end='')
                setf1.add(seq_buffer)
                dictf1[seq_buffer] = list_linesf1[n].rstrip()
            else:
                seq_buffer += list_linesf1[m].rstrip()
        dictf1.pop('bubba', None)
        setf1.remove('bubba')
        list_linesf1 = ''
        n = 0
        seq_buffer = 'bubba'
        for m in range(0,(len(list_linesf2))):
            #print(line, end='')
            if list_linesf2[m][0] == '>':
                #print('bubba')
                setf2.add(seq_buffer)
                dictf2[seq_buffer] = list_linesf2[n].rstrip()
                n = m
                seq_buffer = ''
            elif m == (len(list_linesf2) -1):
                #print('bubba')
                #print(list_linesf2[m], end='')
                setf2.add(seq_buffer)
                dictf2[seq_buffer] = list_linesf2[n].rstrip()
            else:
                seq_buffer += list_linesf2[m].rstrip()
        dictf2.pop('bubba', None)
        setf2.remove('bubba')
    #print(len(setf1))
    #print(len(dictf1))
    #print(len(setf2))
    #print(len(dictf2))
    intersection = setf1.intersection(setf2)
    #print(len(intersection))
    #print('\n##########    Sequence identical by exact string matching  ###########\n')
    for elem in intersection:
        #if elem in dictf1:
            #print(dictf1[elem])
        if elem in dictf2:
            oma_seq_id = ((dictf1[elem]).split('_')[1]) #[1:]
            oma_ortol_type = ((dictf1[elem]).split('_')[2])#split()[0])
            metaphor_seq_id = ((dictf2[elem]).split('_')[0])[1:]
            uniprot_seq_id = (dictf1[elem]).split('_')[3] #'Na'
            species = ((dictf2[elem]).split('_')[4]) + '_' + ((dictf2[elem]).split('_')[5]) #.split(']')[0]
            #if '|' in dictf2[elem]:
                #print(dictf2[elem])
                #metaphor_seq_id = ((dictf2[elem]).split('|')[0]).split('_')[0]
                #uniprot_seq_id = (dictf2[elem]).split('|')[1]
            #print(f1_seq_id, oma_ortol_type)
            #print(metaphor_seq_id + ' ' +  oma_seq_id +  oma_ortol_type + ' '+ uniprot_seq_id + ' ' + species)
            print('>' + species, oma_ortol_type, uniprot_seq_id, metaphor_seq_id, oma_seq_id)
            #print(dictf2[elem])
        print(elem)
    '''print('\n##########    Sequence that are from the same organism, need for manual curation  ###########\n')
    for sequence, header in dictf1.items():
        #print(header)
        #print(sequence)
        #species_name = (header.split('[')[1])[:-1]
        species_name = header.split(' ')[4] + ' ' + header.split(' ')[5]
        #print(species_name)
        for key, value in dictf2.items():
            #print(value)
            if species_name in value and sequence != key:
                print(header)
                print(value)
                print('>>> ', file1, sequence)
                print('>>> ', file2, ' ', key)
'''
if __name__ == "__main__":
    try:
        fasta_file1 = sys.argv[1]
        fasta_file2 = sys.argv[2]
    except:
        print('Program usage: text.py <fasta_file1.txt> <fasta_file2.txt>', file=sys.stderr)
        raise SystemExit
    else:
        intersectioner(fasta_file1, fasta_file2)
            
