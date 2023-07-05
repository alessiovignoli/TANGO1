#!/usr/bin/env python3

import argparse
from python_codebase.tabular import TabularFile

def get_args():

    "get the arguments when using from the commandline"

    parser = argparse.ArgumentParser(description="Read GTEx downlaodable sample annotation file, this file contains description and info for all the samples IDs present in the bulk data. This script pèarse such file to build a dictionary that has tissue names as keys and all sample ids belonging to that tissue as arguments of the corresponding tissue key. It themnn savves such file to a compressed (pickle) file.")
    parser.add_argument("-sa", "--sample_annotations", type=str, required=True, metavar="FILE", help='The sample annotation file. In version v8 this is a tab separeted file (tsv).')
    parser.add_argument("-tp", "--tissue_pos", type=int, required=False, nargs='?', const=6, default=6, metavar="POS", help='the column in the sample annotation file containing the Tissue name. In GTEx v8 is column 7 (6) in python notation. This field follows python notation, first column is 0. Default for this fla is 6.')
    parser.add_argument("-sp", "--sample_pos", type=int, required=False, nargs='?', const=0, default=0, metavar="POS", help='the column in the sample annotation file containing the sampleID. In GTEx v8 is column 1 (0) in python notation. This field follows python notation, first column is 0. Default for this fla is 0.')
    parser.add_argument("-d", "--delimiter", type=str, required=False, nargs='?', const='\t', default='\t', metavar="STR", help='In case it is necessary to specify a different type of separator for the annotation file. Default tab.')

    args = parser.parse_args()
    return args




def main(annot_file, tissue_pos, sammpleID_pos, delimiter):

    # Start by extracting the tissue column and obtaining only unique values
    tabular_obj = TabularFile(annot_file, delimiter)
    column_val = tabular_obj.ExtractColumn(tissue_pos, return_type='set')
    
    # Build the dictionary using the above made set of keys and dump it using pickle (compress)
    grouped_dict = tabular_obj.AggregateFromList(sammpleID_pos, tissue_pos, column_val)
    print(grouped_dict)


if __name__ == "__main__":
    args = get_args()
   
    main(args.sample_annotations, args.tissue_pos, args.sample_pos, args.delimiter) 
