#!/usr/bin/env python3

import argparse
from python_codebase.tabular import TabularFile
from python_codebase.tabular import TabularLine
from python_codebase.file_main import File
from math import log

def get_args():

    "get the arguments when using from the commandline"
    parser = argparse.ArgumentParser(description="This script has three inputs: a bulk GTex data of RNA-seq experiments for example, a python dictionary file saved with pickle containing all the tissue - sample ID relationships, and third a transcript/gene/exon ecc .. ID findable in the bulk data. With all this info the script will compute average expression levels for the given ID for each tissue type found in the dicti, summing all sample from the same Tissue divided by the number encountered for that tissue. Since this script is developed to handle GTex v8 data, it relies on the fact that there is a header line (specifiable) that contains the sample ID present in the Dict, and that the column identified by this sample ID alway contains values associated with it. Basically that column have a meaning related to the samp'le ID.")
    parser.add_argument("-td", "--tissue_dict", type=str, required=True, metavar="FILE", help="The python dictionary containing the relatioshiop between tissues and sample IDs, Tissue names are the keys and the associated samples IDs are the values associated to such key (concatenated in a string separated by commas). This file is expected to be compressed and saved using pickle dump.")
    parser.add_argument("-gt", "--gtex_data", type=str, required=True, metavar="FILE", help="The bulk (usually compressed) data, from which the requested ID has to be extracted. It is assumed that each line contains all the sample values associated to only one gene/transcript ID. Basically every line of this file has a differerent ID per line. It also has to have an header line (specifiable which exactly) that contains the sample-IDs in it as basically column labels/names.")
    parser.add_argument("-id1", "--ID1", type=str, required=True, metavar="STR", help="The first query ID that have to be extracted from the above file.")
    parser.add_argument("-id2", "--ID2", type=str, required=True, metavar="STR", help="The second query ID that have to be extracted from the above file.")
    parser.add_argument("-o", "--out_name", type=str, required=True, metavar="FILE", help='The path to where to save the output.')
    parser.add_argument("-idp", "--id_pos", type=int, required=False, nargs='?', const=0, default=0, metavar="POS", help='The column position (first = 0) that identifies where the query ID is. Default 0.')
    parser.add_argument("-d", "--delimiter", type=str, required=False, nargs='?', const='\t', default='\t', metavar="STR", help='In case it is necessary to specify a different type of separator for the bulk data file. Default tab.')
    parser.add_argument("-hl", "--header_line", type=int, required=False, nargs='?', const=1, default=1, metavar="POS", help='The very impiortant line in bulk data that contains the sample_IDs present as values of the dictionary. From this line the values found in the line of the query Id are correlated (through  the dictionary) to the tissue. first line = 0, default 1 ')

    args = parser.parse_args()
    return args




def main(tissue_dict, gtex_data, ID1, ID2, out_name, id_pos, delimiter, header_line):

    # Load dict
    file_obj = File(tissue_dict)
    loaded_dict = file_obj.PickleLoad()
    
    # Open and extract the query ID lines and place them in a list
    tabular_obj = TabularFile(gtex_data, delimiter)
    query_lines = tabular_obj.GrepLine([ID1, ID2])			# this function returns a list 
    tabline1_obj = TabularLine(query_lines[0], delimiter)
    query_list1 = tabline1_obj.ExtractAllButField(position=id_pos, return_type='list')
    tabline2_obj = TabularLine(query_lines[1], delimiter)
    query_list2 = tabline2_obj.ExtractAllButField(position=id_pos, return_type='list')
    
    # Get header line
    opened_file = tabular_obj.OpenRead(uncompress=True)
    header_list = tabular_obj.ReturnHeader(opened_file, header_lines=header_line)		# returns all the lines up to the one asked
    tabline_obj = TabularLine(header_list[ ( header_line - 1 ) ], delimiter)
    header_elems = tabline_obj.ExtractAllButField(position=id_pos, return_type='list')

    # now that all is loaded the actual look up strategy
    expr_dict = {}
    for i, expr_value in enumerate(query_list1): 
        for tissue, samples in loaded_dict.items():
            if header_elems[i] in samples:
                
                # Deal with divisions by zero when computing log ratio, and log of zero
                if float(query_list2[i]) == 0 or float(expr_value) == 0:
                    continue

                # Update values and counts of existing tissues in final dict
                elif tissue in expr_dict:
                        expr_dict[tissue][0] += log( float(expr_value)/float(query_list2[i]) )
                        expr_dict[tissue][1] += float(expr_value)/float(query_list2[i])
                        expr_dict[tissue][2] += 1

                # To create a new key in expr dict
                else:
                    expr_dict[tissue] = [log( float(expr_value)/float(query_list2[i]) ), float(expr_value)/float(query_list2[i]), 1]
    
    # compute average and write to output
    output_obj = TabularFile(out_name, delimiter)
    opened_out = output_obj.OpenWrite()
    opened_out.write('Tissue' + delimiter + 'average log-ratio' + delimiter + 'average ratio' + delimiter + 'number of samples\n')
    for tissue_name, counts in expr_dict.items():
        opened_out.write( tissue_name + delimiter + str(counts[0]/counts[2]) + delimiter + str(counts[1]/counts[2]) + delimiter +  str(counts[2]) + '\n' )
    



if __name__ == "__main__":
    args = get_args()
   
    main(args.tissue_dict, args.gtex_data, args.ID1, args.ID2, args.out_name, args.id_pos, args.delimiter, args.header_line) 
