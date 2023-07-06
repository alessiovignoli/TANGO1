#!/usr/bin/env python3

import argparse
from python_codebase.tabular import TabularFile
from python_codebase.file_main import File

def get_args():

    "get the arguments when using from the commandline"
    parser = argparse.ArgumentParser(description="This script has three inputs: a bulk GTex data of RNA-seq experiments for example, a python dictionary file saved with pickle containing all the tissue - sample ID relationships, and third a transcript/gene/exon ecc .. ID findable in the bulk data. With all this info the script will compute average expression levels for the given ID for each tissue type found in the dicti, summing all sample from the same Tissue divided by the number encountered for that tissue. Since this script is developed to handle GTex v8 data, it relies on the fact that there is a header line (specifiable) that contains the sample ID present in the Dict, and that the column identified by this sample ID alway contains values associated with it. Basically that column have a meaning related to the samp'le ID.")
    parser.add_argument("-td", "--tissue_dict", type=str, required=True, metavar="FILE", help="The python dictionary containing the relatioshiop between tissues and sample IDs, Tissue names are the keys and the associated samples IDs are the values associated to such key (concatenated in a string separated by commas). This file is expected to be compressed and saved using pickle dump.")
    parser.add_argument("-gt", "--gtex_data", type=str, required=True, metavar="FILE", help="The bulk (usually compressed) data, from which the requested ID has to be extracted. It is assumed that each line contains all the sample values associated to only one gene/transcript ID. Basically every line of this file has a differerent ID per line. It also has to have an header line (specifiable which exactly) that contains the sample-IDs in it as basically column labels/names.")
    parser.add_argument("-id", "--ID", type=str, required=True, metavar="STR", help="The query ID that have to be extracted from the above file.")
    parser.add_argument("-o", "--out_name", type=str, required=True, metavar="FILE", help='The path to where to save the output.')
    parser.add_argument("-idp", "--id_pos", type=int, required=False, nargs='?', const=0, default=0, metavar="POS", help='The column position (first = 0) that identifies where the query ID is. Default 0.')
    parser.add_argument("-d", "--delimiter", type=str, required=False, nargs='?', const='\t', default='\t', metavar="STR", help='In case it is necessary to specify a different type of separator for the bulk data file. Default tab.')
    parser.add_argument("-hl", "--header_line", type=int, required=False, nargs='?', const=1, default=1, metavar="POS", help='The very impiortant line in bulk data that contains the sample_IDs present as values of the dictionary. From this line the values found in the line of the query Id are correlated (through  the dictionary) to the tissue. first line = 0, default 1 ')

    args = parser.parse_args()
    return args




def main(tissue_dict, gtex_data, ID, out_name, id_pos, delimiter, header_line):

    #print(tissue_dict, gtex_data, ID, out_name, id_pos, delimiter, header_line)

    # Load dict
    file_obj = File(tissue_dict)
    loaded_dict = file_obj.PickleLoad()

    # Open and extract the query ID line
    tabular_obj = TabularFile(gtex_data, delimiter)
    query_line = tabular_obj.GrepLine(ID)
    print(query_line)
    #file_obj = File(out_name)
    



if __name__ == "__main__":
    args = get_args()
   
    main(args.tissue_dict, args.gtex_data, args.ID, args.out_name, args.id_pos, args.delimiter, args.header_line) 
