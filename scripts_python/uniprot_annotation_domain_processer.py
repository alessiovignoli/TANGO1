#!/usr/bin/env python3

import sys

def all_annotation_freq_computer(inputfile, outputfile):
    with open(inputfile, 'r') as in_file, open(outputfile, 'w') as out_file:
        list_of_descript = []
        for line in in_file:
            #print(line.split('[')[0])
            for annotation in line.split(' [')[1:]:
                try:                                                #taking care of open sqare brachets in annotation names
                    description = annotation.split(';')[1]
                    #print('description=', description)
                except Exception:
                    continue
                else:
                    if description in list_of_descript:
                        for n in range(0, len(list_of_descript), 2):
                            if description == list_of_descript[n]:      #try and remove inside spaces lines and dashes
                                list_of_descript[(n+1)] += 1
                                break
                    else:
                        list_of_descript.append(description)
                        list_of_descript.append(1)
        #print('\nlist =\n', list_of_descript)
        #print('unique annotations found =', (len(list_of_descript)/2))
        total_count = 0
        for i in range(0, len(list_of_descript), 2):
            #print(list_of_descript[i] + '\t' + str(list_of_descript[(i+1)]))
            out_file.write(str(list_of_descript[(i+1)]) + ' ' + list_of_descript[i] +'\n')
            total_count += list_of_descript[(i+1)]
        print('unique annotations found =', (len(list_of_descript)/2))
        print('total annotation found =', total_count)


if __name__ == "__main__":
    try:
        input_file = sys.argv[1]
        output_file_name =  sys.argv[2]
    except Exception:
        print('Program usage: text.py a domain_info file that has on each line a protein annotation in a special format,\nlike:\nUPI00068AA065 [IPR036514;sgnh hydrolase superfamily;471;634] [IPR002656;acyltransferase 3 domain;5;321] [IPR043968;sgnh domain;400;626]\nwhere the first field is the protein id and then there are as many lists as there have been found annotations for this entry in uniprot, there is the keyword associted with each feature a word descriptive feature and the extremities of the feature, there can be more of them\nto generate this file the script   uniprot_rest_query.py     look into thoat for more details> < the outputfilename that will have for each unique (exact string match more or less)  name of domain bor annotation  the times that such string was found in the file , basically the scripts run through each line (protein annotations) and updates a counter for each time a unique string is found, the final number will be the times the annotation "gfp" has been found in the file > ', file=sys.stderr)
        raise SystemExit
    else:
        all_annotation_freq_computer(input_file, output_file_name)
