#!/usr/bin/env python3

import sys

def all_annotation_freq_computer(inputfile, outputfile):
    with open(inputfile, 'r') as in_file, open(outputfile, 'w') as out_file:
        list_of_descript = []
        for line in in_file:
            #print(line.split('[')[0])
            for annotation in line.split('[')[1:]:
                description = annotation.split(';')[1]
                #print('description=', description)
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
    except:
        print('Program usage: text.py <', file=sys.stderr)
        raise SystemExit
    else:
        all_annotation_freq_computer(input_file, output_file_name)
