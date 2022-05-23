#!/usr/bin/env python3

import sys

def species_grouper(input_file, presence_header=False, print_check=False, out_file=False):
    with open(input_file) as infile:
        if presence_header == 'true' or presence_header == 'True' or presence_header == True:
            infile.readline()
            #print('header present')
        dict_species = {}
        dict_of_counters = {}
        counter_of_lines = 0
        #print(type(dict_species))
        for line in infile:
            counter_of_lines += 1
            #print(line, end='')
            if ',' in line:
                splitted_line = (line.rstrip()).split(',')
                if len(splitted_line) >= 3:
                    print('there are more than one comma separating this line:  ', counter_of_lines, file=sys.stderr)
                    raise SystemExit
                else:
                    if (splitted_line[1]) in dict_species:
                        dict_species[(splitted_line[1])] += (splitted_line[0] + ' ')
                        dict_of_counters[(splitted_line[1])] += 1
                    else:
                        dict_species[(splitted_line[1])] = (splitted_line[0] + ' ')
                        dict_of_counters[(splitted_line[1])] = 1
            elif ':' in line:
                splitted_line = (line.rstrip()).split(':')
                if len(splitted_line) >= 3:
                    print('there are more than one comma separating this line:  ', counter_of_lines, file=sys.stderr)
                    raise SystemExit
                else:
                    if (splitted_line[1]) in dict_species:
                        dict_species[(splitted_line[1])] += (splitted_line[0] + ' ')
                        dict_of_counters[(splitted_line[1])] += 1
                    else:
                        dict_species[(splitted_line[1])] = (splitted_line[0] + ' ')
                        dict_of_counters[(splitted_line[1])] = 1
            elif ' ' in line:
                splitted_line = (line.rstrip()).split(' ')
                if len(splitted_line) >= 3:
                    print('there are more than one comma separating this line:  ', counter_of_lines, file=sys.stderr)
                    raise SystemExit
                else:
                    if (splitted_line[1]) in dict_species:
                        dict_species[(splitted_line[1])] += (splitted_line[0] + ' ')
                        dict_of_counters[(splitted_line[1])] += 1
                    else:
                        dict_species[(splitted_line[1])] = (splitted_line[0] + ' ')
                        dict_of_counters[(splitted_line[1])] = 1
        if print_check == 'true' or print_check == 'True' or print_check == True:
            sorted_dict = sorted(dict_of_counters.items(), reverse=True, key=lambda item: item[1])
            #print(sorted_dict)
            #print(type(sorted_dict))
            if out_file != False:
                with open(out_file, 'w') as outfile:    
                    for elem in sorted_dict:
                        to_be_written = str(elem[0]) + ' : ' + str(elem[1]) + '\n'
                        outfile.write(to_be_written)
            else:
                for elem in sorted_dict:
                    print(elem[0], ':', elem[1])
                print()
                print(len(dict_of_counters))
            print(dict_species)
        return dict_species


if __name__ == "__main__":
    try:
        ids_plus_taxids = sys.argv[1]
        header = sys.argv[2]
        print_details = sys.argv[3]
    except:
        print('Program usage: text.py  < a file containing lines with two columns separated by either tab, comma, column or space. In which the first column is the id of the protein/gene and second column is the organism TaxID or the organism name separated by _ sign> < optionalfield, a boolean for the presence of a one line header> < optional field, a boolean for printing details about the number of species and ammount of ids per specie, true means that it gets print> < optional field, mandatory only if the print deteils must be written to a file>', file=sys.stderr)
        raise SystemExit
    else:
        output_file = False
        if print_details == 'true' or print_details == 'True':
            output_file = sys.argv[4]
        species_grouper(ids_plus_taxids, header, print_details, output_file)
