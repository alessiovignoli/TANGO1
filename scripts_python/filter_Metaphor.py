#!/usr/bin/env python3


import sys


def filterer_Met(infile, outfile):
    with open(infile, 'r') as Met_raw, open(outfile, 'w') as Met_polished:
        first_iteration_counter = 0
        flter_out = 0
        for line in Met_raw:
            if line[0] == '>':
                flter_out = 0
                header_elements = line.split('|')
                len_header_elements = len(header_elements)
                #print(header_elements)
                #print(len_header_elements)
                species_name = (line.split('[')[1]).split(']')[0]
                #print(species_name)
                splitted_speciename = species_name.split(' ')
                if 'strain' in species_name or 'sp.' in species_name:
                        flter_out = 1
                else:
                    if '(' in splitted_speciename[1]:
                        if len_header_elements == 1:
                            new_header = line.split('[')[0] + '| Na | [' + splitted_speciename[0] + ' ' + splitted_speciename[2] + ']\n'
                            if first_iteration_counter == 0:
                                Met_polished.write(new_header)
                                first_iteration_counter += 1
                            else:
                                #print(header_elements)
                                Met_polished.write('\n')
                                Met_polished.write(new_header)
                                #print(line.split('[')[0] + 'Na [' + line.split('[')[1], end='')
                                #continue
                        else:
                            new_header = header_elements[0] + ' | ' + (header_elements[(len_header_elements-1)]).split('[')[0] + '| [' + splitted_speciename[0] + ' ' + splitted_speciename[2] + ']\n'
                            if first_iteration_counter == 0:
                                Met_polished.write(new_header)
                                first_iteration_counter += 1
                            else:
                                #print(header_elements[0], header_elements[(len_header_elements-1)], end='')
                                Met_polished.write('\n')
                                Met_polished.write(new_header)
                    else:
                        if len_header_elements == 1:
                            new_header = line.split('[')[0] + '| Na | [' + splitted_speciename[0] + ' ' + splitted_speciename[1] + ']\n'
                            if first_iteration_counter == 0:
                                Met_polished.write(new_header)
                                first_iteration_counter += 1
                            else:
                                #print(header_elements)
                                Met_polished.write('\n')
                                Met_polished.write(new_header)
                                #print(line.split('[')[0] + 'Na [' + line.split('[')[1], end='')
                                #continue
                        else:
                            new_header = header_elements[0] + ' | ' + (header_elements[(len_header_elements-1)]).split('[')[0] + '| [' + splitted_speciename[0] + ' ' + splitted_speciename[1] + ']\n'
                            if first_iteration_counter == 0:
                                Met_polished.write(new_header)
                                first_iteration_counter += 1
                            else:
                                #print(header_elements[0], header_elements[(len_header_elements-1)], end='')
                                Met_polished.write('\n')
                                Met_polished.write(new_header)
                #continue
            #print((line.rstrip()), end='')
            else:
                if flter_out == 1:
                    continue
                else:
                    Met_polished.write(line.rstrip())

if __name__ == "__main__":
    try:
        input_raw_Met = sys.argv[1]
        output_Met_filename = sys.argv[2]
    except:
        print('Program usage: text.py <Output of MetaPhor search, like freshly downloaded file (it is recommended to have the appendix _raw before .fasta> <filtered and polished file, same name of the input file but without _raw appendix>', file=sys.stderr)
        raise SystemExit
    else:
        filterer_Met(input_raw_Met, output_Met_filename)
