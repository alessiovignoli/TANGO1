#!/usr/bin/env python3


import sys

def filterer_OMA(infile, outfile):
    with open(infile, 'r') as OMA_raw, open(outfile, 'w') as OMA_polished:
        #print(OMA_raw, OMA_polished)
        status_check = 0
        first_iteration_counter = 0
        for line in OMA_raw:
            if line[0] == '>':
                status_check = 0
                splitted_line = line.split(' | ')
                #print(splitted_line)
                if len(splitted_line) > 3:
                    ortol_type = splitted_line[2]
                    #print(ortol_type)
                    if ortol_type == '1:1 ortholog':
                        #print(ortol_type)
                        species_name = splitted_line[3]
                        #print(species_name, end='')
                        if '(strain' in species_name or 'sp.' in species_name or ortol_type == 'self':
                            continue
                        #print(line, end='')
                        status_check = 1
                else:
                    species_name = splitted_line[2]
                    if '(strain' in species_name or 'sp.' in species_name:
                        continue
                    status_check = 1
            if status_check == 0:
                continue
            else:
                #print(line, end='')
                if line[0] == '>':
                    if first_iteration_counter == 0:
                        splitted_line = line.split('[')
                        header = splitted_line[0] + '[' + ((splitted_line[1].split(']'))[0]).split(' ')[0] + ' ' + ((splitted_line[1].split(']'))[0]).split(' ')[1] + ']\n'
                        OMA_polished.write(header)
                        first_iteration_counter += 1
                    else:
                        splitted_line = line.split('[')
                        header = splitted_line[0] + '[' + ((splitted_line[1].split(']'))[0]).split(' ')[0] + ' ' + ((splitted_line[1].split(']'))[0]).split(' ')[1] + ']\n'
                        OMA_polished.write('\n')
                        OMA_polished.write(header)
                else:
                    stripped_line = line.rstrip()
                    OMA_polished.write(stripped_line)


if __name__ == "__main__":
    try:
        input_raw_OMA = sys.argv[1]
        output_OMA_filename = sys.argv[2]
    except:
        print('Program usage: text.py <Output of OMA search, like freshly downloaded file (it is recommended to have the appendix _raw before .fasta> <filtered and polished file, same name of the input file but without _raw appendix>', file=sys.stderr)
        raise SystemExit
    else:
        filterer_OMA(input_raw_OMA, output_OMA_filename)
