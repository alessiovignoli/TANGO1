#!/usr/bin/env python3

# This script is intended to create a distribution of integers (bins) referring to the legths of a specific portion of a set of proteins
# it can work with just one keyword, and that portion reffered from the keyword when found in a protein will be added to a bin of
# legth distribution
# alternatively the script can operate with three keyword to look for a pattern, the best would be that the first keyword 
# refers to the region right before the portion we want to look at, the second keyword would refer to the portion of interest of the protein,
# the third keyword would refer to the region right after the interested part.
# the script has a last option rewferring to the occurence of the interst portion, the first occurence of the pattern or keyword (1) 
# or the second (2), or all the occurrences in a protein (all).

import sys

def distributer(in_phobius_file, occurrence_counter, keyword1, keyword2=None, keyword3=None):
    #print(in_phobius_file, occurrence_counter, keyword1, keyword2, keyword3)
    frequences = {}
    with open(in_phobius_file, 'r') as infile:
        prot_check = False
        occ_instance1 = 0                      # to know on which occurrence of the pattern i am
        target_line = ''
        upper_surraunding_check = False
        target_check = False
        for line in infile:
            #print(line.rstrip(), 'prot', prot_check, 'upper', upper_surraunding_check, 'target', target_check, 'occurr', occ_instance1, keyword3)
            #print(frequences)
            if occurrence_counter != 'all':
                occur_counter = int(occurrence_counter)
                if keyword3 == None:
                    #print(line, end='')
                    if keyword1 in line and prot_check == False:
                        #print(line, end='')
                        occ_instance1 += 1
                    if occ_instance1 == occur_counter:
                        #print(line, end='')
                        occ_instance1 = 0
                        prot_check = True
                        start = int(line[15:21].strip())
                        end = int(line[23:28].strip())
                        length = end - start + 1
                        if length in frequences:
                            frequences[length] += 1
                        else:
                            frequences[length] = 1
                        #print(frequences)
                    elif 'ID' in line:
                        #print(line, end='')
                        occ_instance1 = 0
                        prot_check = False
                else:
                    if 'ID' in line:
                        #print(line, end='')
                        occ_instance1 = 0
                        prot_check = False
                        upper_surraunding_check = False
                        target_check = False
                    else:
                        #print(line, end='')
                        if keyword3 in line and prot_check == False and upper_surraunding_check == True and target_check == True:
                            #print(line, end='')
                            occ_instance1 += 1
                            if occ_instance1 == occur_counter:
                                #print(target_line, end='')
                                occ_instance1 = 0
                                prot_check = True
                                start = int(target_line[15:21].strip())
                                end = int(target_line[23:28].strip())
                                length = end - start + 1
                                if length in frequences:
                                    frequences[length] += 1
                                else:
                                    frequences[length] = 1
                        elif keyword1 in line and prot_check == False:
                            #print(line, end='')
                            upper_surraunding_check = True
                        elif keyword2 in line and prot_check == False and upper_surraunding_check == True:
                            #print(line, end='')
                            target_line = line
                            target_check = True
                            #continue
                        else:
                            #print(line, end='')
                            upper_surraunding_check = False
                            target_check = False
            #print(line.rstrip(), 'prot', prot_check, 'upper', upper_surraunding_check, 'target', target_check, 'occurr', occ_instance1)

            else:
                if keyword3 == None:
                    if keyword1 in line:
                        #print(line, end='')
                        start = int(line[15:21].strip())
                        end = int(line[23:28].strip())
                        length = end - start + 1
                        if length in frequences:
                            frequences[length] += 1
                        else:
                            frequences[length] = 1
                else:
                    if keyword3 in line and prot_check == False and upper_surraunding_check == True and target_check == True:
                            #print(target_line, end='')
                            start = int(target_line[15:21].strip())
                            end = int(target_line[23:28].strip())
                            length = end - start + 1
                            if length in frequences:
                                frequences[length] += 1
                            else:
                                frequences[length] = 1
                    elif keyword1 in line and prot_check == False:
                        #print(line, end='')
                        upper_surraunding_check = True
                    elif keyword2 in line and prot_check == False and upper_surraunding_check == True:
                        #print(line, end='')
                        target_line = line
                        target_check = True
                        #continue
                    else:
                        #print(line, end='')
                        upper_surraunding_check = False
                        target_check = False

    print(frequences)    


if __name__ == "__main__":
    multi_phobius_output_file = None
    occurrence_counter_specifier = None
    first_keyword = None
    try:
        multi_phobius_output_file = sys.argv[1]
        occurrence_counter_specifier = sys.argv[2]
        first_keyword = sys.argv[3]
        second_keyword = sys.argv[4]
        third_keyword = sys.argv[5]
    except:
        if multi_phobius_output_file is not None and occurrence_counter_specifier is not None and first_keyword is not None:
            distributer(multi_phobius_output_file, occurrence_counter_specifier, first_keyword)
        else:
            print('Program usage: text.py <input must be a collection of phobius standard (long) outputs in  a sigle file, just a redirection of the stdout to a file on a phobius command> <the occurrence counter (1,2...integer) or (all) for all the occurrences> <a specified capslock keyword used from phobius to describe a specific partion of the protein> (optional)==> <a specified capslock keyword used from phobius to describe a specific partion of the protein> <a specified capslock keyword used from phobius to describe a specific partion of the protein>', file=sys.stderr)
            raise SystemExit
    else:
        distributer(multi_phobius_output_file, occurrence_counter_specifier, first_keyword, second_keyword, third_keyword)

