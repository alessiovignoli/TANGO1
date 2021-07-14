#!/usr/bin/env python3


import sys

def non_rng_selecter(infile, tot_num, desired_num, outfile):
    if tot_num < desired_num:
        print("the number specified for the total num of entries =", tot_num, file=sys.stderr)
        print("the number specified for num of entries to select =", desired_num, file=sys.stderr)
        print("WARNING the script can not select more entries then what is said to be in the input file", file=sys.stderr)
        raise SystemExit
    with open(infile, 'r') as infasta, open(outfile, 'w') as outfasta:
        #print(infasta, outfasta)
        #print("total num of entries =", tot_num, "   num of entries to select =", desired_num)
        ratio = tot_num // desired_num
        num_iterations_with_add_ratio = tot_num - (ratio * desired_num)
        #print("ratio =", ratio, "   num iterations with ratio+1 =", num_iterations_with_add_ratio) 
        i = 0
        n = 0
        m = 1
        match = False
        if num_iterations_with_add_ratio == 0:
            n = ratio
        else:
            n = (ratio + 1) 
        for line in infasta:
            if line[0] == '>':
                i += 1
                match = False
                if m <= desired_num:
                    if i == n and m < num_iterations_with_add_ratio:
                        #print(" i =", i, "   n=", n, "   m=", m, "   num_iterations_with_add_ratio=", num_iterations_with_add_ratio)
                        #print(line, end='')
                        n += (ratio+1)
                        m += 1
                        match = True
                        outfasta.write(line)
                    elif i == n and m >= num_iterations_with_add_ratio:
                        #print(" i =", i, "   n=", n, "   m=", m, "   num_iterations_with_add_ratio=", num_iterations_with_add_ratio)
                        #print(line, end='')
                        n += ratio
                        m += 1
                        match = True
                        outfasta.write(line)
                else:
                    break
            else:
                if match:
                    #print(line, end='')
                    outfasta.write(line)



if __name__ == "__main__":
    try:
        input_fasta = sys.argv[1]
        total_number_entries = int(sys.argv[2])
        desired_number_entries = int(sys.argv[3])
        output_fasta = sys.argv[4]
    except:
        print('Program usage: text.py <input multifatsa or multifasta like> <total number of entries (headers with ">" at  start) in the input file> <number of the ammount of entries the output file must contain (the entries are chosen one every 3,4,7,42 ecc.) depending on the ratio of the two given numbers (a little more complicated than that)> <output filename>', file=sys.stderr)
        raise SystemExit
    else:
        non_rng_selecter(input_fasta, total_number_entries, desired_number_entries, output_fasta)
