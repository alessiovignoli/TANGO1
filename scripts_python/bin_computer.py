#!/usr/bin/env python3



import sys
from decimal import Decimal

def creator_of_bins(in_file, column=0, num_bins=10, range_nums="0,1"):
    #print(in_file, column, num_bins, range_nums)
    with open(in_file,'r') as infile:
        dict_bins = {}
        interval_key = (Decimal(range_nums.split(',')[1]) - Decimal(range_nums.split(',')[0])) / num_bins
        #print(interval_key, type(interval_key))
        tmp = Decimal(range_nums.split(',')[0])
        for i in range(1, (num_bins+1)):
            tmp += Decimal(interval_key)
            dict_bins[tmp] = 0
        #print(dict_bins)
        for line in infile:
            try:
                desired_value = Decimal((line.strip()).split()[column])
            except:
                print("unable to acquire column", column,"\nat line", line, "check whether is space/tab separated, or if the column is the desired one, remember first column is 0, or even worst the Decimal function o that number", file=sys.stderr)
            else:
                #print(desired_value)
                prev_key = interval_key
                #print("prev_key", prev_key)
                #if desired_value >= 0.75:
                #    print(line.strip().split()[1], line.strip().split()[2])
                for decim_key in dict_bins:
                    #print(decim_key, end=' ')
                    if desired_value > decim_key:
                        prev_key = dict_bins
                        continue
                    else:
                        #print(desired_value, decim_key)
                        dict_bins[decim_key] += 1
                        #print(dict_bins)
                        break
        print(dict_bins)


if __name__ == "__main__":
    try:
        input_file = sys.argv[1]
        column_keyword = int(sys.argv[2])
        number_of_bins = int(sys.argv[3])
        range_of_nums = sys.argv[4]
    except:
        if input_file != None:
            creator_of_bins(input_file)
        else:
            print('Program usage: text.py < a space or tab separated column file, where the selected column holds numbers> < optional arg, the column to be selected as integer, default equal to 0, that means first column, second column is 1 ecc.> < optional arg, number of bin to be created as integer, default is 10, otherwise automatically defines the intervals based on the number of bins and on the next parameter> < optional arg, range of the numbers, mining minimun value, maximum value, default is 0,1>', file=sys.stderr)
            raise SystemExit
    else:
        creator_of_bins(input_file, column_keyword, number_of_bins, range_of_nums)
