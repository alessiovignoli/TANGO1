#!/usr/bin/env python3

import sys

def freq_comp(input_file):
    with open(input_file, 'r') as infile:
        freq_dict = {}
        total_counter = 0
        for line in infile:
            stripped_line = line.rstrip()
            for digit in stripped_line:
                #print(digit)
                total_counter += 1
                if digit in freq_dict:
                    freq_dict[digit] += 1
                else:
                    freq_dict[digit] = 1
        print(freq_dict)
        print(total_counter)
        for digit_key in freq_dict:
            freq_dict[digit_key] = freq_dict[digit_key]/total_counter
        print(freq_dict)


if __name__ == "__main__":
    try:
        one_row_input_file = sys.argv[1]
    except:
        print('Program usage: text.py <a one line file, that have presumably charachters in it. This script will create a dictionary out of every type of digit found in the string and compute frequencies based on the total number of digits in the line.>', file=sys.stderr)
        raise SystemExit
    else:
        freq_comp(one_row_input_file)
