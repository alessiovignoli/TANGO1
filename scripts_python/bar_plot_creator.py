#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

def from_input_to_numpy_array(txt_file, max_iter=1000000, x_field=None, y_field=None):
    with open(txt_file, 'r') as intxt:
        #print(intxt, x_field, y_field)
        i = 0
        x_array = np.empty([max_iter], dtype='U100')            # setting maximum length of dtring type to add to array otherwise str uses only first letter
        y_array = np.empty([max_iter], dtype=int)
        for txtline in intxt:
            if i == max_iter:
                #print('if section')
                return x_array, y_array
            elif txtline[0] == '#':
                continue
            else:
                #print('else section')
                tmp = txtline.split()[x_field]
                for parola in txtline.split()[(x_field + 1):]:
                    tmp += ' '
                    tmp += parola
                x_array[i] = tmp                                    # this will take the sting argument from field x till end of line
                y_array[i] = int(txtline.split()[y_field])
            i += 1
        #print(i)
        return x_array[0:i], y_array[0:i]


if __name__ == "__main__":
    try:
        intxtfile1 = sys.argv[1]
        #intxtfile2 = sys.argv[2]
        max_number_iterations = int(sys.argv[2])
    except:
        raise SystemExit
    else:
        x_ar_1, y_ar_1 = from_input_to_numpy_array(intxtfile1,  max_number_iterations, 1, 0)
        #print(x_ar_1)
        #print(y_ar_1)
        plt.barh(x_ar_1, y_ar_1)
        plt.show()
