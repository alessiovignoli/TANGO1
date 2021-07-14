#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

# sometimes x_field and y_field have to be changed
def from_txt_to_numpy_array(txt_file, max_iter=1000, x_field=0, y_field=1):
    with open(txt_file, 'r') as intxt:
        #print(intxt, x_field, y_field)
        i = 0
        x_array = np.empty([1, max_iter], dtype=float)
        y_array = np.empty([1, max_iter], dtype=float)
        for txtline in intxt:
            if i == max_iter:
                #print('if section')
                return x_array, y_array
            else:
                #print('else section')
                x_array[0][i] = float(txtline.split()[x_field])
                y_array[0][i] = float(txtline.split()[y_field])
            i += 1
        #print(i)

if __name__ == "__main__":
    try:
        intxtfile1 = sys.argv[1]
        intxtfile2 = sys.argv[2]
        max_number_iterations = int(sys.argv[3])
    except:
        raise SystemExit
    else:
        x_ar_1, y_ar_1 = from_txt_to_numpy_array(intxtfile1, max_number_iterations, 1, 0)
        x_ar_2, y_ar_2 = from_txt_to_numpy_array(intxtfile2, max_number_iterations)         # in this case x and y are changed like comment above
        #print(x_ar_1)
        #print(y_ar_1)
        fig, ax = plt.subplots()
        ax.scatter(x_ar_1[0][0], 0.69, c='#d62728', s=11,  alpha=0.22, label='special helix segments with contribution from normal helix')
        ax.scatter(x_ar_1, y_ar_1, c='#d62728', s=11, alpha=0.26)
        ax.scatter(x_ar_2[0][0], 0.69, c='#1f77b4', s=11,  alpha=0.22, label='normal helix segments with contribution from special helix')
        ax.scatter(x_ar_2, y_ar_2, c='#1f77b4', s=11, alpha=0.26)
        #ax.set(xticks=[0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0])
        #ax.set(xticks=[0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        #ax.set(yticks=[0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        ax.set(xticks=[0, 0.5, 1.0])
        ax.set(yticks=[0, 0.5, 1.0])
        ax.set_xlabel('plp average normal helix', fontsize=16)
        ax.set_ylabel('plp average special helix', fontsize=16)
        ax.set_title("Comparison between 20000 normal segments and 20000 special", fontsize=18) #and 20000 normal helix segments", fontsize=18)
        #print(a, type(a))
        line1 = plt.plot([0.0, 0.25, 0.5], [0.5, 0.25, 0.0])
        plt.setp(line1, 'color', 'g', 'linewidth', 2.0, alpha=0.6, label='normal helix plpaverage plus special helix plpaverage = 0.5')
        line2 = plt.plot([0.0, 0.25, 0.75], [0.75, 0.5, 0.0])
        plt.setp(line2, 'color', '#a4eb17', 'linewidth', 2.0, alpha=0.9, label='normal helix plpaverage plus special helix plpaverage = 0.75')
        ax.grid()
        ax.legend(fontsize=16)
        plt.show()
