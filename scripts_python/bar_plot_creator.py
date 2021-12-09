#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from matplotlib.colors import Normalize
#from numpy.ma import masked_array

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

def color_prepare(x_str1, x_str2, x_str3):
    collect_list = [x_str1, x_str2, x_str3]
    dict_of_possible_str = {}
    final_color_list = []
    i = 0
    for lista in collect_list:
        specific_list = []
        for elem in lista:
            if elem in dict_of_possible_str:
                specific_list.append(dict_of_possible_str[elem])
            else:
                dict_of_possible_str[elem] = i
                specific_list.append(i)
                i += 1
        final_color_list.append(specific_list)
    return final_color_list[0], final_color_list[1], final_color_list[2], len(dict_of_possible_str)


if __name__ == "__main__":
    try:
        intxtfile1 = sys.argv[1]
        intxtfile2 = sys.argv[2]
        intxtfile3 = sys.argv[3]
        max_number_iterations = int(sys.argv[4])
    except:
        raise SystemExit
    else:
        x_ar_1, y_ar_1 = from_input_to_numpy_array(intxtfile1,  max_number_iterations, 1, 0)
        print(x_ar_1)
        print(y_ar_1)
        x_ar_2, y_ar_2 = from_input_to_numpy_array(intxtfile2,  max_number_iterations, 1, 0)
        print(x_ar_2)
        print(y_ar_2)
        x_ar_3, y_ar_3 = from_input_to_numpy_array(intxtfile3,  max_number_iterations, 1, 0)
        print(x_ar_3)
        print(y_ar_3)
        # color list prepare this assigns a number to each unique string x argument
        colist1, colist2, colist3, massimo = color_prepare(x_ar_1, x_ar_2, x_ar_3)
        print('colist1:\n', colist1, '\ncolist2:\n', colist2, '\ncolist1:\n', colist3)
        # get a color map here used my own
        top = cm.get_cmap('jet', (massimo))
        inbetween = cm.get_cmap('PiYG', massimo)
        bottom = cm.get_cmap('Pastel2', (massimo*2))
        newcolors = np.vstack((top(np.linspace(0, 1, (massimo))), inbetween(np.linspace(0, 1, massimo)), bottom(np.linspace(0, 1, massimo))))
        newcmp = ListedColormap(newcolors, name='Bubba')
        # get normalize function (takes data in range [vmin, vmax] -> [0, 1])
        my_norm = Normalize(vmin=0, vmax=(massimo))
        # this  is to put multiple plotsin one image
        fig, ax = plt.subplots(1,3)
        #plt.rcParams['font.size'] = '18'
        ax[0].barh(x_ar_1, y_ar_1, color=newcmp(my_norm(colist1)), edgecolor='gray')
        #ax[0].get_xticklabels().set_fontsize(16)
        #plt.subplot(1,3,2)
        ax[1].barh(x_ar_2, y_ar_2, color=newcmp(my_norm(colist2)), edgecolor='gray')
        ax[1].set(yticks=[])
        #plt.subplot(1,3,3)
        ax[2].barh(x_ar_3, y_ar_3, color=newcmp(my_norm(colist3)), edgecolor='gray')
        ax[2].set(yticks=[])
        fig.subplots_adjust(left=0.04, right=0.98, top=0.97, bottom=0.06)
        #plt.tight_layout()
        plt.show()
