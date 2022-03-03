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
        i =  0
        x_array = np.empty([max_iter], dtype='U100')            # setting maximum length of dtring type to add to array otherwise str uses only first letter
        y_array = np.empty([max_iter], dtype=int)
        total_annotations_found= None
        for txtline in intxt:
            if i == max_iter:
                #print('if section')
                return x_array, y_array, total_annotations_found
            elif txtline[0] == '#':
                total_annotations_found = float(txtline.split('=')[1])  ## ADDED JUST IN CASE
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
        return x_array[0:i], y_array[0:i], total_annotations_found

def matcher_to_the_first(reference_numarray, label_numarray2, value_numarray2, label_numarray3, value_numarray3):
    dict_ref = {}
    for label_str in reference_numarray:
        dict_ref[label_str] = [0, 0]                # first zero is for 2 array and second zero is for third array
    len_2 = len(label_numarray2)
    len_3 = len(label_numarray3)
    regulator = None
    if len_2 >= len_3:
        regulator = len_2
    else:
        regulator = len_3
    for i in range(0, regulator):
        if i >= len_2:
            tmp = None
        else:
            for ref_key in dict_ref:
                if label_numarray2[i] == ref_key:
                    #print('hereis a match:', label_numarray2[i], ref_key, value_numarray2[i])
                    dict_ref[ref_key][0] = value_numarray2[i]
        if i >= len_3:
            tmp = None
        else:
            for elem in dict_ref:
                if label_numarray3[i] == elem:
                    #print('hereis a match:', label_numarray2[i], ref_key, value_numarray2[i])
                    dict_ref[elem][1] = value_numarray3[i]
    #print(dict_ref)
    matched_value2 = np.empty([len(reference_numarray)], dtype=int)
    matched_value3 = np.empty([len(reference_numarray)], dtype=int)
    n = 0
    for chiave in dict_ref:
        matched_value2[n] = dict_ref[chiave][0]
        matched_value3[n] = dict_ref[chiave][1]
        n += 1
    return matched_value2, matched_value3



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
        intxtfile1 = sys.argv[1]                #reference the only one on which the max number is applied the other only match to his firlds
        intxtfile2 = sys.argv[2]
        intxtfile3 = sys.argv[3]
        max_number_iterations = int(sys.argv[4])
    except:
        raise SystemExit
    else:
        x_ar_1, y_ar_1, tot1 = from_input_to_numpy_array(intxtfile1,  max_number_iterations, 1, 0)
        #print(x_ar_1)
        #print(y_ar_1)
        #print(tot1)
        x_ar_2, y_ar_2, tot2 = from_input_to_numpy_array(intxtfile2,  999999, 1, 0)
        #print(x_ar_2)
        #print(len(y_ar_2), len(x_ar_2))
        #print(y_ar_2)
        x_ar_3, y_ar_3, tot3 = from_input_to_numpy_array(intxtfile3,  999999, 1, 0)
        #print(x_ar_3)
        #print(len(y_ar_3), len(x_ar_3))
        #print(y_ar_3)

        # Matching of the other to files to the first one necessary to create an order between the different fields of the three files
        correct_y2, correct_y3 = matcher_to_the_first(x_ar_1, x_ar_2, y_ar_2, x_ar_3, y_ar_3)
        #print("ordered second file annotations\n", correct_y2)
        #print("ordered third file annotations\n", correct_y3)

        # used to create the percentages in a strange case scenario
        #print("TANGO1-like adjacent annotation with respect of model confidence")
        print("RIGHT")
        print("plp < 50, 50 <= plp < 75, plp >= 75")
        for n, label_all in enumerate(x_ar_1):
            print(correct_y3[n], '    %.3f%% ,' %(correct_y3[n]/tot3*100), correct_y2[n], '    %.3f%% ,' %(correct_y2[n]/tot2*100),  y_ar_1[n], '     %.3f%% ,' %(y_ar_1[n]/tot1*100), label_all)
        print("unique annotations found")
        print("total annotation found\n", int(tot1), ",", int(tot2), ",", int(tot3))
        """
        
        # color list prepare this assigns a number to each unique string x argument
        #colist1, colist2, colist3, massimo = color_prepare(x_ar_1, x_ar_2, x_ar_3)
        #print('colist1:\n', colist1, '\ncolist2:\n', colist2, '\ncolist1:\n', colist3)

        # get a color map here used my own
        top = cm.get_cmap('jet', 120)
        inbetween = cm.get_cmap('PiYG', 120)
        bottom = cm.get_cmap('Pastel2', 120)
        newcolors = np.vstack((top(np.linspace(0, 1, 120)), inbetween(np.linspace(0, 1, 120)), bottom(np.linspace(0, 1, 120))))
        newcmp = ListedColormap(newcolors, name='Bubba')

        # get normalize function (takes data in range [vmin, vmax] -> [0, 1])
        #my_norm = Normalize(vmin=0, vmax=(y_ar_1[0]))

        # this  is to put multiple plotsin one image
        #fig, ax = plt.subplots(1,3)

        #set some hyperparameters
        fig, ax = plt.subplots()
        plt.rcParams['font.size'] = '22'
        barWidth = 0.28                                                                # the width of the bars
        colist = np.linspace(0.0, 1.0, num=max_number_iterations)

        # Set position of bar on X axis
        br1 = np.arange(len(x_ar_1))
        br2 = [x + barWidth for x in br1]
        br3 = [x + barWidth for x in br2]
        #print(br1, br2, br3)
        #print(len(y_ar_1), len(correct_y2), len(correct_y3))

        
        #bar1 = ax.bar(br1, y_ar_1, color='tomato', width=barWidth, edgecolor='gray',  align='edge', label='TANGO1-like     unique annotations found = 6362     total annotation found = 192522')
        #bar1 = ax.bar(br1, y_ar_1, color='lawngreen', width=barWidth, edgecolor='gray',  align='edge', label='TM non TANGO1     unique annotations found = 7668     total annotation found = 225831')
        bar1 = ax.bar(br1, y_ar_1, color='dodgerblue', width=barWidth, edgecolor='gray',  align='edge', label='non TM     unique annotations found =  10124    total annotation found = 192248')
        #bar2 = ax.bar(br2, correct_y2, color='lawngreen', width=barWidth, edgecolor='gray',  align='edge', label='TM non TANGO1     unique annotations found = 7668     total annotation found = 225831')
        bar2 = ax.bar(br2, correct_y2, color='tomato', width=barWidth, edgecolor='gray',  align='edge', label='TANGO1-like     unique annotations found = 6362     total annotation found = 192522')
        #bar3 = ax.bar(br3, correct_y3, color='dodgerblue', width=barWidth, edgecolor='gray',  align='edge', label='non TM     unique annotations found =  10124    total annotation found = 192248')
        bar3 = ax.bar(br3, correct_y3, color='lawngreen', width=barWidth, edgecolor='gray',  align='edge', label='TM non TANGO1     unique annotations found = 7668     total annotation found = 225831')

        # Adding labels and ticks
        ax.set_title('Domain presence', fontweight ='bold', fontsize = 28)
        #ax.set_xlabel('Domain presence', fontweight ='bold', fontsize = 28)
        ax.set_ylabel('Number of protein with annotation', fontweight ='bold', fontsize = 24)
        #plt.xticks(br3, ['pmt_2', 'nuc bind atp', 'abc transpor', 'HIS kinase', 'tpr repeat', 'basic res', 'ion_trans', 'abc trans t1', 'pas', 'lactamase_b', 'acyl_transf_3', 'mfs', 'pac', 'prot kinase', 'response reg', '4-ASPphosphate', 'ftsx', 'cyclic nuc-bind', 'competence',  'sgnh', 'ank repeat', 'p-loop N3P hydrolase', 'signal transduct HIS kinase', 'ggdef', 'comec/rec2-related prot', 'hsp90-like atpase', 'duf4131', 'tgc', 'ATP bindsite', 'stas', 'zn(2)-c6 fungal', 'ef-hand', 'yfho prot fam'], fontweight ='bold', fontsize = 22)
        #plt.xticks(br3, ['TM helix', 'disorder', 'polar res', 'coiled coil', 'basic acidic res', 'chain', 'signal pept', 'pro res', 'non-term res', 'acid res', 'pmt_2', 'nuc bind atp', 'abc transpor', 'HIS kinase', 'tpr repeat', 'basic res', 'ion_trans', 'abc trans t1', 'pas', 'lactamase_b', 'acyl_transf_3', 'mfs', 'pac', 'prot kinase', 'response reg', '4-ASPphosphate', 'ftsx', 'cyclic nuc-bind', 'competence', 'pas', 'sgnh', 'ank repeat', 'p-loop N3P hydrolase', 'signal transduct HIS kinase', 'ggdef', 'comec/rec2-related prot', 'hsp90-like atpase', 'duf4131', 'tgc', 'ATP bindsite', 'stas', 'zn(2)-c6 fungal', 'ef-hand', 'yfho prot fam', 'ggdef'], fontweight ='bold', fontsize = 22)
        #plt.xticks(br3, ['TM helix', 'disorder', 'polar res', 'coiled coil', 'basic acidic res', 'chain', 'signal pept', 'non-term res', 'pro res', 'HIS kinase', 'prot kinase', 'acid res', 'ATP bindsite', 'pas', 'response reg', 'hamp', '4-ASPphosphate', 'nuc bind atp', 'pac', 'abc transpor', 'basic res', 'tpr repeat', 'ggdef', 'signal transduct HIS kinase', 'hsp90-like atpase', 'pas', 'p-loop N3P hydrolase', 'hsp90-like atpase', 'abc trans t1', 'ggdef', 'disulfide bond', 'eal', 'signal transduct resp reg', 'ftsx', 'ig-like fold', 'meth acc transduct', 'fibronectine iii', 'egf-like', 'prot kinase', 'egf-like', 'mfs', 'eal', 'ig-like', 'fibronectine iii', 'hamp' ], fontweight ='bold', fontsize = 22)
        plt.xticks(br3, ['disorder', 'polar res', 'coiled coil', 'basic acidic res', 'chain', 'signal pept', 'non-term res', 'pro res', 'acid res', 'basic res', 'p-loop N3P hydrolase', 'prot kinase', 'helicase atp-bind', 'helicase c-term', 'tpr repeat', 'HIS kinase', 'pas', 'zn fing c2h2', 'plug', 'carrier', 'o-pantetheine4PhosphoSER',  'tonb_dep_rec', 'pac', 'ATP bindsite', 'c2h2-type', 'response reg', '4-ASPphosphate', 'nuc bind atp', 'repeat wd', 'repeat ank', 'tetratricopeptide-like', 'integrase catalytic', 'ig-like fold', 'pas', 'fibronectine iii', 'reverse transcriptase', 'fibronectine iii', 'phosphopant. bind acp', 'acp-like superfam', 'aaa', 'ring-type', 'por_secre_tail', 'zn(2)-c6 fungal', 'cchc-type', 'pkd' ], fontweight ='bold', fontsize = 22)
        for tick in ax.get_xticklabels():
            tick.set_rotation(90)
        plt.yticks(fontsize = 22)
        ax.set_xlim(left=-0.2, auto=True)                           # adjust to decrease or increase space between y axis and bars
        ax.set_xlim(right=45.5)                                       # adjust to increase or decrease space between bars and end

        fig.subplots_adjust(left=0.08, right=0.98, top=0.96, bottom=0.25)
        #fig.tight_layout()
        plt.legend()
        #plt.savefig("M27-hits_normal_nonTM_55190_comparison_uniprot_annotation.png")
        plt.show()"""
    
