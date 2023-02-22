#!/usr/bin/env python3

import sys
import numpy as np
from matplotlib import pyplot
import pandas as pd
import seaborn as sns

if __name__ == "__main__":
    try:
        intxtfile1 = sys.argv[1]
        intxtfile2 = sys.argv[2]
        #intxtfile2 = sys.argv[3]
    except:
        raise SystemExit
    else:
        in_pd1 = pd.read_csv(intxtfile1, sep='\t')
        in_pd2 = pd.read_csv(intxtfile2, sep='\t')
        volumes_np1 = in_pd1['avg Volume'].to_numpy()
        volumes_np2 = in_pd2['avg Volume'].to_numpy()
        volumes_mean1 = np.mean(volumes_np1)
        volumes_mean2 = np.mean(volumes_np2)
        
        #bins = np.linspace(110, 180, 36)
        
        pyplot.title('Vertebrate helix volume distribution')
        pyplot.xlabel('nm\xb3 * 10\xaf\xb3')
        pyplot.ylabel('n. of proteins')

        sns.histplot(volumes_np1, kde=True, binwidth=2, binrange=[120,160], label='TM')
        sns.histplot(volumes_np2, kde=True, binwidth=2, binrange=[120,160], label='IM')
        pyplot.legend(loc='upper right')

        pyplot.show()
