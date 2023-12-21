# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 14:50:25 2020

@author: hhaque
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl
import glob
from numpy import genfromtxt

def proportions_to_df(master_path, c_limit, HEsize, freqs):
    prospective_df = []
    for i in range(len(freqs)):
        fpath = master_path + freqs[i] + '*'
        fpath = glob.glob(fpath)
        if len(fpath) == 1:
            fpath = fpath[0]
        else:
            print('ERROR: More or less than one matching condition')
            
        tidied_freq = freqs[i].replace('-','.')
        tidied_freq = tidied_freq.replace('_','\n-')
        each_row = []
        
        each_row.append(tidied_freq)
        for j in ['99', '0', '1', '2']:
            AM_specific_path = '{}\\c_limit({})\\HEsizeThld[{}]_AM[{}].csv'.format(fpath, c_limit, HEsize, j)
            AM = genfromtxt(AM_specific_path, delimiter = ';')
            edge_nb = np.shape(np.nonzero(AM)[0])[0]

            if freqs[i] in ['6_8-05', '17-14_20-10']:
                edge_nb = -edge_nb

            each_row.append(edge_nb)
        prospective_df.append(each_row)
    df = pd.DataFrame(prospective_df, columns = ['frequency', 'Shared', 'Shape specific', 'Color specific', 'Location specific'])
    
    return df


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True

sns.set(style="white")

master_folder = 'path_to_extracted_data' # set to directory of extracted data

master_path = '{}\\plot_data\\Figure_4\\'.format(master_folder)
c_limit = '0.06'
HEsize = '5'
freqs = ['6_8-05', '11-25_13-06']        

colors = ['#69CC63', '#EE8549', '#4173CE', '#D25E62']

df = proportions_to_df(master_path, c_limit, HEsize, freqs)
df = df.set_index('frequency')
df = df.div(40000, axis =0)

ax = df.plot.bar(rot=0, color = colors[::-1], width = 0.7)
ax.tick_params(which='both', top=False, right=False)
ax.tick_params(which='major', length = 10, direction='in')
ax.tick_params(which='minor', length = 5, direction='in')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.ylabel('K')
plt.xlabel('Frequency bands')
plt.legend(bbox_to_anchor=(1.0, 1.0), prop={'size': 10}, frameon = False)

plt.title('Shared and specific edges for three features', size = 14)

plt.tight_layout()