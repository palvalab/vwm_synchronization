# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 10:42:36 2017

@author: hhaque
"""

master_folder = 'path_to_extracted_data' # set to directory of extracted data


import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
from matplotlib.ticker import FixedLocator


def get_CI_limits(condition, ratio, engine):
    fpath = '{}\\plot_data\\Figure_7\\{}_python\\'.format(master_folder, engine)
    fname = 'cPLV_jackknife_{}_bothTWs_{}_pos.csv'.format(condition, ratio)
    fpath = fpath + fname

    df = pd.read_csv(fpath, sep = ';')
    
    jackknife_col_names = ['Jackknife_{}'.format(i) for i in range(20)]
    df_jackknife = df.set_index('Frequency')[jackknife_col_names]
    df_jackknife['CI_max'] = df_jackknife.max(axis=1)
    df_jackknife['CI_min'] = df_jackknife.min(axis=1)
    
    return df_jackknife


sns.set_style('white')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True

colors = sns.color_palette('muted')
fig = plt.figure(figsize = (7, 6))

features = ['Shape', 'Color', 'Spatial']

param_dict = {'engine': ['PAC', 'PAC'],
              'ratios': [['1-2'], ['1-3', '1-4', '1-5', '1-6', '1-7', '1-8', '1-9'], ['1-2'], ['1-3', '1-4', '1-5', '1-6', '1-7', '1-8']],
              'y_limit': [0.13, 0.05, 0.06, 0.03],
              'y_locator': [[0.00, 0.05, 0.1], [0.00, 0.02, 0.04], [0.00, 0.02, 0.04, 0.06], [0.00, 0.01, 0.02, 0.03]]}

count = 1
for p_idx in range(len(param_dict['engine'])):
    
    os.chdir('{}\\plot_data\\Figure_7\\{}_python'.format(master_folder, param_dict['engine'][p_idx]))
    threshold = 0.00672
    
    for i in range(3):
        condition = features[i]
    
        ax = fig.add_subplot(4, 3, count)
        ax.tick_params(which='both', top=False, right=False)
        ax.tick_params(which='major', length = 4, direction='in', labelsize = 6, width=0.6)
        ax.tick_params(which='minor', length = 2, direction='in', labelsize = 6, width=0.6)
        ax.set_xscale('log')
        ax.set_xticks([3, 5, 10, 20, 30, 50])
        ax.set_xlim([3, 60])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.xaxis.set_minor_locator(FixedLocator([4, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 40, 60]))
        ax.xaxis.set_minor_formatter(plt.NullFormatter())
        ax.set_ylim([0.0, param_dict['y_limit'][p_idx]])
        ax.yaxis.set_major_locator(FixedLocator(param_dict['y_locator'][p_idx]))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        
        if count in [1, 2, 3]:
            the_title = features[i]
            if features[i] == 'Spatial':
                the_title = 'Location'
            plt.title(the_title, size = 8)
            
        fig.subplots_adjust(hspace= .3)
        
        if count in [1, 4, 7, 10]:
            plt.ylabel('K', size = 7)
        if count in [4, 5, 6, 10, 11, 12]:
            plt.xlabel('Frequency (Hz)', size = 7)
             
        for idr, ratio in enumerate(param_dict['ratios'][p_idx]):
            pos = 'cPLV_jackknife_{}_bothTWs_{}_pos.csv'.format(condition, ratio)
            df = pd.read_csv(pos, sep = ';')
            
            df_jackknife = get_CI_limits(condition, ratio, param_dict['engine'][p_idx])
            
            if ratio == '1-2':
                ax.plot(df.Frequency, df.Property, lw = 1, color = colors[idr], label = ratio.replace('-', ':'))
                plt.fill_between(df.Frequency, df_jackknife.CI_min, df_jackknife.CI_max, color = colors[idr], alpha = 0.25, linewidth = 0.01)
            else:
                ax.plot(df.Frequency, df.Property, lw = 1, color = colors[idr+1], label = ratio.replace('-', ':'))
                plt.fill_between(df.Frequency, df_jackknife.CI_min, df_jackknife.CI_max, color = colors[idr+1], alpha = 0.25, linewidth = 0.01)
                
            if ratio in ['1-2', '1-3']:
                plt.fill_between([3, 60], -threshold, threshold, color = 'gray', alpha = 0.25, linewidth = 0.01)
                
        if count in [1, 4]:
            plt.legend(frameon = False, prop={'size': 6})
        
        count += 1