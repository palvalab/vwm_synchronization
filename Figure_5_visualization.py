# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 16:50:47 2021

@author: hhaque
"""

master_folder = 'path_to_extracted_data' # set to directory of extracted data


import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
from matplotlib.ticker import FixedLocator
import matplotlib.lines as mlines
import numpy as np
from scipy import stats


'''Plot K value of significant correlation'''


os.chdir('{}\\plot_data\\Figure_5\\Correlation'.format(master_folder))

sns.set_style('white')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True

fig = plt.figure(figsize = (7.3, 1.85))

features = ['Shape', 'Color', 'Spatial']
colors = sns.color_palette('muted')

loads = ['02', '04']

threshold = 0.00672

count = 1
for i in range(len(features)):
    for load in loads:
        condition = features[i] + load

        pos = 'Pearson_{}_HR_pos.xls'.format(condition)
        neg = 'Pearson_{}_HR_neg.xls'.format(condition)

        if load == '02':
            ax = fig.add_subplot(1, 3, count)
            ax.tick_params(which='both', top=False, right=False)
            ax.tick_params(which='major', length = 4, direction='in', labelsize = 6, width=0.6)
            ax.tick_params(which='minor', length = 2, direction='in', labelsize = 6, width=0.6)
            ax.set_xscale('log')
            ax.set_xticks([3, 5, 10, 20, 30, 50, 100])
            ax.set_xlim([3, 120])
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.xaxis.set_minor_locator(FixedLocator([4, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 40, 60, 70, 80, 90]))
            ax.xaxis.set_minor_formatter(plt.NullFormatter())
            ax.set_ylim([-0.05, 0.05])
            ax.yaxis.set_major_locator(FixedLocator([-0.04, -0.02, 0.00, 0.02, 0.04]))
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            
            the_title = features[i]
            if features[i] == 'Spatial':
                the_title = 'Location'
    
            plt.title(the_title, size = 8)
            fig.subplots_adjust(hspace= .5)
            
            if count in [1]:
                plt.ylabel('K', size = 7, labelpad = 0.5)
            plt.xlabel('Frequency (Hz)', size = 7)
          
        for excel_file in ['pos', 'neg']:    
            df = pd.read_excel(eval(excel_file))
            df.columns = ['Frequency','Property']
           
            if load == '02':    
                ax.plot(df.Frequency, df.Property, lw = 1, linestyle = '--', color = colors[i], label = '_nolegend_')
            if load == '04':    
                ax.plot(df.Frequency, df.Property, lw = 1, color = colors[i], label = '_nolegend_')
            if load == '04-02':    
                ax.plot(df.Frequency, df.Property, lw = 1, linestyle = (0, (1, 1)), color = colors[i], label = '_nolegend_')
        
        if load == '04':
            count += 1
            
    plt.fill_between(df.Frequency, -threshold, threshold, color = 'gray', alpha = 0.25, linewidth = 0.01)
    
    x = 0.049
    x_neg = 0.045
    
    plt.plot([9.03, 11.25], [x, x], color = 'black', lw = 0.8) # alpha
    plt.plot([6.56, 8.05], [-x_neg, -x_neg], color = 'black', lw = 0.8) # theta
    
    if i == 2:
        dashed_line = mlines.Line2D([], [], color='black', linestyle = '--', lw = 1, label='Load 02')
        solid_line = mlines.Line2D([], [], color='black', lw = 1, label='Load 04')
        dotted_line = mlines.Line2D([], [], color='black', linestyle = (0, (1, 1)), lw = 1, label='Load 04-02')
        plt.legend(handles=[dashed_line, solid_line], frameon = False, prop={'size': 6})
            
plt.tight_layout()
plt.subplots_adjust(wspace=0.15)



'''Plot scatter of graph strength and HR'''


os.chdir('{}\\plot_data\\Figure_5\\Graph_strength'.format(master_folder))
fpath_behavioral = '{}\\metadata\\behavioral_df.csv'.format(master_folder)

fig = plt.figure(figsize = (7.3, 1.85))


band = 'alpha'
features = ['Shape02', 'Color04', 'Spatial04']

#band = 'theta'
#features = ['Shape02', 'Color02', 'Spatial02']


count = 1
for i in range(3):
    condition = features[i]
    feature = features[i][:-2]
    load = features[i][-2:]
    
    df_master = pd.read_csv('{}_{}.csv'.format(condition, band), sep = ';')
    
    df_behavioral = pd.read_csv(fpath_behavioral, sep = ',')
    df_behavioral = df_behavioral.loc[(df_behavioral['Condition'] == feature) & (df_behavioral['Load'] == int(load[1:])) & (np.isnan(df_behavioral['HR']) == False)]
    
    ax = fig.add_subplot(1, 3, count)
    ax.tick_params(which='both', top=False, right=False)
    ax.tick_params(which='major', length = 4, direction='in', labelsize = 6, width=0.6)
    ax.tick_params(which='minor', length = 2, direction='in', labelsize = 6, width=0.6)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylim([-0.042, 0.04])
    ax.yaxis.set_major_locator(FixedLocator([-0.04, -0.02, 0.00, 0.02, 0.04]))
    ax.set_xlim([0.45, 1.0])
    plt.scatter(df_behavioral['HR'], df_master['iPLV'], color = colors[i], s = 5)
    
    m, b = np.polyfit(df_behavioral['HR'], df_master['iPLV'], 1)
    plt.plot(df_behavioral['HR'], m*df_behavioral['HR'] + b, color = 'black', lw = 0.5)
    
    r, prob = stats.pearsonr(df_behavioral['HR'], df_master['iPLV'])
    
    if r < 0:
        y_pos = 0.8
    else:
        y_pos = 0.1
    plt.text(0.7,y_pos, 'r = {:0.2f}\np = {:0.2e}'.format(r, prob), transform=plt.gca().transAxes, size = 6)
    
    plt.xlabel('Hit Rate (HR)', size = 7)
    if count in [1]:
        plt.ylabel('Graph Strength', size = 7)
    
    plt.title('{} {}'.format(feature, load), size = 8)
    count += 1
    
plt.tight_layout()
