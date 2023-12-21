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
from matplotlib.ticker import FixedLocator
import numpy as np
from scipy import stats


loads = ['4', '2']
windows = ['Early', 'Late']

window = 'Early'

subject_list = ['S{:03d}'.format(i) for i in range(1, 21)]

frequencies = ['3.00', '3.28', '3.68', '4.02', '4.29', '4.52', '5.05',
               '5.63', '6.00', '6.56', '7.40', '8.05', '8.57', '9.03',
               '11.25', '12.00', '13.06', '14.89', '16.00',
               '17.14', '18.00', '20.10', '22.50', '24.00', '25.95',
               '30.00', '34.29', '36.00', '40.00', '45.00', '51.43',
               '60.00', '72.00', '80.00', '90.00', '102.86', '120.00']

frequencies = [float(x) for x in frequencies]

colors_first = [(0.9333333333333333, 0.5215686274509804, 0.2901960784313726), 
                (0.2823529411764706, 0.47058823529411764, 0.8156862745098039), 
                (0.8392156862745098, 0.37254901960784315, 0.37254901960784315)]

sig_colors = [(0.8392156862745098, 0.37254901960784315, 0.37254901960784315), '#ED9898']

colors = sns.color_palette('muted')

fig = plt.figure(figsize = (3.75, 1.75))

sns.set_style('white')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True

wm_structures = [['Retention',''], ['Baseline','BL_']]


ax = fig.add_subplot(1, 2, 1)
for idy, structure in enumerate(wm_structures):
    for load in loads:
        df_average = pd.DataFrame()
        for subject in subject_list:
            fpath = '{}\\plot_data\\Figure_6\\Synchrony\\Synchrony_three_features_load{}_{}{}.csv'.format(master_folder, load, structure[1], subject)
            df = pd.read_csv(fpath, sep = ',')
            
            df_average[subject] = df[window]
            
        df_average['avg'] = df_average.mean(axis=1)
        
        #get the CI via jackknife resampling
        df_jackknife_all = pd.DataFrame()
        for subject in subject_list:
            df_jackknife = df_average.drop(['avg'], axis = 1)
            df_jackknife = df_average.drop([subject], axis = 1)
            df_jackknife_all['sans_{}'.format(subject)] = df_jackknife.mean(axis = 1)
        df_jackknife_all['CI_max'] = df_jackknife_all.max(axis = 1)
        df_jackknife_all['CI_min'] = df_jackknife_all.min(axis = 1)
        
        if load == '2':
            ax.plot(frequencies, df_average.avg, lw = 1, color = colors_first[idy], linestyle = '--')
        if load == '4':
            ax.plot(frequencies, df_average.avg, lw = 1, color = colors_first[idy], label = structure[0])
        plt.fill_between(frequencies, df_jackknife_all.CI_min, df_jackknife_all.CI_max, color = colors_first[idy], alpha = 0.25, linewidth = 0.01)
        
        if structure[0] == 'Baseline':
            ax.tick_params(which='both', top=False, right=False)
            ax.tick_params(which='major', length = 4, direction='in', labelsize = 6, width=0.6)
            ax.tick_params(which='minor', length = 2, direction='in', labelsize = 6, width=0.6)
            ax.set_xscale('log')
            ax.set_xticks([3, 5, 10, 20, 30, 50, 100])
            ax.set_xlim([3, 120])
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.xaxis.set_minor_locator(FixedLocator([4, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 40, 60, 70, 80, 90]))
            ax.xaxis.set_minor_formatter(plt.NullFormatter())
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            plt.xlabel('Frequency (Hz)', size = 7)
            ax.set_ylim([0.413, 0.649])
            plt.ylabel('Accuracy', size = 7)
            plt.title(window, size = 8)
            
df_difference = pd.DataFrame()

ax = fig.add_subplot(1, 2, 2)
for idx, load in enumerate(loads):
    for subject in subject_list:
        fpath_retention = '{}\\plot_data\\Figure_6\\Synchrony\\Synchrony_three_features_load{}_{}.csv'.format(master_folder, load, subject)
        fpath_baseline = '{}\\plot_data\\Figure_6\\Synchrony\\Synchrony_three_features_load{}_BL_{}.csv'.format(master_folder, load, subject)
        
        df_retention = pd.read_csv(fpath_retention, sep = ',')
        df_baseline = pd.read_csv(fpath_baseline, sep = ',')
        
        df_difference[subject] = df_retention[window] - df_baseline[window]
    df_difference['avg'] = df_difference.mean(axis=1)
    
    df_jackknife_all = pd.DataFrame()
    for subject in subject_list:
        df_jackknife = df_difference.drop(['avg'], axis = 1)
        df_jackknife = df_difference.drop([subject], axis = 1)
        df_jackknife_all['sans_{}'.format(subject)] = df_jackknife.mean(axis = 1)
    df_jackknife_all['CI_max'] = df_jackknife_all.max(axis = 1)
    df_jackknife_all['CI_min'] = df_jackknife_all.min(axis = 1)
    
    df_stats = df_difference.drop(['avg'], axis = 1)
    stats_positive = np.zeros(len(frequencies))
    stats_negative = np.zeros(len(frequencies))
    for i in range(len(frequencies)):
        differences = df_stats.values[i]
        statistic, prob = stats.wilcoxon(differences, mode = 'approx')
        if np.mean(differences) > 0 and prob < 0.001:
            stats_positive[i] = prob
        if np.mean(differences) < 0 and prob < 0.001:
            stats_negative[i] = prob
    
    if load == '2':
        ax.plot(frequencies, df_difference.avg, lw = 1, color = colors_first[2], label = 'Retention > Baseline', linestyle = '--')
    if load == '4':
        ax.plot(frequencies, df_difference.avg, lw = 1, color = colors_first[2], label = 'Retention > Baseline')
    
    
    plt.fill_between(frequencies, df_jackknife_all.CI_min, df_jackknife_all.CI_max, color = colors_first[2], alpha = 0.25, linewidth = 0.01)
    
    ceiling = plt.ylim()[1]
    for sig_index in np.nonzero(stats_positive)[0]:
        ax.hlines(y=ceiling, xmin=frequencies[sig_index-1], xmax=frequencies[sig_index+1], color = sig_colors[idx], lw = 1)
    
    if load == '2':
        ax.tick_params(which='both', top=False, right=False)
        ax.tick_params(which='major', length = 4, direction='in', labelsize = 6, width=0.6)
        ax.tick_params(which='minor', length = 2, direction='in', labelsize = 6, width=0.6)
        ax.set_xscale('log')
        ax.set_xticks([3, 5, 10, 20, 30, 50, 100])
        ax.set_xlim([3, 120])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.xaxis.set_minor_locator(FixedLocator([4, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 40, 60, 70, 80, 90]))
        ax.xaxis.set_minor_formatter(plt.NullFormatter())
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        floor = -0.027
        plt.fill_between(df_retention.Frequency, floor, 0.0, color = 'gray', alpha = 0.25, linewidth = 0.01)
        plt.ylim(ymin=floor)
        ax.set_ylim([floor, 0.086])
        plt.xlabel('Frequency (Hz)', size = 7)
        plt.ylabel('D - Accuracy', size = 7)
        plt.title(window, size = 8)
        
plt.tight_layout()