# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 11:50:40 2020

@author: hhaque
"""


master_folder = 'path_to_extracted_data' # set to directory of extracted data

source_directory    = '{}\\code\\utils'.format(master_folder)

import sys
sys.path.append(source_directory)
import numpy as np
from numpy import genfromtxt
import pandas as pd
import plot_functions as plots
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
from scipy import stats


'''Plot susbsystem matrices for each condition'''

sns.set_style('white')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

df = pd.read_csv('{}\\metadata\\Patch_Grouping_divisions.csv'.format(master_folder), sep = '\t', header = None)
yeo = list(df[6])
yeo = [int(x/2) for x in yeo]

nb_of_yeo = max(yeo) + 1

conditions = ['Shape', 'Color', 'Location']    #conditions
labels = ['Alpha']
freqs = ['11-25_13-06']

all_subsystems = []

total_edges = np.zeros((nb_of_yeo,nb_of_yeo))
for index in np.ndindex(nb_of_yeo,nb_of_yeo):
    total_edges[index] = yeo.count(index[0]) * yeo.count(index[1])
    
ERM_path = '{}\\metadata\\current_DEM.csv'.format(master_folder)
ERM = np.genfromtxt(ERM_path, delimiter = ';')

for index in np.ndindex(200, 200):
    if ERM[index] == 0:
        total_edges[yeo[index[0]]][yeo[index[1]]] -= 1

for freq in freqs:
    for condition in conditions:
        
        fpath = '{}\\plot_data\\Figure_3\\Single_condition'.format(master_folder)
        fpath = fpath + '\\{}_pos_{}.csv'.format(condition, freq)

        IM = genfromtxt(fpath, delimiter = ';')
        
        max_indices = []
        for i in np.argsort(-IM, axis = None):
            index = np.unravel_index(i, np.shape(IM))
            if IM[index[0]][index[1]] > 0:
                max_indices.append(index)
        
        strongest_edges = []
        for pair in max_indices:
            strongest_edges.append(tuple((yeo[pair[0]], yeo[pair[1]])))    
        
        subsystems = np.zeros((nb_of_yeo,nb_of_yeo))
        for edge in strongest_edges:
            subsystems[edge[0], edge[1]] += 1
            
        subsystems = subsystems / total_edges
        
        subsystems = np.delete(subsystems, [0, 1, 4, 6], 0)
        subsystems = np.delete(subsystems, [0, 1, 4, 6], 1)
        
        new_index = [7, 6, 5, 4, 3, 1, 0, 2]
        subsystems = subsystems[np.ix_(new_index, new_index)]
        
        all_subsystems.append(subsystems)



systems = ['Early', 'Late', 'V4,V8', 'Vent', 'Dors', 'VAN', 'DAN', 'FPN']

titles = ['Shape', 'Color', 'Location']
for i in range(len(titles)):
    titles.append('')

ylabel_text = [''] * 3
    
plots.plot_heatmaps_pruned(all_subsystems, titles = titles, N_cols = 3, xticks = list(range(len(systems))), yticks = list(range(len(systems))),
                            ylabel = ylabel_text, zlabel = 'K', fontsizeL = 13, xticklabels=systems, yticklabels=systems, cmap = 'Reds',
                            xticklab_rot=90)
                           


'''Plot subsystem edge differences'''

band = 'alpha'

fpath = '{}\\plot_data\\Figure_3\\Subsystem_edges\\visual_edges_load_averaged_{}.csv'.format(master_folder, band)

df = pd.read_csv(fpath, sep = ';')

subsystem_pairs = list(df['Edge'].unique())
pair_indices = np.tril_indices(len(subsystem_pairs))

subsystem_names = ['Dors-Dors', 'Dors-Vent', 'Dors-V4,V8', 'Dors-L0',
                   'Dors-V1', 'Vent-Vent', 'Vent-V4,V8', 'Vent-L0',
                   'Vent-V1', 'V4,V8-V4,V8', 'V4,V8-L0', 'V4,V8-V1',
                   'L0-L0', 'L0-V1', 'V1-V1']

f,(ax1,ax2,ax3, axcb) = plt.subplots(1,4, 
            gridspec_kw={'width_ratios':[1,1,1,0.08]}, figsize = (7, 2.52))
ax1.get_shared_y_axes().join(ax2,ax3)

conditions = ['Shape', 'Color', 'Spatial']

count = 1
for condition in conditions:
    df_condition = df.loc[df['Condition'] == condition]
    
    subsystem_pairs_matrix = np.zeros((len(subsystem_pairs), len(subsystem_pairs)))
    sig = np.full((len(subsystem_pairs), len(subsystem_pairs)), '')
    
    for i in range(len(pair_indices[0])):
        idx_edge1 = pair_indices[0][i]
        idx_edge2 = pair_indices[1][i]
        
        if idx_edge1 != idx_edge2:
            cat1 = df_condition[(df_condition['Edge'] == subsystem_pairs[idx_edge1])]
            cat1 = list(cat1['Mean strength'])
            
            cat2 = df_condition[(df_condition['Edge'] == subsystem_pairs[idx_edge2])]
            cat2 = list(cat2['Mean strength'])
        
            iPLV_difference = np.mean(np.array(cat1) - np.array(cat2))
            subsystem_pairs_matrix[idx_edge1][idx_edge2] = iPLV_difference
        
            stat, prob = stats.wilcoxon(cat1, cat2, alternative = 'two-sided')
            #stat, prob = stats.ttest_rel(cat1, cat2, alternative = 'two-sided')
            if prob < 0.05:            
                verdict = 'significant'
                sig[idx_edge1][idx_edge2] = '*'
            else:
                verdict = 'n.s.'
            
            
    cmap = sns.diverging_palette(230, 20, as_cmap=True)
    mask = np.triu(np.ones_like(subsystem_pairs_matrix, dtype=bool))
    
    if condition == 'Shape': 
        sns.heatmap(subsystem_pairs_matrix, cmap = cmap, center = 0, annot = sig, fmt = '',
                         mask = mask, xticklabels = subsystem_names, yticklabels = subsystem_names,
                         annot_kws = {'color': 'black', 'size': 4}, cbar_kws={'label': 'D-iPLV', 'ticks': [-0.008, -0.004, 0.0, 0.004, 0.008]},
                         vmin = -0.008, vmax = 0.008, cbar = False, ax = ax1)
        
        ax1.tick_params(axis='both', labelsize=6)
        ax1.set_ylabel('Subsystem Edges A', fontsize = 7)
        #ax1.set_xlabel('Subsystem Edges B', fontsize = 7)
        ax1.set_title('{}'.format(condition), fontsize = 8)
        
    if condition == 'Color':
        sns.heatmap(subsystem_pairs_matrix, cmap = cmap, center = 0, annot = sig, fmt = '',
                         mask = mask, xticklabels = subsystem_names, yticklabels = False,
                         annot_kws = {'color': 'black', 'size': 4}, cbar_kws={'label': 'D-iPLV', 'ticks': [-0.008, -0.004, 0.0, 0.004, 0.008]},
                         vmin = -0.008, vmax = 0.008, cbar = False, ax = ax2)
        
        ax2.tick_params(axis='both', labelsize=6)
        ax2.set_xlabel('Subsystem Edges B', fontsize = 7)
        ax2.set_title('{}'.format(condition), fontsize = 8)
        
    if condition == 'Spatial':
        sns.heatmap(subsystem_pairs_matrix, cmap = cmap, center = 0, annot = sig, fmt = '',
                         mask = mask, xticklabels = subsystem_names, yticklabels = False,
                         annot_kws = {'color': 'black', 'size': 4}, cbar_kws={'label': 'Edge A > Edge B (D-iPLV)', 'ticks': [-0.008, -0.004, 0.0, 0.004, 0.008]},
                         vmin = -0.008, vmax = 0.008, ax = ax3, cbar_ax=axcb)
        
        ax3.tick_params(axis='both', labelsize=6)
        #ax3.set_xlabel('Subsystem Edges B', fontsize = 7)
        ax3.set_title('Location', fontsize = 8)
        
        axcb.tick_params(labelsize=6)
        axcb.yaxis.label.set_size(7)
        
    count += 1
    
plt.tight_layout()