# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 16:50:47 2021

@author: hhaque
"""

master_folder = 'path_to_extracted_data' # set to directory of extracted data

source_directory    = '{}\\code\\utils'.format(master_folder)

import sys
sys.path.append(source_directory)
import vwm_core
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statistics
from scipy import stats
import seaborn as sns


def threshold_by_q(IM_group, ERM, q):
    IM_group = IM_group * ERM
    
    all_sig = np.count_nonzero(IM_group)
    threshold_index = all_sig - round(IM_group.size * q)
    
    if threshold_index < 0:
        threshold_index = 0
    
    flaty = IM_group.flatten()
    indy = np.argsort(flaty)[::-1]
    threshold_value = flaty[indy[threshold_index]]
     
    return threshold_value

def weird_division(n, d):
    return n / d if d else 0

fpath_master = '{}\\'.format(master_folder)
morphing_op_path = '{}\\metadata\\Morphing_operators'.format(master_folder)
fpath_behavioral = '{}\\metadata\\behavioral_df.csv'.format(master_folder)
output_fpath = '{}\\plot_data\\Figure_5'.format(master_folder)

subject_list = ['S{:03d}'.format(i) for i in range(1, 21)]

load = '04'
sign = 'pos'

band = 'alpha'
frequencies = ['9.03', '10.09', '11.25']
    
ERM_path = '{}\\metadata\\current_DEM.csv'.format(master_folder)
ERM = np.genfromtxt(ERM_path, delimiter = ';')

df_master = pd.DataFrame(columns = ['Subject', 'Feature', 'iPLV'])

q = 0.00672

IM_all_masks = []

features = ['Location']
the_color = sns.color_palette('muted')[0]

for TW in ['1', '2']:
    for frequency in frequencies:
        IM_mask_features = np.zeros((200, 200))
        for feature in features:
            condition = feature + load + '_{}'.format(sign)
            fpath_group = '{}\\metadata\\__group_statistics_csv\\3x2_Feature_x_Obj EdgeD  Pearson_{}\\Phase-Phase 1-1 Original No-Surrogates cPLV Lag_1-0 Lag=1.000 Low = {}hi={} parc2009_200AFSsign-stat.csv'.format(master_folder, condition, frequency, frequency)
            
            IM_group = np.genfromtxt(fpath_group, delimiter = ';')
            if TW == '1':
                IM_group = IM_group[200:400, :]
            if TW == '2':
                IM_group = IM_group[400:600, :]
                
            IM_group = np.nan_to_num(IM_group)
            threshold_value = threshold_by_q(IM_group, ERM, q)
            
            IM_mask = np.where(IM_group > 0, 1, 0)
            IM_mask = IM_mask * ERM
            print('For {} and {}, the nb of significant edges is {} and K is {}'.format(frequency, condition, np.count_nonzero(IM_mask), np.count_nonzero(IM_mask) / 40000))
            
            IM_mask_features = IM_mask_features + IM_mask
        IM_mask_features = np.where(IM_mask_features == len(features), 1, 0)
        print('Combined mask has {} edges and K is {}'.format(np.count_nonzero(IM_mask_features), np.count_nonzero(IM_mask_features) / 40000))
        print('-----------------------------------------------------')
        IM_all_masks.append(IM_mask_features)
        
        
    for feature in features:
        if feature == 'Location':
            feature = 'Spatial'
        for subject in subject_list:
            iplv_band = []
            for idx, frequency in enumerate(frequencies):
                
                fpath = fpath_master + '{}\\GraphData'.format(subject)
                iplv_dir = [f.path for f in os.scandir(fpath) if f.is_dir()]
                iplv_dir = [x for x in iplv_dir if 'Feature_x_Obj' in x][0]
                fpath_tdms = '{}\\EdgeDPhase-Phase 1-1 Original No-Surrogates cPLV Lag_1-0 Lag=1.000 Low = {} {}{}_HIT 3x2_Feature_x_Obj {}Hz.tdms'.format(iplv_dir, frequency, feature, load[1:], frequency.replace('.','-'))
                
                IM_baseline = vwm_core.read_tdms(fpath_tdms, 400, '0')
                IM_baseline = vwm_core.morphing_edge_data(IM_baseline, subject, morphing_op_path)
                
                IM = vwm_core.read_tdms(fpath_tdms, 400, TW)
                IM = vwm_core.morphing_edge_data(IM, subject, morphing_op_path)
                
                IM = IM - IM_baseline
                
                IM = IM * IM_all_masks[idx]
                
                iplv_band.append(weird_division(np.sum(IM),np.count_nonzero(IM)))
            
            if TW == '1':
                row_dict = {'Subject': subject, 'Feature': feature, 'iPLV': statistics.mean(iplv_band)}
                df_master = df_master.append(row_dict, ignore_index=True)
            if TW == '2':
                both_iPLV = (statistics.mean(iplv_band) + df_master.loc[df_master['Subject'] == subject, 'iPLV'].item()) / 2
                df_master.loc[df_master['Subject'] == subject, 'iPLV'] = both_iPLV
        

df_behavioral = pd.read_csv(fpath_behavioral, sep = ',')
df_behavioral = df_behavioral.loc[(df_behavioral['Condition'] == feature) & (df_behavioral['Load'] == int(load[1:])) & (np.isnan(df_behavioral['HR']) == False)]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.tick_params(which='both', top=False, right=False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.scatter(df_behavioral['HR'], df_master['iPLV'], color = the_color)

m, b = np.polyfit(df_behavioral['HR'], df_master['iPLV'], 1)
plt.plot(df_behavioral['HR'], m*df_behavioral['HR'] + b, color = 'black', lw = 0.5)

r, prob = stats.pearsonr(df_behavioral['HR'], df_master['iPLV'])
plt.text(0.875,-0.03, 'r = {:0.2f}\np = {:0.2e}'.format(r, prob))

plt.xlabel('Hit Rate (HR)')
plt.ylabel('iPLV')

plt.title('{} {}'.format(feature, load))

output_fname = '{}{}_{}'.format(feature, load, band)
df_master.to_csv('{}\\{}.csv'.format(output_fpath, output_fname), sep = ';', index = False)