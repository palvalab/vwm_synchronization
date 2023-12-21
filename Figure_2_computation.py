# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 11:23:38 2023

@author: hhaque
"""


master_folder = 'path_to_extracted_data' # set to directory of extracted data


source_directory    = '{}\\code\\utils'.format(master_folder)

import sys
sys.path.append(source_directory)
import vwm_core
import pandas as pd
import numpy as np


def compute_and_append_GS(df, frequencies, mask_fname, band, subject_list, condition, time_windows, master_folder, morphing_op_path, DEM_path, graphfolder):
    mask = np.genfromtxt(mask_fname, delimiter=';')
    mask = np.where(mask != 0, 1, 0)
    surviving_edges = np.count_nonzero(mask)

    all_IMs = vwm_core.get_averaged_IMs(subject_list, frequencies, master_folder, morphing_op_path, DEM_path, graphfolder, [condition], time_windows, sub_bl=False)

    for IM in all_IMs:
        IM = IM * mask
        gs = np.sum(IM) / surviving_edges

        df.loc[len(df)] = [condition, gs]

    return df


frequencies = ['3.00', '3.28', '3.68', '4.02', '4.29', '4.52', '5.05',
               '5.63', '6.00', '6.56', '7.40', '8.05', '8.57', '9.03',
               '10.09', '11.25', '12.00', '13.06', '14.89', '16.00',
               '17.14', '18.00', '20.10', '22.50', '24.00', '25.95',
               '30.00', '34.29', '36.00', '40.00', '45.00', '51.43',
               '60.00', '72.00', '80.00', '90.00', '102.86', '120.00']

DEM_path = '{}\\metadata\\current_DEM.csv'.format(master_folder)
time_windows = ['1', '2']
subject_list = ['S{:03d}'.format(i) for i in range(1, 21)]
morphing_op_path = '{}\\metadata\\Morphing_operators'.format(master_folder)
output_folder = '{}\\plot_data\\Figure_2\\Wilcoxon'.format(master_folder)
graphfolder = '3x2_Feature_x_Obj'
analysis_string = '1-1_PS'


'''Perform statistics on retention networks'''

#Shape
list_of_conditions = [['Shape2'], ['Shape4'], ['Shape4', 'Shape2']]
for conditions in list_of_conditions:
    vwm_core.jackknife_resampling(subject_list, frequencies, master_folder, morphing_op_path, DEM_path, graphfolder, conditions, time_windows, output_folder, analysis_string)


#Color
list_of_conditions = [['Color2'], ['Color4'], ['Color4', 'Color2']]
for conditions in list_of_conditions:
    vwm_core.jackknife_resampling(subject_list, frequencies, master_folder, morphing_op_path, DEM_path, graphfolder, conditions, time_windows, output_folder, analysis_string)


#Location
list_of_conditions = [['Spatial2'], ['Spatial4'], ['Spatial4', 'Spatial2']]
for conditions in list_of_conditions:
    vwm_core.jackknife_resampling(subject_list, frequencies, master_folder, morphing_op_path, DEM_path, graphfolder, conditions, time_windows, output_folder, analysis_string)



'''Obtain individual graph strength'''

conditions = ['Shape', 'Color', 'Spatial']
bands = ['Theta', 'Alpha']
list_of_frequencies = [['6.00', '6.56', '7.40', '8.05'], ['11.25', '12.00', '13.06']]

# Inputs
fpath = '{}\\output\\graph_strength'.format(master_folder)

graphfolder = '3x2_Load_avg'
mask_path = '{}\\metadata\\difference_masks'.format(master_folder)


for idb, band in enumerate(bands):
    df = pd.DataFrame(columns=['Condition', 'Mean_GS'], dtype=object)
    frequencies = list_of_frequencies[idb]
    mask_fname = '{}\\{}.csv'.format(mask_path, band)
    for condition in conditions:
        df = compute_and_append_GS(df, frequencies, mask_fname, band, subject_list, condition, time_windows, master_folder, morphing_op_path, DEM_path, graphfolder)

    df_fname = '{}\\GS_{}.csv'.format(fpath, band)
    df.to_csv(df_fname, sep = ';', index = False)
