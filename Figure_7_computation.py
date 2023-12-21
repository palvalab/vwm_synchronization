# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 12:56:30 2020

@author: hhaque
"""

master_folder = 'path_to_extracted_data' # set to directory of extracted data

source_directory    = '{}\\code\\utils'.format(master_folder)
import sys
sys.path.append(source_directory)
import vwm_core


subject_list = ['S{:03d}'.format(i) for i in range(1, 21)]

frequencies = ['6.00', '6.56', '7.40', '8.05', '8.57', '9.03',
               '10.09', '11.25', '12.00', '13.06', '14.89', '16.00',
               '17.14', '18.00', '20.10', '22.50', '24.00', '25.95',
               '30.00', '34.29', '36.00', '40.00', '45.00', '51.43',
               '60.00', '72.00', '80.00', '90.00', '102.86', '120.00']


morphing_op_path = '{}\\metadata\\Morphing_operators'.format(master_folder)
DEM_path = '{}\\metadata\\current_DEM.csv'.format(master_folder)
graphfolder = '3x2_Load_avg_PAC'
time_windows = ['1', '2']
output_folder = '{}\\plot_data\\Figure_7\\PAC'.format(master_folder)
analysis_string = 'PAC'
ratio = '1-2'


list_of_conditions = [['Shape'], ['Color'], ['Spatial']]

for conditions in list_of_conditions:
    vwm_core.jackknife_resampling(subject_list, frequencies, master_folder, morphing_op_path, DEM_path, graphfolder, conditions, time_windows, output_folder, analysis_string, ratio)