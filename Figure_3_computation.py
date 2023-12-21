# -*- coding: utf-8 -*-
"""
Created on Fri May 14 09:39:21 2021

@author: hhaque
"""

master_folder = 'path_to_extracted_data' # set to directory of extracted data

source_directory    = '{}\\code\\utils'.format(master_folder)

import sys
sys.path.append(source_directory)
import vwm_core
import pandas as pd
import numpy as np
import itertools


subject_list = ['S{:03d}'.format(i) for i in range(1, 21)]

#frequencies = ['6.00', '6.56', '7.40', '8.05'] #theta
frequencies = ['11.25', '12.00', '13.06'] #alpha

time_windows = ['1', '2']

system_list = ['Dorsal', 'Ventral', 'V4-V8', 'L0', 'V1-V4']
system_indices = np.triu_indices(len(system_list))


list_of_conditions = [['Shape'], ['Color'], ['Spatial']]

morphing_op_path = '{}\\metadata\\Morphing_operators'.format(master_folder)
DEM_path = '{}\\metadata\\current_DEM.csv'.format(master_folder)
graphfolder = '3x2_Load_avg'


df = pd.DataFrame(columns = ['Subject', 'Condition', 'Edge', 'Mean strength']) # structure of dataframe that would contain all edge weights

for idx_cond, conditions in enumerate(list_of_conditions):
    all_IMs = vwm_core.get_averaged_IMs(subject_list, frequencies, master_folder, morphing_op_path, DEM_path, graphfolder, conditions, time_windows)
    
    ROI_path = '{}\\metadata\\Functional_boundaries.csv'.format(master_folder)
    roi_df = pd.read_csv(ROI_path, sep = '\t', header = None)
    roi_list = list(roi_df[0])
    
    for idx_subj, subject in enumerate(subject_list): #loop through each subject
        subject_IM = all_IMs[idx_subj]
        np.fill_diagonal(subject_IM, 0)
        
        for i in range(len(system_indices[0])):
            idx_source = system_indices[0][i]
            idx_target = system_indices[1][i]
            
            source_indices = [i for i, e in enumerate(roi_list) if e == idx_source] #get the indices within the source subsystem
            target_indices = [i for i, e in enumerate(roi_list) if e == idx_target] #ditto
                
            subsystem_pair_edges = [] #all edge strength between parcels of the two subsystems
            
            if source_indices != target_indices:
                for source_parcel in source_indices: #loop through parcels in source subsystem
                    for target_parcel in target_indices: #loop through parcels in parcel subsystem
                        edge = subject_IM[source_parcel][target_parcel]
                        if edge != 0:
                            subsystem_pair_edges.append(edge) #get edge for 'source parcel - target parcel' and append
            
            else: #to prevent edges from being selected twice when looking at within subsystem connections
                index_pairs = list(set(itertools.combinations(source_indices, 2)))
                for index_pair in index_pairs:
                    edge = subject_IM[index_pair[0]][index_pair[1]]
                    if edge != 0:
                        subsystem_pair_edges.append(edge) #get edge for 'source parcel - target parcel' and append
                        
            mean_strength = np.mean(subsystem_pair_edges)
            
            edge_name = '{}-{}'.format(system_list[idx_source], system_list[idx_target])
            df.loc[len(df)] = [subject, conditions[0], edge_name, mean_strength]