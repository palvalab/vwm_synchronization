# -*- coding: utf-8 -*-
"""
@author: hamed
"""

import sys
sys.path.append('L:\\nttk-data3\\palva\\VWM2\\_Population\\npTDMS')

from nptdms import TdmsFile
import numpy as np
import os
from scipy import stats
import pandas as pd
from numpy import genfromtxt
#import pingouin as pg

''' list of functions: 
    
    get_IM_from_group_csv: describe function
    
'''


def get_IM_from_group_csv(group_csv_path, frequency, start_index, stop_index):   
    '''
    description    
    
    argument:       description of argument  
    '''
    filename = ' Phase-Phase 1-1 Original No-Surrogates cPLV Lag_1-0 Lag=1.000 Low = {}hi={} parc2009_200AFSsign-stat.csv'.format(frequency, frequency)
    fpath = group_csv_path + filename
    IM = np.genfromtxt(fpath, delimiter = ';')    
    IM = IM[start_index:stop_index, :]
    all_NaNs = np.isnan(IM)
    IM[all_NaNs] = 0
    
    return IM


def read_tdms(fpath, parcel_nbs, TW):
    graphdata = TdmsFile.read(fpath)
    
    for i in range(parcel_nbs):
        channel_name = 'Edge Data (CSGreim)§{}__imag'.format(i)
        
        group = graphdata[TW]
        channel = group[channel_name]
        col = channel[:]
        if i == 0:
            IM = col
        else:
            IM = np.column_stack((IM, col))
    if np.shape(IM)[0] != parcel_nbs:
        print('Given parcel nb not same as in Graph Data')
    IM = np.abs(IM)
    
    return IM


def discard_false_discoveries(stats, alpha, parcel_nbs, edge_nbs, is_symmetric = True, gdm_python = True):
    FD_nbs = int((alpha / 2) * edge_nbs) #number of p values to remove
    sig_data = [] #list that will contain all p values
    
    if gdm_python:
        stats = np.reshape(stats, (parcel_nbs, parcel_nbs)) # reshape from flattened back to 2D
    if is_symmetric:
        i_lower = np.tril_indices(parcel_nbs, -1) # lower triangular indices
        stats[i_lower] = stats.T[i_lower] # copy the upper triangular data to the empty lower triangular section
       
    stats = np.ndarray.flatten(stats)
    for edge in range(edge_nbs):
        if stats[edge] != 0: #this statistical matrix is thresholded so only p values below alpha remain and rest are zero
            sig_data.append(stats[edge])
    
    sig_data.sort(reverse = True)
    if FD_nbs < len(sig_data):
        threshold_p = sig_data[FD_nbs] #obtain the p value for which values greater would be removed
    else:
        threshold_p = min(sig_data)
    for edge in range(edge_nbs):
        if stats[edge] >= threshold_p:
            stats[edge] = 0
    stats = np.reshape(stats, (parcel_nbs, parcel_nbs))
    
    return stats


def create_FDR_mask(stat_matrix, alpha, nb_of_tests):
    stat_matrix = abs(stat_matrix) #abs because includes pos and neg tail
    FD_nbs = int(alpha * nb_of_tests) #get nb of pvals to discard
    
    pval_indices = np.nonzero(stat_matrix) #get indices of the p values
    pvals = stat_matrix[pval_indices] #get a 1D array of the p values
    
    if pvals.shape[0] < FD_nbs: #if there are less p values than threshold
        fdr_mask = np.zeros(stat_matrix.shape) #create a mask of only zeros
    else:
        pvals = np.sort(pvals)[::-1] #perform reverse sorting to get highest p values first
        threshold = pvals[FD_nbs] #get p value that cuts off required number of large p valaues
        fdr_mask = np.where(stat_matrix < threshold, 1, 0) #create the mask based on threshold
        
    return fdr_mask


def indices_and_weights(target_index, morphing_op):
    weights =  morphing_op[:,target_index]
    indices = np.nonzero(weights)[0]
    weights = weights[weights != 0]
    
    zipped = list(zip(indices, weights))
    return zipped


def morphing_edge_data(IM, subject, morphing_op_path):
    morphing_op_path = '{}\\{}.csv'.format(morphing_op_path, subject)
    morphing_op = genfromtxt(morphing_op_path, delimiter = ';')
    
    target_data = np.zeros(shape=(200,200))
    
    for target_i in range(np.shape(target_data)[0]):
        sources_i = indices_and_weights(target_i, morphing_op)
        for target_j in range(np.shape(target_data)[1]):
            sources_j = indices_and_weights(target_j, morphing_op)
            
            area_weighted_average = []
            for source_i in sources_i:
                for source_j in sources_j:
                    area_weighted = IM[source_i[0]][source_j[0]] * source_i[1] * source_j[1]
                    area_weighted_average.append(area_weighted)
                    
            target_data[target_i][target_j] = sum(area_weighted_average)
    return target_data
    

def jackknife_resampling(subject_list, frequencies, master_folder, morphing_op_path, DEM_path, graphfolder, conditions, time_windows, output_folder, analysis_string, ratio = '1-1'):
    if ratio != '1-1':
        low_frequencies, high_frequencies = prune_frequencies(frequencies, ratio, master_folder)
    else:
        low_frequencies = frequencies
        high_frequencies = frequencies
    
    all_K_positives, all_K_negatives = gdm_python(subject_list, high_frequencies, master_folder, morphing_op_path, DEM_path, graphfolder, conditions, time_windows, ratio)
    float_freq = [float(f) for f in low_frequencies]
    
    df_pos = pd.DataFrame({'Frequency': float_freq, 'Property': all_K_positives})
    
    all_K_negatives = [-k for k in all_K_negatives]
    df_neg = pd.DataFrame({'Frequency': float_freq, 'Property': all_K_negatives})
    
    for i in range(len(subject_list)):
        jackknife_sample = list(subject_list)
        jackknife_sample.pop(i)
        
        print(subject_list[i])
        jackknife_pos, jackknife_neg = gdm_python(jackknife_sample, high_frequencies, master_folder, morphing_op_path, DEM_path, graphfolder, conditions, time_windows, ratio)
        
        column_header = 'Jackknife_{}'.format(i)
        df_pos[column_header] = jackknife_pos
        
        jackknife_neg = [-k for k in jackknife_neg]
        df_neg[column_header] = jackknife_neg
    
    condition = "_".join(conditions)    
    filename = output_folder + '\\cPLV_jackknife_' + condition + '_' + analysis_string + '_' + ratio + '_pos.csv'
    df_pos.to_csv(filename, sep = ';', index = False)
    filename =  output_folder + '\\cPLV_jackknife_' + condition + '_' + analysis_string + '_' + ratio + '_neg.csv'
    df_neg.to_csv(filename, sep = ';', index = False)
    

def prune_frequencies(frequencies, ratio, master_folder):
    cf_matrix = pd.read_csv('{}\\metadata\\cf_matrix_high_full.csv'.format(master_folder), sep = ';', dtype = 'object')
    
    low_frequencies = []
    high_frequencies = []
    for idx, freq in enumerate(frequencies):
        low_f = cf_matrix.loc[cf_matrix['high'] == freq, ratio].iloc[0]
        if isinstance(low_f, str) == True:
            low_frequencies.append(low_f)
            high_frequencies.append(freq)
            
    return low_frequencies, high_frequencies


def get_low_frequency(freq, ratio, master_folder):
    cf_matrix = pd.read_csv('{}\\metadata\\cf_matrix_high_full.csv'.format(master_folder), sep = ';', dtype = 'object')
    low_f = cf_matrix.loc[cf_matrix['high'] == freq[0], ratio].iloc[0]
    low_f = [low_f, low_f.replace('.', '-')]
    
    return low_f
    
    
def collapse_conditions(fpath, right_graphdata_folder, freq, conditions, graphfolder, original_parcel_nbs, time_window, master_folder, ratio = '1-1', merge = 'subtract', sub_bl = True):
    low_freq = get_low_frequency(freq, ratio, master_folder) if ratio != '1-1' else freq
    engine = 'Phase' if 'PAC' not in right_graphdata_folder else 'Amplitude'
    
    IMs = []
    
    for condition in conditions:
        fpath_condition = fpath + '\\' + right_graphdata_folder + '\\EdgeD{}-Phase {} Original No-Surrogates cPLV Lag_1-0 Lag=1.000 Low = {} {}_HIT {} {}Hz.tdms'.format(engine, ratio, low_freq[0], condition, graphfolder, freq[1])
        
        if sub_bl:
            IM_condition = read_tdms(fpath_condition, original_parcel_nbs, TW=time_window) - read_tdms(fpath_condition, original_parcel_nbs, TW='0')
        else:
            IM_condition = read_tdms(fpath_condition, original_parcel_nbs, TW=time_window)
        
        IMs.append(IM_condition)
        
    if merge == 'average':
        IM = sum(IMs) / len(conditions)
    else:
        IM = IMs[0] - sum(IMs[1:])
    
    return IM


def gdm_python(subject_list, frequencies, master_folder, morphing_op_path, DEM_path, graphfolder, conditions, time_windows, ratio):
    frequencies_dash = [f.replace('.', '-') for f in frequencies]                    
    zip_freq = list(zip(frequencies, frequencies_dash))
    
    original_parcel_nbs = 400
    alpha = 0.05
    
    all_K_positives = [0] * len(frequencies)
    all_K_negatives = [0] * len(frequencies)
    
    for time_window in time_windows:            
        count = 0
        for freq in zip_freq:
            all_IMs = []
            for subject in subject_list:              
                fpath = '{}\\{}\\GraphData'.format(master_folder, subject)
                directories = [dI for dI in os.listdir(fpath) if os.path.isdir(os.path.join(fpath,dI))]
                for directory in directories:
                    if graphfolder in directory:
                        right_graphdata_folder = directory
                
                IM = collapse_conditions(fpath, right_graphdata_folder, freq, conditions, graphfolder, original_parcel_nbs, time_window, master_folder, ratio)
                IM = morphing_edge_data(IM, subject, morphing_op_path)
                morphed_parcel_nbs = np.shape(IM)[0]
                edge_nbs = morphed_parcel_nbs * morphed_parcel_nbs
                
                if np.allclose(IM, IM.T):
                    IM = np.triu(IM)
                    is_symmetric = True
                else:
                    is_symmetric = False
                    
                flatten_IM = np.ndarray.flatten(IM)
                all_IMs.append(flatten_IM)
            
            stats_positive = np.zeros(edge_nbs,)
            stats_negative = np.zeros(edge_nbs,)
            
            for edge in range(np.shape(stats_positive)[0]):
                all_edges = []
                for subject_IM in all_IMs:
                    all_edges.append(subject_IM[edge])
                if all(e == 0 for e in all_edges) == False:
                    t, prob = stats.wilcoxon(all_edges)    
                    if np.mean(all_edges) > 0 and prob < alpha:
                        stats_positive[edge] = prob
                    if np.mean(all_edges) < 0 and prob < alpha:
                        stats_negative[edge] = prob
            
            weee_positive = discard_false_discoveries(stats_positive, alpha, morphed_parcel_nbs, edge_nbs, is_symmetric)
            weee_negative = discard_false_discoveries(stats_negative, alpha, morphed_parcel_nbs, edge_nbs, is_symmetric)
            
            DEM = np.genfromtxt(DEM_path, delimiter = ';')
            weee_positive = weee_positive * DEM
            weee_negative = weee_negative * DEM
            surviving_edges = np.shape(np.nonzero(DEM)[1])[0]
            
            K_positive = np.shape(np.nonzero(weee_positive)[1])[0] / surviving_edges
            K_negative = np.shape(np.nonzero(weee_negative)[1])[0] / surviving_edges
            
            all_K_positives[count] = all_K_positives[count] + K_positive
            all_K_negatives[count] = all_K_negatives[count] + K_negative
            print(freq[0])
            count += 1
        
    all_K_positives = [x/len(time_windows) for x in all_K_positives]
    all_K_negatives = [x/len(time_windows) for x in all_K_negatives]
    return all_K_positives, all_K_negatives
   
    
def k_to_excel(k_list, output_folder, frequencies, conditions, sign, ratio = '1-1'):
    condition = "_".join(conditions)
    float_freq = [float(f) for f in frequencies]
    df = pd.DataFrame({'Frequency': float_freq, 'Property': k_list})
    filename = '{}\\cPLV_python_{}_{}_{}.xlsx'.format(output_folder, ratio, condition, sign)
    df.to_excel(filename, index = False)
    
    
def get_averaged_IMs(subject_list, frequencies, master_folder, morphing_op_path, DEM_path, graphfolder, conditions, time_windows, sub_bl = True):
    frequencies_dash = [f.replace('.', '-') for f in frequencies]                    
    zip_freq = list(zip(frequencies, frequencies_dash))
    original_parcel_nbs = 400
    all_IMs = []
    
    for subject in subject_list:
        all_IMs_subject = []
        
        for idx_F, freq in enumerate(zip_freq):
            all_IMs_freq = []
            
            for time_window in time_windows:             
                fpath = '{}\\{}\\GraphData'.format(master_folder, subject)
                directories = [dI for dI in os.listdir(fpath) if os.path.isdir(os.path.join(fpath,dI))]
                for directory in directories:
                    if graphfolder in directory and 'PAC' not in directory and 'CFS' not in directory:
                        right_graphdata_folder = directory
                        
                IM = collapse_conditions(fpath, right_graphdata_folder, freq, conditions, graphfolder, original_parcel_nbs, time_window, master_folder, sub_bl = sub_bl)
                IM = morphing_edge_data(IM, subject, morphing_op_path)
                
                all_IMs_freq.append(IM)
                
            IM_avg_freq = np.mean(all_IMs_freq, axis=0)
            all_IMs_subject.append(IM_avg_freq)
        
        IM_avg_subject = np.mean(all_IMs_subject, axis=0)
        all_IMs.append(IM_avg_subject)
        print('Subject: {}'.format(subject))
        
    return all_IMs