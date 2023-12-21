# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 11:55:43 2023

@author: hhaque
"""

import numpy as np
import os
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut


def get_complex(a, b):
    c = complex(a,b)
    return c
vcomplex = np.vectorize(get_complex)

def bin_to_npy(subject, condition, condition_set_name, frequency):
    fpath = '{}\\__trial_connectomes_csv\\{}\\{}  {} Hz\\'.format(subject, condition_set_name, condition, frequency)
    nb_of_trials = int(len(os.listdir(fpath)) / 2)
    for trial in range(nb_of_trials):
        filename_re = '{}  {}  {} Hz trial_{}  re.bin'.format(condition_set_name, condition, frequency, trial)
        filename_im = '{}  {}  {} Hz trial_{}  im.bin'.format(condition_set_name, condition, frequency, trial)
        fpath_re = fpath + filename_re
        fpath_im = fpath + filename_im
        IM_real = np.fromfile(fpath_re, dtype='>f')
        IM_imag = np.fromfile(fpath_im, dtype='>f')
        trial_IM = vcomplex(IM_real, IM_imag)
        baseline_IM = trial_IM[0:79800]
        firstTW = trial_IM[79800:159600]
        secondTW = trial_IM[159600:]
        
        if trial == 0:
            all_IM = np.stack([baseline_IM, firstTW, secondTW])
            all_IM = np.expand_dims(all_IM, axis = 0)
        if trial != 0:
            IM = np.stack([baseline_IM, firstTW, secondTW])
            IM = np.expand_dims(IM, axis = 0)
            
            all_IM = np.vstack([all_IM, IM])
            
    return all_IM
    
def get_iPLV_matrix(all_matrices, TW):
    IMs = all_matrices[:, TW, :] #extracting single TW from 3D (results in 2D)
    average = np.mean(IMs, axis = 0) #averages across all trials (results in 1D)
    iPLV_matrix = abs(average.imag) #taking iPLV of each element
    return iPLV_matrix

def get_binary_mask(all_matrices, TW, tie = 1000):
    baseline_iPLV = get_iPLV_matrix(all_matrices, 0) #getting 1D matrix of average iPLV values for baseline
    retention_iPLV = get_iPLV_matrix(all_matrices, TW) #getting 1D matrix of average iPLV values for a retention window
    difference_iPLV = np.subtract(retention_iPLV, baseline_iPLV) #difference between the two
    
    max_indices = np.argsort(difference_iPLV)[::-1]
    
    binary_mask = np.zeros(np.shape(difference_iPLV))
    for i in range(tie):
        binary_mask[max_indices[i]] = 1
    return binary_mask

def get_1D_sync_vector(input_matrix, binary_mask):
    input_matrix = input_matrix.imag
    input_matrix = np.where(input_matrix == 0, 1e-12, input_matrix)
    input_matrix = input_matrix * binary_mask #taking imag ONLY and then masking it
    output_matrix = input_matrix[input_matrix != 0] #removing elements with 0 and vector gets shorter
    return output_matrix
    
def get_feature_matrix(npy_file, TW, binary_mask):
    nb_trials = np.shape(npy_file)[0]
    for i in range(nb_trials):
        trial_matrix = npy_file[i, TW, :]
        trial_vector = get_1D_sync_vector(trial_matrix, binary_mask)
        if i == 0:
            feature_matrix = trial_vector
        else:
            feature_matrix = np.vstack((feature_matrix, trial_vector))
    return feature_matrix



'''Import trial IM and create feature array'''

subject_list = ['S{:03d}'.format(i) for i in range(1, 21)]

condition_set_name = '3x2_Feature_x_Obj'

load = '2'
conditions = ['Shape{}_HIT'.format(load), 'Color{}_HIT'.format(load), 'Spatial{}_HIT'.format(load)]

frequencies = ['3.00', '3.28', '3.68', '4.02', '4.29', '4.52', '5.05',
               '5.63', '6.00', '6.56', '7.40', '8.05', '8.57', '9.03',
               '10.09', '11.25', '12.00', '13.06', '14.89', '16.00',
               '17.14', '18.00', '20.10', '22.50', '24.00', '25.95',
               '30.00', '34.29', '36.00', '40.00', '45.00', '51.43',
               '60.00', '72.00', '80.00', '90.00', '102.86', '120.00']


for frequency in frequencies:
    frequency = frequency.replace('.', '-')
    print(frequency)
    all_matrices = []
    for subject in subject_list:
        for idx, condition in enumerate(conditions):
            npy_matrix = bin_to_npy(subject, condition, condition_set_name, frequency)
            
            if idx == 0:
                subject_matrices = npy_matrix
            if idx > 0:
                subject_matrices = np.vstack((subject_matrices, npy_matrix))
                
        all_matrices.append(subject_matrices)
        print('{} added to all_matrices for {}'.format(subject, frequency))
        
    for TW in [1, 2]:
        binary_mask = get_binary_mask(all_matrices, TW)
        
        for id_s, subject in enumerate(subject_list):
            subject_IMs_allTW = all_matrices[id_s]
            
            subject_IMs = subject_IMs_allTW[:, TW, :]
            subject_IMs_baseline = subject_IMs_allTW[:, 0, :]
            
            trial_nbs = int(subject_IMs.shape[0] / len(conditions))
            
            output_fpath = '\\Single_trial_matrices\\{}\\{}'.format(subject, condition_set_name)
            if not os.path.exists(output_fpath):
                os.makedirs(output_fpath)
            
            for id_c, condition in enumerate(conditions):
                if id_c == 0:
                    condition_IM = subject_IMs[0:trial_nbs, :]
                    condition_IM_baseline = subject_IMs_baseline[0:trial_nbs, :]
                if id_c == 1:
                    condition_IM = subject_IMs[trial_nbs:(trial_nbs*2), :]
                    condition_IM_baseline = subject_IMs_baseline[trial_nbs:(trial_nbs*2), :]
                if id_c == 2:
                    condition_IM = subject_IMs[(trial_nbs*2):, :]
                    condition_IM_baseline = subject_IMs_baseline[(trial_nbs*2):, :]
            
                feature_matrix = get_feature_matrix(condition_IM, binary_mask)
                
                feature_matrix_baseline = get_feature_matrix(condition_IM_baseline, binary_mask)
                feature_matrix_baseline = np.mean(feature_matrix_baseline, axis = 0)
                
                feature_matrix = feature_matrix - feature_matrix_baseline
            
                feature_matrix_output = '{}\\{}_TW{}_{}_feature_array.npy'.format(output_fpath, condition, TW, frequency)
                np.save(feature_matrix_output, feature_matrix)



'''Perform random forest classification'''

fpath = '\\Single_trial_matrices'

conditions = ['Shape4_HIT', 'Color4_HIT', 'Spatial4_HIT']

for subject in subject_list:
    print(subject)
    output_name = 'Synchrony_three_features_load4_{}'.format(subject)
    
    list_of_accuracies = []
    
    for TW in [1, 2]:
        classification_accuracy = []
        
        for frequency in frequencies:
            frequency = frequency.replace('.', '-')
            fpath = 'Single_trial_matrices\\{}\\{}'.format(subject, condition_set_name)
            
            for idx, condition in enumerate(conditions):
                npy_file = '{}\\{}_TW{}_{}_feature_array.npy'.format(fpath, condition, TW, frequency)
                npy_matrix = np.load(npy_file)
                
                y_condition = [idx] * np.shape(npy_matrix)[0]
                
                if idx == 0:
                    X = npy_matrix
                    y = y_condition
                if idx > 0:
                    X = np.concatenate([X, npy_matrix])
                    y.extend(y_condition)
                    
            y = np.array(y)
               
            corr_predict = 0
            wrong_predict = 0
            
            loo = LeaveOneOut()
            for train_index, test_index in loo.split(X):
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = y[train_index], y[test_index]
                
                clf = RandomForestClassifier(random_state=0)
                clf.fit(X_train, y_train)
                X_predicted = clf.predict(X_test)
                
                if np.array_equal(X_predicted, y_test) == True:
                    corr_predict += 1
                else:
                    wrong_predict += 1
            
            if (corr_predict + wrong_predict) == np.shape(X)[0]:
                print('{} is done'.format(frequency))
            else:
                print('Hmm, sth went wrong')
                
            classification_accuracy.append(corr_predict / np.shape(X)[0])
            
        list_of_accuracies.append(classification_accuracy)
        
    df = pd.DataFrame(list_of_accuracies)
    df = df.transpose()
    df.columns = ['Early', 'Late']
    df['Frequency'] = frequencies
    df = df[['Frequency', 'Early', 'Late']]
    
    output_fpath = 'Classification_output\\'
    output_fpath = output_fpath + output_name + '.csv'
    df.to_csv(output_fpath, index = False, sep = ',')
