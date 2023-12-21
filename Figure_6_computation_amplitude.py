# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 12:03:21 2023

@author: hhaque
"""

import numpy as np
import os
import copy
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut


'''Import trial amplitudes and create feature array'''

subject_list = ['S{:03d}'.format(i) for i in range(1, 21)]

condition_set_name = 'Amplitude_1D'
conditions = ['Shape2_HIT', 'Shape4_HIT', 'Color2_HIT', 'Color4_HIT', 'Spatial2_HIT', 'Spatial4_HIT']

frequencies = ['3.00', '3.28', '3.68', '4.02', '4.29', '4.52', '5.05',
               '5.63', '6.00', '6.56', '7.40', '8.05', '8.57', '9.03',
               '10.09', '11.25', '12.00', '13.06', '14.89', '16.00',
               '17.14', '18.00', '20.10', '22.50', '24.00', '25.95',
               '30.00', '34.29', '36.00', '40.00', '45.00', '51.43',
               '60.00', '72.00', '80.00', '90.00', '102.86', '120.00']

for frequency in frequencies:
    frequency = frequency.replace('.', '-')
    
    for subject in subject_list:
        for idx, condition in enumerate(conditions):
            fpath = '{}\\__trial_connectomes_csv\\{}\\{}  {} Hz\\'.format(subject, condition_set_name, condition, frequency)
            nb_of_trials = int(len(os.listdir(fpath)))
            for trial in range(nb_of_trials):
                filename = '{}  {}  {} Hz trial_{}.bin'.format(condition_set_name, condition, frequency, trial)
                trial_path = fpath + filename
                amplitudes = np.fromfile(trial_path, dtype='>f')
                baselineTW = amplitudes[0:400]
                firstTW = amplitudes[400:800]
                secondTW = amplitudes[800:]
                
                if trial == 0:
                    firstTW_all = copy.deepcopy(firstTW)
                    secondTW_all = copy.deepcopy(secondTW)
                    baselineTW_all = copy.deepcopy(baselineTW)
                    
                if trial != 0:
                    firstTW_all = np.vstack([firstTW_all, firstTW])
                    secondTW_all = np.vstack([secondTW_all, secondTW])
                    baselineTW_all = np.vstack([baselineTW_all, baselineTW])
                    
            output_fpath = 'Single_trial_matrices\\{}\\{}'.format(subject, condition_set_name)
            if not os.path.exists(output_fpath):
                os.makedirs(output_fpath)
            
            np.save('{}\\{}_TW0_{}_feature_array.npy'.format(output_fpath, condition, frequency), baselineTW_all)
            np.save('{}\\{}_TW1_{}_feature_array.npy'.format(output_fpath, condition, frequency), firstTW_all)
            np.save('{}\\{}_TW2_{}_feature_array.npy'.format(output_fpath, condition, frequency), secondTW_all)
            
        print(subject)
    print(frequency)
    
    
    
'''Perform random forest classification'''
    
fpath = '\\Single_trial_matrices'

load = '4'
conditions = ['Shape{}_HIT'.format(load), 'Color{}_HIT'.format(load), 'Spatial{}_HIT'.format(load)]

TWs = [0]
TW_names = ['Baseline']
output_fpath = '\\Amplitude_1D'


for subject in subject_list:
    print(subject)
    output_name = 'amplitude_three_features_load{}_{}'.format(load, subject)
    
    list_of_accuracies = []
    
    for TW in TWs:
        classification_accuracy = []
        
        for frequency in frequencies:
            frequency = frequency.replace('.', '-')
            fpath = '\\Single_trial_matrices\\{}\\{}'.format(subject, condition_set_name)
            
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
    df.columns = TW_names
    df['Frequency'] = frequencies
    TW_names.insert(0, 'Frequency')
    df = df[TW_names]
    
    output_fpath = output_fpath + output_name + '.csv'
    df.to_csv(output_fpath, index = False, sep = ',')