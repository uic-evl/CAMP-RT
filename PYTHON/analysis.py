# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 10:08:52 2019

@author: Andrew
"""
from PatientSet import PatientSet, Rankings
from Patient import Patient
from Constants import Constants
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def show_pca_scatterplot(db, reduction_fun = Rankings.pca, on = 'dose'):
    residual = db.evaluate()
    if on == 'distance':
        x = db.gen_tumor_distance_matrix()[0]
    elif on == 'error':
        x = np.absolute(residual['differences'])
    else:
        on = 'dose'
        x = db.doses
    components = reduction_fun(x)
    scale = np.mean(np.absolute(residual['differences']), axis = 1)**2
    colors = db.get_class_list()
    plt.scatter(components[:,0], components[:,1], scale, c = colors, cmap='Accent')
    plt.title(on + '-distribution PCA')
    plt.xlabel('PC 1')
    plt.ylabel('PC 2')

def correlate_avg(db, key1, key2):
    avg = db.get_average_patient_data()
    return(np.correlate(avg[key1], avg[key2]))

def organ_droput_test(db, weights = np.ones((Constants.num_organs,))):
    from copy import copy
    errors = {}
    for x in range(Constants.num_organs):
        if weights[x] == 0:
            continue
        new_weights = copy(weights)
        new_weights[x] = 0
        result = db.evaluate(weights = new_weights)
        print('dropping ', Constants.organ_list[x], ' ', base_error - result['mean_error'])
        errors[Constants.organ_list[x]] = base_error - result['mean_error']
    return errors

db = PatientSet(patient_set = None, root = 'data\\patients_v*\\',
                outliers = Constants.v2_bad_entries + Constants.v3_bad_entries, class_name = None)
boost = [ (1 if p.neck_boost else 0) for p in db.get_patients()]
ht = [ (1 if (p.group == 2) else 0) for p in db.get_patients()]
#db.export(patient_data_file = 'data\\patient_dataset_v23.json', score_file = 'data\\all_ssim_scores_v23.csv')
#result = db.evaluate()
#print(result['mean_error'])
