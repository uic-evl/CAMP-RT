# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 10:08:52 2019

@author: Andrew
"""
from PatientSet import PatientSet
from ErrorChecker import ErrorChecker
from Constants import Constants
import numpy as np
import matplotlib.pyplot as plt

#def show_pca_scatterplot(db, reduction_fun = Rankings.pca, on = 'dose'):
#    residual = db.evaluate()
#    if on == 'distance':
#        x = db.gen_tumor_distance_matrix()[0]
#    elif on == 'error':
#        x = np.absolute(residual['differences'])
#    else:
#        on = 'dose'
#        x = db.doses
#    components = reduction_fun(x)
#    scale = np.mean(np.absolute(residual['differences']), axis = 1)**2
#    colors = db.get_class_list()
#    plt.scatter(components[:,0], components[:,1], scale, c = colors, cmap='Accent')
#    plt.title(on + '-distribution PCA')
#    plt.xlabel('PC 1')
#    plt.ylabel('PC 2')


db = PatientSet(root = 'data\\patients_v3*\\',
                outliers = [],
                class_name = None, 
                use_distances = False)
ErrorChecker().check_missing_organs(db)
#db.export(patient_data_file = 'data\\patient_dataset_v3.json', score_file = 'data\\all_ssim_scores_v3.csv')
#result = db.evaluate()
#print(result['mean_error'])
