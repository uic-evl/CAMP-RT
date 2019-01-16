# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 09:19:09 2019

@author: Andrew
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json

patients_file = "data\\patients_SSIM_wDoses_wDists.json"
ssim_scores_file = "latest_results\\all_ssim_scores.csv"

with open(patients_file) as file:
    data = json.load(file)
data_map = {entry["ID_internal"]-1:entry for entry in data}
ssim_scores = {entry["ID_internal"]: entry["scores_ssim"] for entry in data}
ssim_rankings = {entry["ID_internal"]: entry['similarity'] for entry in data}
organ_list = list(data[0]['organData'].keys())
ssim_scores = pd.read_csv(ssim_scores_file,index_col=0)

def gen_dose_matrix(data):
    #takes the data from the json file an generates a num_patients by num_organs matrix of
    #mean radiation dosages. idexes are the internal id - 1
    dose_matrix = np.empty((len(data), len(organ_list)))
    for idx in range(0,len(data)):
        patient = data[idx - 1]
        organ_idx = 0
        for organ in organ_list:
            try:
                dose_matrix[idx, organ_idx] = patient['organData'][organ]['meanDose']
            except:
                dose_matrix[idx, organ_idx] = 0
            organ_idx += 1
    return dose_matrix

dose_matrix = gen_dose_matrix(data)
