# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 09:19:09 2019

@author: Andrew
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
import ssim
import random

def load_patient_data(patients_file = "data\\patients_SSIM_wDoses_wDists.json"):
    with open(patients_file) as file:
        data = json.load(file)
    return data

def get_ssim_scores(entry1, entry2):
    #takes two entries (formatted as dictionaries) and returns an array of the ssim scores
    scores = np.empty((4,))
    scores[0] = ssim.compute_ssim(entry1['matrix_ssim'], entry2['matrix_ssim'])
    scores[1] = ssim.compute_ssim(entry1['matrix_ssim_dist'], entry2['matrix_ssim_dist'])
    scores[2] = ssim.compute_ssim(entry1['matrix_ssim_vol'], entry2['matrix_ssim_vol'])
    scores[3] = 1 if (entry1['laterality'] == entry2['laterality']) else 0
    return scores

def gen_dose_matrix(data):
    #takes the data from the json file an generates a num_patients by num_organs matrix of
    #mean radiation dosages. indexes are the internal id - 1
    organ_list = list(data[0]['organData'].keys())
    dose_matrix = np.empty((len(data), len(organ_list)))
    for idx in range(0,len(data)):
        patient = data[idx]
        organ_idx = 0
        for organ in organ_list:
            try:
                dose_matrix[idx, organ_idx] = patient['organData'][organ]['meanDose']
            except:
                dose_matrix[idx, organ_idx] = 0
            organ_idx += 1
    return dose_matrix

def generate_dose_estimates(ranks, doses, num_matches = 6):
    estimates = np.zeros(doses.shape)
    for patient_idx in range(0, len(doses)):
        #get index of the scores in ascending order
        top_matches = np.argsort(-ranks[patient_idx,:]) #ranks is negative so the result is in decending order
        top_matches = top_matches[0:num_matches]
        scores = ranks[patient_idx, top_matches] #scores, can be used later for figuring out scaling
        matched_dosages = doses[top_matches, :]
        for match_idx in range(0, num_matches):
            matched_dosages[match_idx, :] = scores[match_idx]*matched_dosages[match_idx, :]
        estimates[patient_idx, :] = np.mean(matched_dosages, axis = 0)/np.mean(scores)
    return estimates

def load_matrix_file(matrix_file = "latest_results\\matrices.json"):
    with open(matrix_file) as file:
        ssim_matrices = json.load(file)
    ssim_matrices = {
                        int(key): {
                            'laterality': value['laterality'],
                            'matrix_ssim': np.asarray(value['matrix_ssim']), #doses
                            'matrix_ssim_dist': np.asarray(value['matrix_ssim_dist']), #distances
                            'matrix_ssim_vol': np.asarray(value['matrix_ssim_vol']) #volume
                        }
                        for key, value in ssim_matrices.items()
                     }
    return ssim_matrices

def rank_by_ssim(ssim_matrices):
    #creates a num_pateintsxnum_patients array of ssim scores.  diagonal is zero here instead of 1 like before
    ssim_score_matrix = np.zeros((len(ssim_matrices), len(ssim_matrices)))
    for row in range(0, len(ssim_matrices)):
        for col in range(row + 1, len(ssim_matrices)): #matrix should be semetric so we take an upper-triangular matrix?
            person1 = ssim_matrices[row + 1]
            person2 = ssim_matrices[col + 1]
            scores = get_ssim_scores(person1, person2) #returns a 4x0 array so I can dot it with scalars later
            ssim_score_matrix[row, col] = np.mean(scores)
    ssim_score_matrix = ssim_score_matrix + np.transpose(ssim_score_matrix)
    return ssim_score_matrix

def rank_by_laterality(ssim_matrices):
    #creates a num_pateintsxnum_patients array of ssim scores.  diagonal is zero here instead of 1 like before
    ssim_score_matrix = np.zeros((len(ssim_matrices), len(ssim_matrices)))
    for row in range(0, len(ssim_matrices)):
        for col in range(row + 1, len(ssim_matrices)): #matrix should be semetric so we take an upper-triangular matrix?
            person1 = ssim_matrices[row + 1]
            person2 = ssim_matrices[col + 1]
            ssim_score_matrix[row, col] = 1 if person1['laterality'] == person2['laterality'] else 0
    ssim_score_matrix = ssim_score_matrix + np.transpose(ssim_score_matrix)
    return ssim_score_matrix

def rank_randomly(ssim_matrices):
    #creates a num_pateintsxnum_patients array of ssim scores.  diagonal is zero here instead of 1 like before
    ssim_score_matrix = np.zeros((len(ssim_matrices), len(ssim_matrices)))
    for row in range(0, len(ssim_matrices)):
        for col in range(row + 1, len(ssim_matrices)): #matrix should be semetric so we take an upper-triangular matrix?
            ssim_score_matrix[row, col] = random.random()
    ssim_score_matrix = ssim_score_matrix + np.transpose(ssim_score_matrix)
    return ssim_score_matrix

def save_rank(ranks):
    np.savetxt('data\\ssim_rank_matrix.csv', ranks, delimiter = ',')

def load_rank():
    ranks = np.loadtxt('data\\ssim_rank_matrix.csv', delimiter=',')
    return ranks

def generate_ssim_rank_csv():
    ssim_matrices = load_matrix_file()
    ranks = rank_by_ssim(ssim_matrices)
    save_rank(ranks)
    
def run_with_metric(rank_metric):
    data = load_patient_data()
    ssim_matrices = load_matrix_file()
    ranks= rank_metric(ssim_matrices)
    doses = gen_dose_matrix(data)
    mse_hist = []
    variance_hist = []
    for count in range(1,len(data) + 1):
        dose_estimates = generate_dose_estimates(ranks, doses, num_matches = count)
        differences = dose_estimates - doses
        variance = np.var(differences)
        mse = np.sqrt(np.mean((differences)**2))
        mse_hist.append(mse)
        variance_hist.append(variance)
    return mse_hist

def plot():
    hist1 = run_with_metric(rank_by_laterality)
    hist2 = run_with_metric(rank_by_ssim)
    hist3 = run_with_metric(rank_randomly)
    print("laterality min: ", np.min(hist1))
    print('ssim min: ', np.min(hist2))
    print('random min: ', np.min(hist3))
    plt.plot( range(0, len(hist1)), hist1, range(0, len(hist2)), hist2, range(0, len(hist3)), hist3)
    plt.show()

plot()