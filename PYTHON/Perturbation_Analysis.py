# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 09:19:09 2019

@author: Andrew
"""
from IPython import get_ipython
get_ipython().magic('reset -sf') #clears variables always

import numpy as np
import matplotlib.pyplot as plt
import json
import skimage.measure
import random

def load_patient_data(patients_file = "data\\patients_SSIM_wDoses_wDists.json"):
    with open(patients_file) as file:
        data = json.load(file)
    organs = data[0]['organData'].keys()
    for entry in data:
        for organ in organs:
            if organ not in entry['organData'] and organ != 'GTVp':
                #print('entry ', data.index(entry), ' is missing organ ', organ)
                break
    return data

def get_ssim_scores(entry1, entry2):
    #takes two entries (formatted as dictionaries) and returns an array of the ssim scores
    scores = np.empty((3,))
#    scores[3] = skimage.measure.compare_ssim(entry1['matrix_ssim'], entry2['matrix_ssim'],
#          gaussian_weights = True)
    scores[0] = skimage.measure.compare_ssim(entry1['matrix_ssim_dist'], entry2['matrix_ssim_dist'],
          gaussian_weights = True)
    scores[1] = skimage.measure.compare_ssim(entry1['matrix_ssim_vol'], entry2['matrix_ssim_vol'],
          gaussian_weights = True)
    scores[2] = 1 if (entry1['laterality'] == entry2['laterality']) else 0
    return scores

def generate_dose_estimates(ranks, doses, num_matches = 6):
    estimates = np.zeros(doses.shape)
    for patient_idx in range(0, len(doses)):
        #get index of the scores in ascending order
        top_matches = np.argsort(-ranks[patient_idx,:]) #ranks is negative so the result is in decending order
        top_matches = top_matches[0:num_matches]
        scores = ranks[patient_idx, tuple(top_matches)] #scores, can be used later for figuring out scaling
        matched_dosages = doses[tuple(top_matches), :]
        #scale based on scores, I don't feel like this does much
        #does this actually make sense for mse?
        score_ratios = scores/scores.max()
        for match_idx in range(0, num_matches):
            matched_dosages[match_idx, :] = score_ratios[match_idx]*matched_dosages[match_idx, :]
        estimates[patient_idx, :] = np.mean(matched_dosages, axis = 0)/np.mean(score_ratios)
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

def load_features(matrix_file = "latest_results\\features.json"):
    with open(matrix_file) as file:
        ssim_matrices = json.load(file)
    ssim_matrices = {
                        int(key): {
                            'laterality': value['laterality'], #'L', 'R', or whatever middle is
                            'organ_distances': np.asarray(value['organ_distances']), #45x45 matrix
                            'tumor_volumes': value['tumor_volumes'], #single value
                            'total_doses': value['total_doses'], #single value
                            'tumor_distances': np.asarray(value['tumor_distances']) #len 45 array
                        }
                        for key, value in ssim_matrices.items()
                     }
    return ssim_matrices


def rank_by_ssim(ssim_matrices):
    #creates a num_pateintsxnum_patients array of ssim scores.  diagonal is zero here instead of 1 like before
    ssim_score_matrix = np.zeros((len(ssim_matrices), len(ssim_matrices)))
    for row in range(0, len(ssim_matrices)):
        for col in range(row + 1, len(ssim_matrices)): #matrix should be semetric so we take an upper-triangular matrix?
            person1 = ssim_matrices[row + 1] #ids start at one for some reason so they're offset
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
    for count in range(2,len(data)//2):
        dose_estimates = generate_dose_estimates(ranks, doses, num_matches = count)
        differences = dose_estimates - doses
        mse = np.sqrt(np.mean(differences**2))
        mse_hist.append(mse)
    print('min of: ', min(mse_hist), "at ", np.argmin(mse_hist))
    return mse_hist

def run_old_plot():
    hist2 = run_with_metric(rank_by_ssim)
    hist3 = run_with_metric(rank_randomly)
    print('ssim min: ', np.min(hist2))
    print('random min: ', np.min(hist3))
    return((hist2, hist3))

def gen_dose_matrix(data):
    #takes the data from the json file an generates a num_patients by num_organs matrix of
    #mean radiation dosages. indexes are the internal id - 1
    organ_list = list(data[0]['organData'].keys())
    dose_matrix = np.empty((len(data), len(organ_list)))
    bad_entries = []
    for idx in range(0,len(data)):
        patient = data[idx]
        organ_idx = 0
        for organ in organ_list:
            try:
                dose_matrix[idx, organ_idx] = patient['organData'][organ]['meanDose']
            except:
                dose_matrix[idx, organ_idx] = np.mean(dose_matrix[:idx, organ_idx])
                bad_entries.append((idx, organ_idx))
            organ_idx += 1
        #normalize to total dose, can we assume this is a given?
        dose_matrix[idx, :] = (patient['total_Dose']/np.sum(dose_matrix[idx, :]))*dose_matrix[idx, :]
    #replace missing entries with the mean value?
    for (idx, organ_idx) in bad_entries:
        dose_matrix[idx, organ_idx] = np.mean(dose_matrix[:, organ_idx])
    return dose_matrix

def main(data, features,rank_function = None, weights = np.ones((3,))):
    #ranks features with a weighting
    ranks = np.zeros((len(features),len(features)))
    #organ_weights = weights[3:]
    if rank_function == None or rank_function == 'mse':
        rank_function = lambda x, y: -1*skimage.measure.compare_mse(x,y)
    if rank_function == 'ssim':
        rank_function = lambda x, y: skimage.measure.compare_ssim(x,y, gaussian_weights = True)
    percent_difference = lambda x, y: 1- abs((x-y)/(x+y+.000001))
    upper_triangle_indicies = np.triu_indices(45)
    for row in range(0, len(features)):
        for col in range(row, len(features)):
            if(col == row):
                ranks[row, col] = -np.inf
                continue
            person1 = features[row + 1]
            person2 = features[col + 1]
            get_distances = lambda x: (
                    (x['organ_distances'] + np.diag(x['tumor_distances']))[upper_triangle_indicies]
                    ) #* organ_weights
            #get_distances = lambda x: x['organ_distances'] + np.diag(x['tumor_distances'])
            scores = np.array([
                    rank_function(get_distances(person1), get_distances(person2)),
                    percent_difference(person1['tumor_volumes'], person2['tumor_volumes']),
                    (1 if person1['laterality'] == person2['laterality'] else 0)
                ])
            ranks[row, col] = np.mean(scores*weights[0:3])
    ranks = ranks + np.transpose(ranks)
    doses = gen_dose_matrix(data)
    mse_hist = []
    for count in range(2, len(data)//2):
        dose_estimates = generate_dose_estimates(ranks, doses, num_matches = count)
        differences = dose_estimates - doses
        mse = np.sqrt(np.mean(differences**2))
        mse_hist.append((mse))
    print('min of: ', min(mse_hist), " at ", np.argmin(mse_hist) + 2)
    return(mse_hist)

data = load_patient_data()
features = load_features()
mse_hist = main(data, features, rank_function = 'ssim', weights = np.array([1,.1,.8]))
ssim_hist = run_with_metric(rank_by_ssim)
rand_hist = run_with_metric(rank_randomly)
x = np.linspace(1, len(mse_hist), len(mse_hist))
plt.plot(x, mse_hist[:len(x)], x, ssim_hist[:len(x)], x, rand_hist[:len(x)])
plt.legend(['mse_new','ssim_old','random'])
plt.xlabel('Matches')
plt.ylabel('Average Mean-Squared Error')
plt.title('Error vs Number of Matches')
