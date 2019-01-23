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
import copy

#rows in the sorted matrix where the mean error is > 8 for the original-ish algorithm
outlier_rows_v1 = [0, 8, 12, 15, 22, 27, 31, 37, 42, 43, 51, 56, 59, 62, 67, 83, 86, 88, 91, 96]
#rows where error > 8 for the new-ish algorithm (ssim, weights 1, .05, 1)
outlier_rows_v2 = [0, 8, 9, 12, 22, 27, 37, 42, 43, 56, 62, 67, 87, 89, 91, 96]
#outliers > mean + 2*std for both variants above
consistent_outliers = [12, 22, 37, 88]

def load_patient_data(patients_file = "data\\patients_SSIM_wDoses_wDists.json"):
    with open(patients_file) as file:
        data = json.load(file)
    data = remove_outliers(data)
#    organs = data[0]['organData'].keys()
#    for entry in data:
#        for organ in organs:
#            if organ not in entry['organData'] and organ != 'GTVp':
#                #print('entry ', data.index(entry), ' is missing organ ', organ)
#                break
    return data

def remove_outliers(data, outliers = [37]):
    data_list = copy.copy(data)
    num_removed = 0
    for idx in outliers:
        try:
            del data_list[idx - num_removed]
            num_removed += 1
        except:
            print('list does not contain entry ', idx)
    return data_list

def get_ssim_scores(entry1, entry2):
    #takes two entries (formatted as dictionaries) and returns an array of the ssim scores
    scores = np.empty((3,))
    scores[0] = skimage.measure.compare_ssim(entry1['matrix_ssim_dist'], entry2['matrix_ssim_dist'],
          gaussian_weights = True)
    scores[1] = skimage.measure.compare_ssim(entry1['matrix_ssim_vol'], entry2['matrix_ssim_vol'],
          gaussian_weights = True)
    scores[2] = 1 if (entry1['laterality'] == entry2['laterality']) else 0
    return scores

def generate_dose_estimates(ranks, doses, min_matches = 5, min_rank = .9999999):
    estimates = np.zeros(doses.shape)
    for patient_idx in range(0, len(doses)):
        #get index of the scores in ascending order
        #use either the minimum number of matches for ones above a minimum rank
        num_matches = max([len(np.where(ranks[patient_idx,:] > min_rank)[0]), min_matches])
        top_matches = np.argsort(-ranks[patient_idx,:]) #ranks is negative so the result is in decending order
        top_matches = top_matches[0:num_matches] 
        scores = ranks[patient_idx, tuple(top_matches)] #scores, can be used later for figuring out scaling
        matched_dosages = doses[tuple(top_matches), :]
        #scale based on scores, I don't feel like this does much
        #does the ratio so that non-normalized metrics can be used
        score_ratios = scores/scores.max()
        for match_idx in range(0, num_matches):
            matched_dosages[match_idx, :] = score_ratios[match_idx]*matched_dosages[match_idx, :]
        estimates[patient_idx, :] = np.mean(matched_dosages, axis = 0)/np.mean(score_ratios)
    return estimates

def load_matrix_file(matrix_file = "latest_results\\matrices.json"):
    with open(matrix_file) as file:
        ssim_matrices = json.load(file)
    #everything is a list now to make it easier to get outlier removal to work and be consistent with the 
    #data format used before, should be a dictionary (or array which isn't a thing) if redone for faster lookup which happens a lot
    ssim_matrices = [
                        {
                            'internal_id': int(key),
                            'laterality': value['laterality'],
                            'matrix_ssim': np.asarray(value['matrix_ssim']), #doses
                            'matrix_ssim_dist': np.asarray(value['matrix_ssim_dist']), #distances
                            'matrix_ssim_vol': np.asarray(value['matrix_ssim_vol']) #volume
                        }
                        for key, value in ssim_matrices.items()
                     ]
    ssim_matrices = sorted(ssim_matrices, key = lambda x: x['internal_id'])
    ssim_matrices = remove_outliers(ssim_matrices)
    return ssim_matrices

def load_features(matrix_file = "latest_results\\features.json"):
    with open(matrix_file) as file:
        features = json.load(file)
    features = [
                    {
                        'internal_id': int(key),
                        'laterality': value['laterality'], #'L', 'R', or whatever middle is
                        'organ_distances': np.asarray(value['organ_distances']), #45x45 matrix
                        'tumor_volumes': value['tumor_volumes'], #single value
                        'total_doses': value['total_doses'], #single value
                        'tumor_distances': np.asarray(value['tumor_distances']) #len 45 array
                    }
                    for key, value in features.items()
                ]
    features = sorted(features, key = lambda x: x['internal_id'])
    features = remove_outliers(features)
    return features


def rank_by_ssim(ssim_matrices):
    #creates a num_pateintsxnum_patients array of ssim scores.  diagonal is zero here instead of 1 like before
    ssim_score_matrix = np.zeros((len(ssim_matrices), len(ssim_matrices)))
    for row in range(0, len(ssim_matrices)):
        for col in range(row + 1, len(ssim_matrices)): #matrix should be semetric so we take an upper-triangular matrix?
            person1 = ssim_matrices[row ] #ids start at one for some reason so they're offset
            person2 = ssim_matrices[col]
            scores = get_ssim_scores(person1, person2) #returns a 4x0 array so I can dot it with scalars later
            ssim_score_matrix[row, col] = np.mean(scores)
    ssim_score_matrix = ssim_score_matrix + np.transpose(ssim_score_matrix)
    return ssim_score_matrix

def rank_by_laterality(ssim_matrices):
    #creates a num_pateintsxnum_patients array of ssim scores.  diagonal is zero here instead of 1 like before
    ssim_score_matrix = np.zeros((len(ssim_matrices), len(ssim_matrices)))
    for row in range(0, len(ssim_matrices)):
        for col in range(row + 1, len(ssim_matrices)): #matrix should be semetric so we take an upper-triangular matrix?
            person1 = ssim_matrices[row]
            person2 = ssim_matrices[col]
            ssim_score_matrix[row, col] = 1 if person1['laterality'] == person2['laterality'] else 0
    ssim_score_matrix = ssim_score_matrix + np.transpose(ssim_score_matrix)
    return ssim_score_matrix

def rank_randomly(ssim_matrices):
    #creates a num_pateintsxnum_patients array of ssim scores.  diagonal is zero here instead of 1 like before
    ssim_score_matrix = np.zeros((len(ssim_matrices), len(ssim_matrices)))
    for row in range(0, len(ssim_matrices)):
        for col in range(row + 1, len(ssim_matrices)): #matrix should be semetric so we take an upper-triangular matrix?
            ssim_score_matrix[row, col] = random.random()
        #should I normalize the max score to be 1?
    ssim_score_matrix = ssim_score_matrix + np.transpose(ssim_score_matrix)
    return ssim_score_matrix

def modified_rank_by_ssim(dummy, rank_function = None, weights = np.array([1,.05,1])):
    #ranks features with a weighting
    features = load_features()
    ranks = np.zeros((len(features),len(features)))
    #organ_weights = weights[3:]
    if rank_function == 'mse':
        rank_function = lambda x, y: -1*skimage.measure.compare_mse(x,y)
    if rank_function == 'ssim' or rank_function == None:
        rank_function = lambda x, y: skimage.measure.compare_ssim(x,y)
    percent_difference = lambda x, y: 1- abs((x-y)/(x+y+.000001))
    upper_triangle_indicies = np.triu_indices(45)
    for row in range(0, len(features)):
        for col in range(row, len(features)):
            if(col == row):
                ranks[row, col] = -np.inf
                continue
            person1 = features[row]
            person2 = features[col]
            get_distances = lambda x: (
                    (x['organ_distances'] + np.diag(x['tumor_distances']))[upper_triangle_indicies]
                    ) #* organ_weights
            scores = np.array([
                    rank_function(get_distances(person1), get_distances(person2)),
                    percent_difference(person1['tumor_volumes'], person2['tumor_volumes']),
                    (1 if person1['laterality'] == person2['laterality'] else 0),
                ])
            ranks[row, col] = np.mean(scores*weights)/np.mean(weights)
    ranks = ranks + np.transpose(ranks)
    return(ranks)

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
    diff_hist = []
    for count in range(1,len(data)//2):
        dose_estimates = generate_dose_estimates(ranks, doses, min_matches = count)
        differences = np.abs(dose_estimates - doses)
        mse = np.mean(differences)#np.sqrt(np.mean(differences**2))
        mse_hist.append(mse)
        diff_hist.append(differences)
    print('min of: ', min(mse_hist), "at ", np.argmin(mse_hist) + 1)
    return((mse_hist, diff_hist[np.argmin(mse_hist)]))

def run_with_ranks(ranks, min_matches = 1):
    data = load_patient_data()
    doses = gen_dose_matrix(data)
    mse_hist = []
    diff_hist = []
    for min_rank in np.linspace(.999, .85, 100):
        dose_estimates = generate_dose_estimates(ranks, doses, min_matches = min_matches, min_rank = min_rank)
        differences = np.abs(dose_estimates - doses)
        mse = np.mean(differences)#np.sqrt(np.mean(differences**2))
        mse_hist.append((mse))
        diff_hist.append(differences)
    print('min of: ', min(mse_hist), " at ", np.argmin(mse_hist) + 1)
    #returns histor of error and difference for the lowest point
    return((mse_hist, diff_hist[np.argmin(mse_hist)]))

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
        #shouldn't this actually not do anything?  It changes it all drastically
        #dose_matrix[idx, :] = (patient['total_Dose']/np.sum(dose_matrix[idx, :]))*dose_matrix[idx, :]
    #replace missing entries with the mean value?
    for (idx, organ_idx) in bad_entries:
        dose_matrix[idx, organ_idx] = np.mean(dose_matrix[:, organ_idx])
    return dose_matrix

(mse_hist, diffs1) = run_with_metric(modified_rank_by_ssim)
(ssim_hist, diff2) = run_with_metric(rank_by_ssim)
(rand_hist, diff3) = run_with_metric(rank_randomly)
x = np.linspace(1, len(mse_hist), len(mse_hist))
plt.plot(x, mse_hist[:len(x)], x, ssim_hist[:len(x)], x, rand_hist[:len(x)])
plt.legend(['ssim_unstructured','ssim_old','random'])
plt.xlabel('Matches')
plt.ylabel('Mean error')
plt.title('Organ-Wise Mean Error vs Number of Matches With 16 Outliers Removed')
