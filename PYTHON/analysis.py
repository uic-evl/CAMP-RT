# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 10:08:52 2019

@author: Andrew
"""
from numpy.random import seed
seed(1)
from tensorflow import set_random_seed
set_random_seed(2)

from PatientSet import PatientSet
from ErrorChecker import ErrorChecker
from Constants import Constants
from Models import *
import numpy as np
import json
import pandas as pd
from collections import OrderedDict
import matplotlib.pyplot as plt
import copy
#import metric_learn
from preprocessing import *
from Metrics import *
import re
from sklearn.manifold import TSNE, MDS



def export(data_set, patient_data_file = 'data\\patient_dataset.json', score_file = 'scores.csv',
           model = None, estimator = None, similarity = None):
    if model is None:
        model = TJaccardModel()
    if estimator is None:
        estimator = KnnEstimator(match_type = 'clusters')
    if similarity is None:
        similarity = model.get_similarity(data_set) #similarity scores
    predicted_doses = estimator.predict_doses(similarity, data_set)
    error = estimator.get_error(predicted_doses, data_set.doses) #a vector of errors
    n_patients = data_set.get_num_patients()
    disimilarity = 1- np.round(similarity[:n_patients, :n_patients], 3)
    if n_patients < similarity.shape[0]:
        similarity = np.maximum(similarity[:n_patients, :n_patients], similarity[:n_patients, n_patients:])
    similar_patients = estimator.get_matches(similarity, data_set)
    dose_pca = pca(data_set.doses)
    distance_tsne = TSNE(perplexity = 60, init = 'pca').fit_transform(data_set.tumor_distances) #pca(data_set.tumor_distances)
    similarity_embedding = MDS(dissimilarity='precomputed', random_state = 1).fit_transform(disimilarity)
    export_data = []
    for x in range(data_set.get_num_patients()):
        entry = OrderedDict()
        entry['ID'] = data_set.ids[x]
        entry['ID_internal'] = x+1
        matches = similar_patients[x].tolist()
        local_similarity = sorted( similarity[x, :], key = lambda x: -x)
        entry['similarity_scores'] = local_similarity[:len(matches)]
        entry['similar_patients'] = matches
        entry['mean_error'] = round(error[x], 4)

        entry['cluster'] = data_set.classes[x]
        entry['laterality'] = data_set.lateralities[x]
        entry['tumorSubsite'] = data_set.subsites[x]
        entry['total_Dose'] = int(data_set.prescribed_doses[x])
        entry['dose_pca'] = dose_pca[x, :].tolist()
        entry['distance_pca'] = distance_tsne[x, :].tolist()
        entry['similarity_embedding'] = similarity_embedding[x,:].tolist()

        organ_data = OrderedDict()
        organ_centroids = data_set.centroids[x, :, :]
        for idx in range(Constants.num_organs):
            organ_entry = OrderedDict()
            organ_name = Constants.organ_list[idx]
            centroid = organ_centroids[idx, :]
            organ_entry['x'] = centroid[0]
            organ_entry['y'] = centroid[1]
            organ_entry['z'] = centroid[2]
            organ_entry['volume'] = data_set.volumes[x, idx]
            organ_entry['meanDose'] = data_set.doses[x, idx]
            organ_entry['minDose'] = data_set.min_doses[x, idx]
            organ_entry['maxDose'] = data_set.max_doses[x, idx]
            organ_entry['estimatedDose'] = predicted_doses[x, idx]
            organ_data[organ_name] = organ_entry

        tumors = data_set.gtvs[x]
        entry['tumorVolume'] = np.sum([tumor.volume for tumor in tumors])
        entry['gtvp_volume'] = tumors[0].volume
        entry['gtvn_volume'] = tumors[1].volume
        for tumor_idx in range(len(tumors)):
            tumor = tumors[tumor_idx]
            tumor_entry = OrderedDict()
            tumor_entry['x'] = tumor.position[0]
            tumor_entry['y'] = tumor.position[1]
            tumor_entry['z'] = tumor.position[2]
            tumor_entry['volume'] = tumor.volume
            tumor_entry['meanDose'] = tumor.doses[1]
            tumor_entry['minDose'] = tumor.doses[0]
            tumor_entry['maxDose'] = tumor.doses[2]
            organ_data[tumor.name] = tumor_entry
        entry['organData'] = organ_data
        export_data.append(entry)
        #save the vast dictionary of data for the front-end
    try:
        json.encoder.FLOAT_REPR = lambda o: format(o, '.2f')
        def default(o):
            if isinstance(o, np.int32):
                return int(o)
        with open(patient_data_file, 'w+') as f:  # generate JSON
            json.dump( export_data, f, indent=4, default = default)
        print('successfully save patient data to ', patient_data_file)
        #save a labeled matrix of similarity scores for other people
    except:
        print('error exporting patient data to json')
    try:
        scaled_similarity = minmax_scale(similarity)
        for i in range(scaled_similarity.shape[0]):
            scaled_similarity[i,i] = 1
        score_df = pd.DataFrame(scaled_similarity, index = data_set.ids, columns = data_set.ids)
        score_df.to_csv(score_file)
        print('successfully saved similarity score matrix to ', score_file)
    except:
        print('error saving ssim score matrix')
    return

def threshold_grid_search(db, similarity, start_k = .4, max_matches = 20, 
                          print_out = True, n_itters = 20, get_model = False):
    best_score = 100 #this is percent error at time of writing this
    best_threshold = 0
    best_min_matches = 0
    for k in np.linspace(start_k, 1, n_itters):
        for m in range(1, max_matches):
            result = KnnEstimator(match_threshold = k, min_matches = m).evaluate(similarity, db)
            if result.mean() < best_score:
                best_score = result.mean()
                best_threshold = copy.copy(k)
                best_min_matches = copy.copy(m)
    if print_out:
        print('Score-', round(100*best_score,2) , ': Threshold-', round(best_threshold,2) , ': Min matches-', best_min_matches)
    if get_model:
        return KnnEstimator(match_type = 'threshold',
                            match_threshold = best_threshold,
                            min_matches = best_min_matches)
    else:
        return((best_score, best_threshold, best_min_matches))       

def organ_selection(organ_list, db, similarity_function = None, 
                    use_classes = False):
    def tsim(x):
        model = TsimModel(organs = [Constants.organ_list.index(o) for o in x],
                                   similarity_function = similarity_function,
                                   use_classes = use_classes)
        return model.get_similarity(db)
    distance_similarity = tsim(organ_list)
    baseline = threshold_grid_search(db, distance_similarity)[0]
    optimal = (organ_list, baseline)
    bad_organs = []
    best_score = 100
    for organ in organ_list:
        organ_subset = copy.copy(organ_list)
        organ_subset.remove(organ)
        distance_subset_sim = tsim(organ_subset)
        best_score, best_threshold, best_min_matches = threshold_grid_search(db, distance_subset_sim, print_out = False)
        if best_score < baseline:
            bad_organs.append((organ, best_score, best_threshold, best_min_matches))
            if best_score < optimal[1]:
                optimal = (organ_subset, best_score)
                print(set(Constants.organ_list) - set(optimal[0]), best_score)
    return optimal

def optimal_organ_search(db, similarity_function = None, use_classes = False):
    optimal_organs = []
    organ_set = Constants.organ_list
    best_score = None
    while True:
        optimal_organs, best = organ_selection(organ_set, db, 
                                               similarity_function = similarity_function, 
                                               use_classes = use_classes)
        if len(optimal_organs) == len(organ_set):
            break
        organ_set = optimal_organs
        best_score = best
    return optimal_organs, best_score

#db = PatientSet(root = 'data\\patients_v*\\',
#                use_distances = False)
distance_sim = TJaccardModel().get_similarity(db, augment = True)
export(db, similarity = distance_sim)
result = KnnEstimator().evaluate(distance_sim, db)
print(result.mean())
#print(result[db.tumorcount_patients()].mean())
#print(result[db.tumorcount_patients(4)].mean())
#print(result[np.argwhere(db.classes < 3)].mean())
    

#distances = db.get_all_tumor_distances()
#distances = Denoiser(normalize = False, noise = .5).fit_transform(distances, lr = .0001)

#o_centroids, t_centroids = db.get_transformed_centroids()
#t_centroids = np.vstack(t_centroids)

#from sklearn.cluster import AffinityPropagation, SpectralClustering, MeanShift
#
#clusterer = AffinityPropagation(max_iter = 400, damping = .96)
#clusterer = SpectralClustering(n_clusters = 6)
#clusterer = MeanShift()

#dist_pca = pca(distances, 3)
#c_pca = pca(t_centroids, 2)
#dist_tsne = TSNE( perplexity = 100).fit_transform(db.tumor_distances)
#dist_mds = MDS().fit_transform(db.tumor_distances)
#x = dist_mds
#
#clusters = clusterer.fit_predict(x)
#plt.scatter(x[:,0], x[:,1], c=clusters)
#print(len(np.unique(clusters)))

#tumor_sets = np.zeros((db.get_num_patients(), Constants.num_organs, 2))
#for p in range(db.get_num_patients()):
#    gtvs = db.gtvs[p]
#    left = np.inf*np.ones((Constants.num_organs,))
#    right = copy.copy(left)
#    #position[0] > 0 is left side 
#    for gtv in gtvs:
#        if gtv.position[0] > 0:
#            left = np.minimum(left, gtv.dists)
#        else:
#            right = np.minimum(right, gtv.dists)
#    tumor_sets[p, :, 0] = left
#    tumor_sets[p, :, 1] = right

            
def tsim(d1, d2, adjacency):
    scores = []
    if d2.sum() == np.inf:
        return 0
    for organ_set in adjacency:
        scores.append(jaccard_distance(d1[organ_set], d2[organ_set]))
    return np.mean(scores)

def symmetric_similarity(db):
    flip_args = get_flip_args()
    adjacency = TJaccardModel().get_adjacency_lists(db.organ_distances, np.arange(Constants.num_organs))
    normal_distances = db.tumor_distances
    flipped_distances = db.tumor_distances[:, flip_args]
    flipped_doses = db.doses[:, flip_args]
    dose_predictions = np.zeros((db.get_num_patients(), Constants.num_organs))
    for p1 in range(db.get_num_patients()):
        matches = []
        for p2 in range(0, db.get_num_patients()):
            if p1 == p2:
                continue
            base_similarity = tsim(normal_distances[p1], normal_distances[p2], adjacency)
            flipped_similarity = tsim(normal_distances[p1], flipped_distances[p2], adjacency)
            if base_similarity > flipped_similarity:
                match = (base_similarity, db.doses[p2])
            else:
                match = (flipped_similarity, flipped_doses[p2])
            matches.append(match)
        matches = sorted(matches, key = lambda x: -x[0])
        n_matches = int(round(np.sqrt( len( np.where(db.classes == db.classes[p1])[0] ) )))
        prediction = np.array([x[1] for x in matches[0:n_matches]])
        weights = np.array([x[0] for x in matches[0:n_matches]]).reshape(-1,1)
        if weights.mean() <= 0:
            print(weights, p1, [x[0] for x in matches])
        dose_predictions[p1,:] = np.mean(prediction*weights, axis = 0)/np.mean(weights)
    
    return(dose_predictions)