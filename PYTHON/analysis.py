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
    disimilarity = (disimilarity + disimilarity.T)/2
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

def dose_similarity(dose_predictions, distance_metric = None):
    if distance_metric is None:
        distance_metric = mse
    n_patients = dose_predictions.shape[0]
    dists = np.zeros((n_patients, n_patients))
    for p1 in range(n_patients):
        d1 = dose_predictions[p1]
        for p2 in range(p1+1, n_patients):
            d2 = dose_predictions[p2]
            dists[p1,p2] = distance_metric(d1, d2)
    dists += dists.transpose()
    similarity = dist_to_sim(dists)
    return similarity
    
def get_tumor_organ_vectors(db):
    o_centroids, t_centroids = db.get_transformed_centroids()
    vectors = np.zeros((o_centroids.shape))
    for p in range(db.get_num_patients()):
        o_centers = o_centroids[p]
        t_centers = t_centroids[p]
        distances = np.stack([g.dists for g in db.gtvs[p] if g.volume > 0], axis = 1)
        new_vectors = np.zeros(o_centers.shape)
        for organ in range(new_vectors.shape[0]):
            nearest_tumor_arg = np.argmin(distances[organ])
            nearest_tumor = t_centers[nearest_tumor_arg]
            organ_tumor_vector = nearest_tumor - o_centers[organ]
            new_vectors[organ] = organ_tumor_vector/np.linalg.norm(organ_tumor_vector)
        vectors[p] = new_vectors
    return vectors

def tumor_cosine_similarity(p1, p2, t_o_vectors, adjacency):
    vects1 = t_o_vectors[p1]
    vects2 = t_o_vectors[p2]
    dist = []
    for organ in range(t_o_vectors.shape[1]):
        overlap = np.linalg.norm(np.dot(vects1[organ], vects2[organ]))
        dist.append(overlap)
    return np.mean(dist)

#db = PatientSet(root = 'data\\patients_v*\\',
#                use_distances = False)

from sklearn.ensemble import RandomForestRegressor
class_densities = [len(np.argwhere(db.classes == c))/len(db.classes) for c in sorted(np.unique(db.classes))]
class_args = [np.argwhere(db.classes == c).ravel() for c in sorted(np.unique(db.classes))]
generator = RandomForestRegressor(n_estimators = 10)
patients_to_generate = 200
for c in range(len(class_args)):
    args = class_args[c]
    o_centroids = db.centroids[c]
    tumor_centroids = []
    tumor_volumes = []
    tumor_distances = []
    training_organ_centroids = []
    for arg in args:
        tumorset = db.gtvs[arg]
        for tumor in tumorset:
            if tumor.volume > 0:
                tumor_centroids.append(tumor.position)
                tumor_volumes.append(tumor.volume)
                tumor_distances.append(tumor.dists)
                training_organ_centroids.append(db.centroids[arg].ravel())
    tumor_centroids = np.vstack(tumor_centroids).astype('float32')
    tumor_volumes = np.vstack(tumor_volumes).astype('float32')
    t_centroid_stats =  (tumor_centroids.mean(axis = 0), tumor_centroids.std(axis = 0 ))
    volume_stats = (tumor_volumes.mean(), tumor_volumes.std())
    y = np.hstack([np.vstack(tumor_distances), np.vstack(training_organ_centroids)])
    x = np.hstack([tumor_centroids, tumor_volumes])
    generator.fit(x,y)
    generated_tumor_distance = np.zeros((patients_to_generate, Constants.nunum_organs))
#t_o_vectors = get_tumor_organ_vectors(db)
#adjacency = TsimModel().get_adjacency_lists(db.tumor_distances)
#cosine_dist = lambda d,x,y: tumor_cosine_similarity(x,y, t_o_vectors, adjacency)
#cosine_sim = get_sim(db, cosine_dist)
#distance_sim = TJaccardModel().get_similarity(db, augment = False)
#vol_sim = dist_to_sim(get_sim(db, gtv_volume_dist))
#total_dose_sim = dist_to_sim(get_sim(db, lambda d,x,y: np.abs(db.prescribed_doses[x] - db.prescribed_doses[y])))
#
#similarity_list = [cosine_sim, distance_sim, vol_sim, total_dose_sim]
#fused_similarity = SimilarityBooster().get_similarity(db, similarity_list)
#export(db, similarity = fused_similarity)
#threshold_grid_search(db, fused_similarity, n_itters = 10)

#from sklearn.ensemble import AdaBoostClassifier, AdaBoostRegressor
#tree = AdaBoostRegressor(n_estimators = 2*len(similarity_list),
#                          learning_rate = 1,
#                          random_state = 1)
#true_similarity = dose_similarity(db.doses)
#similarity_stack = np.hstack(similarity_list)
#n_measures = len(similarity_list)
#for p in range(db.get_num_patients()):
#    




