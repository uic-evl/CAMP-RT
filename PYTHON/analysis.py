# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 10:08:52 2019

@author: Andrew
"""
from numpy.random import seed
seed(1)
from tensorflow.compat.v1 import set_random_seed
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
from copy import copy
#import metric_learn
from preprocessing import *
from Metrics import *

from SyntheticDataGenerator import *
from sklearn.manifold import TSNE, MDS
from sklearn.cluster import KMeans


def export(data_set = None,
           patient_data_file = '../data/patient_dataset.json',
           score_file = 'data/scores.csv',
           method = 'tanimoto',
#           model = None,
#           estimator = None,
#           similarity = None,
#           predicted_doses = None,
           clusterer=None):
    if data_set is None:
        data_set = PatientSet(root = 'data/patients_v*/',
                use_distances = False)
    if method == 'tsim':
        estimator = KnnEstimator()
        similarity = tsim_similarity(data_set)
        predicted_doses = tsim_prediction(data_set, similarity)
    else:
        if method != 'tanimoto':
            print('error, unknown prediction method given')
        estimator = TreeKnnEstimator()
        similarity = default_similarity(data_set)
        predicted_doses = default_rt_prediction(data_set, [similarity])

    if clusterer == 'default':
        from sklearn.cluster import KMeans
        clusterer = KMeans(n_clusters = 3)
    error = estimator.get_error(predicted_doses, data_set.doses) #a vector of errors
    print('error: ', error.mean(), '%')
    n_patients = data_set.get_num_patients()
    disimilarity = 1- np.round(similarity[:n_patients, :n_patients], 3)
    disimilarity = (disimilarity + disimilarity.T)/2
    flipped = np.ones(similarity.shape).astype('int32')
    if n_patients < similarity.shape[0]:
        print('using augmented similarities')
        where_flipped = np.where(similarity[:n_patients, n_patients:] > similarity[:n_patients, :n_patients])
        flipped[where_flipped] = -1
        similarity = np.maximum(similarity[:n_patients, :n_patients], similarity[:n_patients, n_patients:])
    similar_patients = estimator.get_matches(similarity, data_set)
    dose_pca = pca(data_set.doses)
    distance_tsne = TSNE(perplexity = 60, init = 'pca').fit_transform(data_set.tumor_distances)
    similarity_embedding = MDS(dissimilarity='precomputed', random_state = 1).fit_transform(disimilarity)
    if clusterer is not None:
        clusters = clusterer.fit_predict(disimilarity).ravel()
        clusters = (clusters - clusters.min() + 1).astype('int32')
    else:
        clusters = data_set.classes
    clusters = clusters.astype('int32') - clusters.min() + 1
    export_data = []
    for x in range(data_set.get_num_patients()):
        entry = OrderedDict()
        entry['ID'] = data_set.ids[x]
        entry['ID_internal'] = x+1
        matches = similar_patients[x].tolist()
        simsort_args = np.argsort(-similarity[x,:]).ravel()
        local_similarity = similarity[x, simsort_args]*flipped[x, simsort_args]
        entry['similarity_scores'] = (local_similarity[:len(matches)]).tolist()
        entry['similar_patients'] = matches
        entry['mean_error'] = round(error[x], 4)

        entry['cluster'] = clusters[x]
        entry['laterality'] = data_set.lateralities[x]
        entry['tumorSubsite'] = data_set.subsites[x]
        entry['total_Dose'] = int(data_set.prescribed_doses[x])
        entry['dose_pca'] = dose_pca[x, :].tolist()
        entry['distance_pca'] = distance_tsne[x, :].tolist()
        entry['similarity_embedding'] = similarity_embedding[x,:].tolist()
        entry['toxicity'] = 1 if (data_set.feeding_tubes[x] + data_set.aspiration[x]) > 0 else 0
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
        tvols = [tumor.volume for tumor in tumors]
        entry['tumorVolume'] = np.sum(tvols)
        entry['gtvp_volume'] = tumors[0].volume
        if len(tumors) > 1:
            entry['gtvn_volume'] = np.sum(tvols[1:])
        else:
            entry['gtvn_volume'] = 0
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
                best_threshold = copy(k)
                best_min_matches = copy(m)
    if print_out:
        print('Score-', round(100*best_score,2) , ': Threshold-', round(best_threshold,2) , ': Min matches-', best_min_matches)
    if get_model:
        return KnnEstimator(match_type = 'threshold',
                            match_threshold = best_threshold,
                            min_matches = best_min_matches)
    else:
        return((best_score, best_threshold, best_min_matches))

def tsim_similarity(db):
    return TsimModel().get_similarity(db, augment = True)

def tsim_prediction(db, sim = None):
    sim = tsim_similarity(db) if sim is None else sim
    return KnnEstimator().predict_doses(sim, db)

def default_similarity(db):
    discrete_dists = discretize(-db.tumor_distances)
    return augmented_sim(discrete_dists, jaccard_distance)

def default_rt_prediction(db, similarity = None):
    similarity = [default_similarity(db)] if similarity is None else similarity
    estimator = TreeKnnEstimator()
    return estimator.predict_doses(similarity, db)

