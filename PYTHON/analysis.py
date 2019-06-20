# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 10:08:52 2019

@author: Andrew
"""
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
from sklearn.cluster import KMeans
from NCA import NeighborhoodComponentsAnalysis
from Metrics import *



def export(data_set, patient_data_file = 'data\\patient_dataset.json', score_file = 'scores.csv',
           model = None, estimator = None, similarity = None):
    if model is None:
        model = TsimModel()
    if estimator is None:
        estimator = KnnEstimator()
    if similarity is None:
        similarity = model.get_similarity(data_set) #similarity scores
    predicted_doses = estimator.predict_doses(similarity, data_set)
    similar_patients = estimator.get_matches(similarity, data_set)
    error = estimator.get_error(predicted_doses, data_set.doses) #a vector of errors
    dose_pca = pca(data_set.doses)
    distance_pca = pca(data_set.tumor_distances)

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
        entry['distance_pca'] = distance_pca[x, :].tolist()

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
        scaled_similarity = (similarity - similarity.min())/(similarity.max() - similarity.min())
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

def get_all_features(data, num_pca_components = 10):
    num_patients = data.get_num_patients()
    tumor_volumes = np.zeros((num_patients, 2))
    tumor_count = np.array([len(gtv) for gtv in data.gtvs]).reshape(-1,1)
    for i in range(num_patients):
        gtvs = data.gtvs[i]
        gtvp_volume = gtvs[0].volume
        gtvn_volume = 0
        for gtvn in gtvs[1:]:
            gtvn_volume += gtvn.volume
        tumor_volumes[i, :] = (gtvp_volume, gtvn_volume)
    laterality = data.lateralities.reshape(num_patients, 1)
    laterality = np.vectorize(Constants.laterality_map.__getitem__)(laterality)
    subsites = data.subsites.reshape(num_patients, 1)
    subsites = np.vectorize(Constants.subsite_map.__getitem__)(subsites)
    total_doses = data.prescribed_doses.reshape(num_patients, 1)
    ages = data.ages.reshape(-1,1)
    features = np.hstack([tumor_volumes, total_doses, subsites,
                          laterality, tumor_count,
                          data.tumor_distances, data.volumes])
    pca_features = pca(features, features.shape[1])
    return (pca_features - pca_features.mean(axis = 0))/pca_features.std(axis = 0)

def get_input_tumor_features(data, num_pca_components = 10):
    num_patients = data.get_num_patients()
    tumor_volumes = np.zeros((num_patients, 2))
    tumor_count = np.array([len(gtv) for gtv in data.gtvs]).reshape(-1,1)
    for i in range(num_patients):
        gtvs = data.gtvs[i]
        gtvp_volume = gtvs[0].volume
        gtvn_volume = 0
        for gtvn in gtvs[1:]:
            gtvn_volume += gtvn.volume
        tumor_volumes[i, :] = (gtvp_volume, gtvn_volume)
    laterality = data.lateralities.reshape(num_patients, 1)
    laterality = np.vectorize(Constants.laterality_map.__getitem__)(laterality)
    subsites = data.subsites.reshape(num_patients, 1)
    subsites = np.vectorize(Constants.subsite_map.__getitem__)(subsites)
    total_doses = data.prescribed_doses.reshape(num_patients, 1)
    features = np.hstack([tumor_volumes, total_doses, 
                          subsites, tumor_count,
                          data.ajcc8.reshape(-1,1)])
    return copy.copy(features)

def get_input_organ_features(data):
    return copy.copy(data.volumes)

def get_input_distance_features(data, num_pca_components = 10):
    num_patients = data.get_num_patients()
    tumor_volumes = np.zeros((num_patients, 2))
    tumor_count = np.array([len(gtv) for gtv in data.gtvs]).reshape(-1,1)
    for i in range(num_patients):
        gtvs = data.gtvs[i]
        gtvp_volume = gtvs[0].volume
        gtvn_volume = 0
        for gtvn in gtvs[1:]:
            gtvn_volume += gtvn.volume
        tumor_volumes[i, :] = (gtvp_volume, gtvn_volume)
    laterality = data.lateralities.reshape(num_patients, 1)
    laterality = np.vectorize(Constants.laterality_map.__getitem__)(laterality)
    subsites = data.subsites.reshape(num_patients, 1)
    subsites = np.vectorize(Constants.subsite_map.__getitem__)(subsites)
    total_doses = data.prescribed_doses.reshape(num_patients, 1)
    features = np.hstack([tumor_volumes, 
                          total_doses, 
                          tumor_count, 
                          laterality,
                          data.ajcc8.reshape(-1,1),
                          data.tumor_distances])
    return copy.copy(features)

def get_input_lymph_features(data, num_pca_components = 10):
    return copy.copy(data.lymph_nodes)

def get_dose_clusters(doses):
    kmeans = KMeans(n_clusters = 5, random_state = 0)
    bad_patients = set(ErrorChecker().get_data_outliers(doses))
    good_patients = [i for i in np.arange(doses.shape[0]) if i not in bad_patients]
    kmeans_clusters = kmeans.fit_predict(doses[good_patients])
    return kmeans_clusters, good_patients

def get_nca_features(features, doses, min_components = 5, lmnn = False, k = 4, reg = .25):
    n_patients = doses.shape[0]
    if lmnn is False:
        n_components = min([min_components, features.shape[1]])
    else:
        n_components = features.shape[1]
    output = np.zeros((n_patients,n_components))
    for p in range(n_patients):
        feature_subset = np.delete(features, p, axis = 0)
        dose_subset = np.delete(doses, p , axis = 0)
        nca = get_fitted_nca(feature_subset, dose_subset,
                                          n_components = n_components,
                                          lmnn = lmnn, k = k,
                                          reg = reg)
        tf = nca.transform(features[p,:].reshape(1,-1))
#        output[p,:] = (tf - tf.min())/(tf.max() - tf.min())
        output[p,:] = tf/np.linalg.norm(tf)
        print(output[p,:])
    return output

def get_fitted_nca(features, doses, n_components = 5, lmnn = False, k = 4, reg = .25):
    if lmnn:
        nca = metric_learn.lmnn.python_LMNN(k = k, use_pca = False,
                                            regularization = reg)
    else:
        nca = NeighborhoodComponentsAnalysis(n_components = n_components,
                                         max_iter = 300,
                                         init = 'pca',
                                         random_state = 0)
    kmeans_clusters, good_clusters = get_dose_clusters(doses)
    for col in range(features.shape[1]):
        feature = features[:, col]
        if feature.std() > .0001:
            features[:, col] = (feature - feature.mean())/feature.std()
        else:
            features[:, col] = feature - feature.mean()
    sampler = SMOTE(k_neighbors  = 2)
    resampled_features, resampled_clusters = sampler.fit_resample(
            features[good_clusters, :], kmeans_clusters)
    nca.fit(resampled_features, resampled_clusters)
#    features = nca.transform(features)
    return nca

def get_nca_similarity(db, feature_type = 'tumors', min_nca_components = 4, lmnn = False, k = 4, reg = .25):
    doses = db.doses
    n_patients = doses.shape[0]
    if feature_type == 'tumors':
        input_features = get_input_tumor_features(db)
    elif feature_type in ['distance', 'distances']:
        input_features = get_input_distance_features(db)
    elif feature_type in ['lymph', 'lymph nodes']:
        input_features = get_input_lymph_features(db)
    elif feature_type in ['organ', 'organs']:
        input_features = get_input_organ_features(db)
    nca_features = get_nca_features(input_features, doses, 
                                    min_components = min_nca_components,
                                    lmnn = lmnn, k = k, reg = reg)
    similarity = np.zeros((n_patients, n_patients))
    max_similarities = set([])
    mixed_laterality = set(['R','L']) 
    for p1 in range(n_patients):
        x1 = nca_features[p1, :]
        for p2 in range(p1+1, n_patients):
            if set(db.lateralities[[p1,p2]]) == mixed_laterality:
                continue
            x2 = nca_features[p2, :]
            if np.linalg.norm(x1 - x2) < .001:
                max_similarities.add((p1,p2))
                continue
            similarity[p1, p2] = 1/np.linalg.norm(x1 - x2)
    similarity += similarity.transpose()
    similarity = .99*(similarity - similarity.min())/(similarity.max() - similarity.min())
    for pair in max_similarities:
        similarity[pair[0], pair[1]] = 1
    for i in range(n_patients):
        similarity[i,i] = 0 
    return similarity

def get_test_tumor_similarity(db):
    n_patients = db.get_num_patients()
    new_distance_similarity = np.zeros((n_patients, n_patients))
    for i in range(n_patients):
        gtvs1 = db.gtvs[i]
        for ii in range(n_patients):
            gtvs2 = db.gtvs[ii]
            new_distance_similarity[i,ii] = get_max_tumor_ssim(gtvs1, gtvs2)
    new_distance_similarity += new_distance_similarity.transpose()
    return new_distance_similarity

def centroid_based_tumor_organ_pairs(db):
    for p in range(db.get_num_patients()):
        p_centroids = db.centroids[p,:,:]
        gtv = db.gtvs[p]
        for t_ind in range(len(gtv)):
            t = gtv[t_ind]
            organ = Constants.organ_list[0]
            min_dist = np.linalg.norm(t.position - p_centroids[0])
            for o in range(1, Constants.num_organs):
                dist = np.linalg.norm(t.position - p_centroids[o])
                if dist < min_dist:
                    min_dist = dist
                    organ = Constants.organ_list[o]
            gtv[t_ind] = GTV(t.name, t.volume, t.position, t.doses, t.dists, organ)
            

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
class_similarity = get_sim(db, lambda d,x,y: 1 if db.classes[x] == db.classes[y] else 0)    

t_centroids, t_tumor_centroids = db.get_transformed_centroids()

tjaccard_similarity = TJaccardModel().get_similarity(db)

knn = threshold_grid_search(db, tjaccard_similarity, start_k = .8, n_itters = 10, get_model = True)
export(db, similarity = tjaccard_similarity, model = knn)
total_dose_distance = lambda d,p1,p2: np.abs(db.prescribed_doses[p1] - db.prescribed_doses[p2])
total_dose_similarity = dist_to_sim(get_sim(db, total_dose_distance))

estimator = SimilarityFuser()
similarity = estimator.get_similarity(db, [tjaccard_similarity*total_dose_similarity])
threshold_grid_search(db,similarity)
