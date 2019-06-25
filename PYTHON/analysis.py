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
import re



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
#                use_distances = False, denoise = False)
#
#distances = db.get_all_tumor_distances()
#distances = Denoiser(normalize = False, noise = .5).fit_transform(distances, lr = .0001)

o_centroids, t_centroids = db.get_transformed_centroids()
t_centroids = np.vstack(t_centroids)
from sklearn.cluster import AffinityPropagation, SpectralClustering, MeanShift
#clusterer = AffinityPropagation(max_iter = 400, damping = .96)
#clusterer = SpectralClustering(n_clusters = 6)
clusterer = MeanShift()

dist_pca = pca(distances, 3)
c_pca = pca(t_centroids, 2)

x = dist_pca

clusters = clusterer.fit_predict(dist_pca)
plt.scatter(x[:,0], x[:,1], c=clusters)
print(len(np.unique(clusters)))

tumor_sets = np.zeros((db.get_num_patients(), Constants.num_organs, 2))
for p in range(db.get_num_patients()):
    gtvs = db.gtvs[p]
    left = np.inf*np.ones((Constants.num_organs,))
    right = copy.copy(left)
    #position[0] > 0 is left side 
    for gtv in gtvs:
        if gtv.position[0] > 0:
            left = np.minimum(left, gtv.dists)
        else:
            right = np.minimum(right, gtv.dists)
    tumor_sets[p, :, 0] = left
    tumor_sets[p, :, 1] = right

def get_flip_args(organ_list = None):
    if organ_list is None:
        organ_list = Constants.organ_list
    flip_args = np.arange(len(organ_list))
    for organ in organ_list:
        for pattern in [('Rt_', 'Lt_'), ('Lt_', 'Rt_')]:
            if re.match(pattern[0], organ) is not None:
                other_organ = re.sub(pattern[0], pattern[1], organ)
                idx1 = organ_list.index(organ)
                idx2 = organ_list.index(other_organ)
                flip_args[idx1] = idx2
    return flip_args
            
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
    
    print(KnnEstimator().get_error(dose_predictions, db.doses).mean())
    return(dose_predictions)

#best_val = np.inf
#best_min_matches = 0
#best_max_error= 0
#for min_matches in range(2,20):
#    for max_error in np.arange(0, .2, 40):      
#        estimator = SimilarityFuser(min_matches = min_matches, max_error = max_error)
#        similarity = estimator.get_similarity(db, [tjaccard_similarity*total_dose_similarity])
#        error, min_matches, whatever= threshold_grid_search(db,similarity)
#        if error < best_val:
#            best_val = error
#            best_min_matches = min_matches
#            best_max_error = max_error
    
