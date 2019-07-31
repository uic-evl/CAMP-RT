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
from NCA import NeighborhoodComponentsAnalysis
import re
from SyntheticDataGenerator import *
from sklearn.manifold import TSNE, MDS
from sklearn.cluster import KMeans
from sklearn.preprocessing import KBinsDiscretizer

from sklearn.naive_bayes import ComplementNB, BernoulliNB, GaussianNB
from sklearn.cluster import KMeans
from sklearn.preprocessing import OneHotEncoder, quantile_transform
from sklearn.model_selection import cross_validate, cross_val_predict, LeaveOneOut
from sklearn.metrics import accuracy_score, recall_score, roc_auc_score, roc_curve
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier, VotingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.feature_selection import mutual_info_classif, f_classif, SelectKBest
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC

def export(data_set, patient_data_file = 'data\\patient_dataset.json', score_file = 'scores.csv',
           model = None, estimator = None, similarity = None, predicted_doses = None, clusterer=None):
    if model is None:
        model = TJaccardModel()
    if similarity is None:
        similarity = model.get_similarity(data_set) #similarity scores
    if estimator is None:
        if isinstance(similarity, list):
            estimator = TreeKnnEstimator()
        else:
            estimator = KnnEstimator(match_type = 'clusters')
    if predicted_doses is None:
        predicted_doses = estimator.predict_doses(similarity, data_set)
        if isinstance(similarity, list):
            similarity = dose_similarity(predicted_doses)
    if clusterer == 'default':
        from sklearn.cluster import KMeans
        clusterer = KMeans(n_clusters = 3)
    error = estimator.get_error(predicted_doses, data_set.doses) #a vector of errors
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
        organ_subset = copy(organ_list)
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

def classify_by_clusters(db, clusterer=None, inplace = True):
    if clusterer is None:
        clusterer = KMeans(n_clusters = 5)
    outliers = ErrorChecker().get_data_outliers(db.doses)
    innlier_doses = np.delete(db.doses, list(outliers), axis = 0)
    clusterer.fit(innlier_doses)
    clusters = clusterer.predict(db.doses)
    clusters = (clusters - clusters.min() + 1).ravel().astype('int32')
    old_classes = np.copy(db.classes)
    if inplace:
        db.classes = clusters
    return (clusters, old_classes)

def get_bayes_features(db, num_bins = 5):
    features_to_use = [db.ajcc8, db.prescribed_doses, db.tumor_distances,
                       np.array([np.sum([g.volume for g in gtv]) for gtv in db.gtvs]),
                       db.volumes.sum(axis = 1)]
    discretizer = KBinsDiscretizer(n_bins = num_bins,
                                   encode = 'ordinal',
                                   strategy = 'kmeans')
    formated_features = []
    for feature in features_to_use:
        if len(feature.shape) == 1:
            feature = feature.reshape(-1,1)
        if len(np.unique(feature)) > num_bins:
            feature = discretizer.fit_transform(feature)
        feature = feature - feature.min()
        formated_features.append(feature)
    return(np.hstack(formated_features))

def select_by_info(xtrain, ytrain, xtest, threshold = 0):
    mutual_info = mutual_info_classif(xtrain, ytrain)
    good_features = np.argwhere(mutual_info > threshold).ravel()
    if xtest.ndim >1:
        xout = xtest[:,good_features]
    else:
        xout = xtest[good_features]
    return xtrain[:,good_features], xout

def select_by_fscore(xtrain, ytrain, xtest, threshold = .33):
    if threshold >= 1:
        return xtrain, xtest
    pvals = f_classif(xtrain, ytrain)[1]
    good_features = np.argwhere(pvals < threshold).ravel()
    if len(good_features) <= 1:
        return select_by_fscore(xtrain, ytrain, xtest, threshold + .1)
    if xtest.ndim >1:
        xout = xtest[:,good_features]
    else:
        xout = xtest[good_features]
    return xtrain[:,good_features], xout

def recall_based_model(x, y, model, score_threshold = .99, selector = None):
    loo = LeaveOneOut()
    loo.get_n_splits(x)
    y_out = np.zeros(y.shape)
    if len(np.argwhere(y > 0)) <= 1:
        return y_out
    for train_index, test_index in loo.split(x):
        xtrain, ytrain = x[train_index], y[train_index]
        if len(set(ytrain)) <= 1:
            y_out[test_index] = ytrain[0]
            continue
        xtest = x[test_index]
        if selector is not None:
            if selector == 'info':
                xtrain, xtest = select_by_info(xtrain, ytrain, xtest)
            elif selector == 'correlation':
                xtrain, xtest = select_by_fscore(xtrain, ytrain, xtest)
                print(xtest.shape)
            else:
                fitted = selector.fit_transform(xtrain, ytrain)
                if fitted.ndim > 1 and fitted.shape[1] > 1:
                    xtrain = fitted
                    xtest = selector.transform(xtest)
        model.fit(xtrain, ytrain)
        yfit_pred = model.predict_proba(xtrain)
        sorted_scores = sorted(yfit_pred[:, 1], key = lambda x: -x)
        threshold_i = 0
        y_pred = yfit_pred[:,1] >= sorted_scores[threshold_i]
        while recall_score(y[train_index], y_pred) < score_threshold:
            threshold_i = threshold_i + 1
            y_pred = yfit_pred[:,1] >= sorted_scores[threshold_i]
        y_out[test_index] = model.predict_proba(xtest)[:,1] >= sorted_scores[threshold_i]
    return y_out

def ensemble_recall_based_model(x, y, model_list, score_threshold = 1, selector = None):
    if not isinstance(model_list, list):
        return recall_based_model(x, y, model_list, score_threshold, selector)
    labels = np.zeros(y.shape)
    for model in model_list:
        labels = labels + recall_based_model(x, y, model,
                                             score_threshold,
                                             selector = selector)
    return labels == len(model_list)

def nca_cv(x, y, n_components = 15, quantile = False):
    nca = NeighborhoodComponentsAnalysis(n_components = n_components, max_iter= 300)
    loo = LeaveOneOut()
    loo.get_n_splits(x)
    nca_dist = np.zeros((len(y), len(y)))
    for train_index, test_index in loo.split(x):
        nca.fit(x[train_index], y[train_index])
        xfit = nca.transform(x)
        for p in range(len(y)):
            nca_dist[test_index, p] = np.linalg.norm(xfit[test_index] - xfit[p])
    nca_sim = dist_to_sim(nca_dist)
    if quantile:
        nca_sim = quantile_transform(nca_sim, axis = 1)
    return nca_sim

def get_model_auc(x, y, model):
    ypred = cross_val_predict(model, x, y, cv = LeaveOneOut(), method = 'predict_proba')
    ypred = ypred[:,1]
    roc_score = roc_auc_score(y, ypred)
    fpr, tpr, thresholds = roc_curve(y, ypred)
    plt.plot(fpr, tpr)
    return fpr, tpr, thresholds, roc_score

def cv_feature_weights(x, y, scale = 1):
    n_patients = len(y)
    weights = np.zeros(x.shape)
    for p in range(n_patients):
        xtrain = np.delete(x, p, axis = 0)
        ytrain = np.delete(y, p, axis = 0)
        weights[p,:] = mutual_info_classif(xtrain, ytrain)
        weights[p,:] = minmax_scale(weights[p,:]**scale)
    return weights

def cv_feature_selection(x, y, percentile = 80, k = 40):
    n_patients = len(y)
    x_out = []
    k = np.min([int(np.ceil(percentile*x.shape[1]/100)), k])
    print(k)
    selector = SelectKBest(k = k)
    for p in range(n_patients):
        xtrain = np.delete(x, p, axis = 0)
        ytrain = np.delete(y, p, axis = 0)
        selector.fit(xtrain, ytrain)
        xfit = selector.transform(x[p,:].reshape(1,-1)).ravel()
        x_out.append( xfit )
    return np.vstack(x_out)

def feature_matrix(db):
    discrete_dists = discretize(-db.tumor_distances, n_bins = 15, strategy='uniform')
    t_volumes = np.array([np.sum([g.volume for g in gtvs]) for gtvs in db.gtvs]).reshape(-1,1)
    discrete_volumes = discretize(t_volumes, n_bins = 15, strategy='uniform')
    x = np.hstack([
        discrete_dists,
        discrete_volumes,
        db.prescribed_doses.reshape(-1,1),
        db.dose_fractions.reshape(-1,1),
        db.has_gtvp.reshape(-1,1),
        OneHotEncoder(sparse = False).fit_transform(db.lateralities.reshape(-1,1)),
        OneHotEncoder(sparse = False).fit_transform(db.subsites.reshape(-1,1)),
               ])
    return x

def discretize(x, n_bins = 9, encode = 'ordinal', strategy = 'kmeans'):
    discretizer = KBinsDiscretizer(n_bins = n_bins,
                                   encode = encode,
                                   strategy = strategy)
    return discretizer.fit_transform(x)

def cluster_with_model(db, similarity, toxicity, model, selector = None):
    cluster_stats = ClusterStats()
    features = feature_matrix(db)
    features = cv_feature_selection(features, toxicity)
    labels = ensemble_recall_based_model(features, toxicity,
                                         model, selector = selector)
    if isinstance(model, list):
        print(recall_score(toxicity, labels))
    else:
        print(get_model_auc(features, toxicity, model)[-1])
    print(cluster_stats.get_contingency_table(labels, toxicity))
    print(cluster_stats.fisher_exact_test(labels, toxicity))

    predicted_doses = TreeKnnEstimator().predict_doses([similarity], db)
    cluster_results = cluster_stats.get_optimal_clustering(predicted_doses, toxicity,
                              patient_subset = np.argwhere(labels).ravel(), subset = True)
    print(cluster_results[1])
    print(KnnEstimator().get_error(predicted_doses, db.doses).mean())
    print(cluster_stats.get_contingency_table(cluster_results[0], toxicity))

    return cluster_results[0], labels

def get_labels(x, y, models,depth,
               selector = None,
               min_size = 10,
               min_ratio = .1):
    x_subset = cv_feature_selection(x, y)
    labels = ensemble_recall_based_model(x_subset, y,
                                         models,
                                         selector = selector).astype('int32')
    positives = np.argwhere(labels > 0).ravel()
    negatives = np.argwhere(labels < .0001).ravel()
    n_pos, n_neg = len(positives), len(negatives)
    if depth < 4 and min([n_pos, n_neg]) > min_ratio * len(y):
        if n_pos > min_size:
            labels[positives] += get_labels(x[positives], y[positives],
                  models, depth+1, selector, min_size).astype('int32')
        if n_neg > min_size:
            labels[negatives] += get_labels(x[negatives], y[negatives],
                  models, depth+1, selector, min_size).astype('int32')
    print(labels)
    return labels

def get_rankings(x, y, tiers, models):
    if isinstance(models, list):
        estimators = [('estimator' + str(k), models[k]) for k in range(len(models))]
        models = VotingClassifier(estimators, voting = 'soft')
    classes = np.unique(tiers)
    loo = LeaveOneOut()
    x = cv_feature_selection(x, y)
    scores = cross_val_predict(models, x, y, cv = loo, method='predict_proba')[:,1]
    for c in classes:
        class_args = np.argwhere(tiers == c).ravel()
        scores[class_args] += c
    return minmax_scale(scores)

def tiered_model(db, toxicity, models, selector = None):
    cluster_stats = ClusterStats()
    features = feature_matrix(db)
    new_classes = get_labels(features, toxicity, models, 0,
                             selector = selector)

    new_classes = discretize(new_classes.reshape(-1,1),
                             n_bins = len(set(new_classes)),
                             strategy = 'uniform').ravel()
    rankings = get_rankings(features, toxicity, new_classes, models)
    print(cluster_stats.get_contingency_table(new_classes, toxicity))
    print(cluster_stats.fisher_exact_test(new_classes, toxicity))
    print(roc_auc_score(toxicity, rankings))
    return new_classes, rankings


db = PatientSet(root = 'data\\patients_v*\\',
                use_distances = False, use_clean_subset=False,
                denoise = False, outliers = Constants.all_bad_entries)

dists = TsimModel(use_classes = True).get_similarity(db)
export(db, similarity = dists, estimator = KnnEstimator(match_type = 'clusters'), patient_data_file = 'paper_data_file_165.json')





