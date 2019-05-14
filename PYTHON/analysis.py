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
    dose_pca = Rankings.pca(data_set.doses)
    distance_pca = Rankings.pca(data_set.tumor_distances)

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
        entry['total_Dose'] = float(data_set.prescribed_doses[x])
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
        entry['tumorVolume'] = max([tumor.volume for tumor in tumors])
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

def get_lymph_similarity(db, file = 'data/spatial_lymph_scores.csv'):
    lymph_df = pd.read_csv(file, index_col = 0)
    all_patients = set(lymph_df.index)
    similarity = np.zeros((db.get_num_patients(), db.get_num_patients()))
    for p1 in range( db.get_num_patients() ):
        name1 = 'Patient ' + str(db.ids[p1])
        if name1 not in all_patients:
            print(name1, 'Not in Lymph Data', db.n_categories[p1], db.therapy_type[p1], db.gtvs[p1][1].volume)
            continue
        for p2 in range(p1 + 1, db.get_num_patients()):
            name2 = 'Patient ' + str(db.ids[p2])
            if name2 not in all_patients:
                continue
            similarity[p1, p2] = lymph_df.loc[name1, name2]
    return similarity + similarity.transpose()

def get_sim(db, similarity_function):
    num_patients = db.get_num_patients()
    similarity_matrix = np.zeros((num_patients, num_patients))
    for p1 in range(num_patients):
        for p2 in range(p1 + 1, num_patients):
            similarity_matrix[p1,p2] = similarity_function(db, p1, p2)
    similarity_matrix += similarity_matrix.transpose()
    return similarity_matrix

def n_category_sim(db, p1, p2):
    n_categories = db.n_categories.astype('int32')
    n1 = n_categories[p1]
    n2 = n_categories[p2]
    normalized_difference = abs(n1 - n2)/n_categories.max()
    return 1 - normalized_difference

t_category_map = {'Tis': 0, 'Tx': 1, 'T1': 2, 'T2': 3, 'T3': 4, 'T4': 5}
def t_category_sim(db,p1,p2):
    t1 = db.t_categories[p1]
    t2 = db.t_categories[p2]
    normalized_difference = (abs(t_category_map.get(t1, 0) - t_category_map.get(t2, 0))/5)
    return 1 - normalized_difference


def gtv_volume_sim(db,p1,p2):
    gtvns1 = db.gtvs[p1]
    gtvns2 = db.gtvs[p2]
    vol1 = sorted([gtv.volume for gtv in gtvns1], key = lambda x: -x)
    vol2 = sorted([gtv.volume for gtv in gtvns2], key = lambda x: -x)
    if max([vol1, vol2]) == 0:
        return 1
    return np.abs(np.sum(vol1) - np.sum(vol2))/(np.sum(vol1) + np.sum(vol2))

def gtv_organ_sim(db,p1,p2):
    def vectorify(p):
        v = np.zeros((Constants.num_organs,))
        for gtv in db.gtvs[p]:
            if gtv.organ in Constants.organ_list:
                pos = Constants.organ_list.index(gtv.organ)
                v[pos] = 1
        return v
    v1 = vectorify(p1)
    v2 = vectorify(p2)
    return Rankings.jaccard_distance(v1,v2)

def gtv_count_sim(db,p1,p2):
    gtvs1 = db.gtvs[p1]
    gtvs2 = db.gtvs[p2]
    count1 = 0
    count2 = 0
    for gtv in gtvs1:
        if gtv.volume > 0:
            count1 += 1
    for gtv in gtvs2:
        if gtv.volume > 0:
            count2 += 1
    return min([count1, count2])/max([count1, count2])

def gtv_volume_jaccard_sim(db,p1,p2):
    vols1 = [gtv.volume for gtv in db.gtvs[p1]]
    vols2 = [gtv.volume for gtv in db.gtvs[p2]]
    vector_len = np.max([len(vols2), len(vols1)])
    volume_array1 = np.zeros((vector_len,))
    volume_array2 = np.zeros((vector_len,))
    volume_array1[0:len(vols1)] = vols1
    volume_array2[0:len(vols2)] = vols2
    return Rankings.jaccard_distance(volume_array1, volume_array2)

def get_all_features(data, num_pca_components = 10):
    pca = lambda x: Rankings.pca(x, num_pca_components)
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
    laterality = np.vectorize(TreeEstimator.laterality_map.__getitem__)(laterality)
    subsites = data.subsites.reshape(num_patients, 1)
    subsites = np.vectorize(TreeEstimator.subsite_map.__getitem__)(subsites)
    total_doses = data.prescribed_doses.reshape(num_patients, 1)
    ages = data.ages.reshape(-1,1)
    features = np.hstack([tumor_volumes, total_doses, subsites,
                          laterality, tumor_count,
                          data.tumor_distances, data.volumes])
    pca_features = Rankings.pca(features, features.shape[1])
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
    laterality = np.vectorize(TreeEstimator.laterality_map.__getitem__)(laterality)
    subsites = data.subsites.reshape(num_patients, 1)
    subsites = np.vectorize(TreeEstimator.subsite_map.__getitem__)(subsites)
    total_doses = data.prescribed_doses.reshape(num_patients, 1)
    features = np.hstack([tumor_volumes, total_doses, subsites,
                          tumor_count, laterality,
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
    laterality = np.vectorize(TreeEstimator.laterality_map.__getitem__)(laterality)
    subsites = data.subsites.reshape(num_patients, 1)
    subsites = np.vectorize(TreeEstimator.subsite_map.__getitem__)(subsites)
    total_doses = data.prescribed_doses.reshape(num_patients, 1)
    features = np.hstack([tumor_volumes, total_doses, subsites,
                          laterality, tumor_count, 
                          data.ajcc8.reshape(-1,1),
                          data.tumor_distances])
    return copy.copy(features)

def get_input_lymph_features(data, num_pca_components = 10):
    return copy.copy(data.lymph_nodes)

from sklearn.cluster import KMeans, DBSCAN
from NCA import NeighborhoodComponentsAnalysis
from imblearn.over_sampling import RandomOverSampler, SMOTE, ADASYN

def get_dose_clusters(doses):
    kmeans = KMeans(n_clusters = 5, random_state = 0)
    dbscan = DBSCAN(eps = 80, min_samples = 2)
    db_clusters = dbscan.fit_predict(doses)
    good_clusters = np.where(db_clusters >= 0)[0]
    bad_patients = set(np.where(db_clusters == -1)[0])
    kmeans_clusters = kmeans.fit_predict(doses)
    return kmeans_clusters, good_clusters

def get_nca_features(features, doses, min_nca_components = 5):
    nca = NeighborhoodComponentsAnalysis(n_components = min([min_nca_components, features.shape[1]]),
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
            features[good_clusters, :], kmeans_clusters[good_clusters])
    nca.fit(resampled_features, resampled_clusters)
    features = nca.transform(features)
    return features

def get_nca_similarity(db, feature_type = 'tumors', min_nca_components = 4):
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
                                    min_nca_components = min_nca_components)
    similarity = np.zeros((n_patients, n_patients))
    max_similarities = set([])
    for p1 in range(n_patients):
        x1 = nca_features[p1, :]
        for p2 in range(p1+1, n_patients):
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

def get_max_tumor_ssim(patient1, patient2):
    options = set()
    scores = -np.ones((len(patient1),))
    similarities = []
    for p1 in range(len(patient1)):
        t1 = patient1[p1]
        for p2 in range(len(patient2)):
            t2 = patient2[p2]
            options.add(p2)
            similarity = Rankings.local_ssim(t1.dists, t2.dists, t1.volume, t2.volume)
            similarities.append((similarity, p1, p2))
    similarities = sorted(similarities, key = lambda x: -x[0])
    for (s, p1, p2) in similarities:
        if scores[p1] == -1 and p2 in options:
            scores[p1] = s
            options.remove(p2)
    t_volumes = np.array([t.volume for t in patient1])
    max_similarity = np.mean((scores*t_volumes)/t_volumes.sum())
    return max_similarity

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

#db = PatientSet(root = 'data\\patients_v*\\',
#                use_distances = False)
#db.change_classes(class_name = 'rtward4', class_file = 'data\\clusters2.csv')

#from sklearn.ensemble import RandomForestClassifier
#from sklearn.linear_model import LogisticRegression
#from sklearn.model_selection import cross_validate, cross_val_predict
#x = get_all_features(db)
#db.change_classes()
#dose_clusters = get_dose_clusters(db.doses)[0]
#rf = RandomForestClassifier(n_estimators = 20,
#                             max_depth = 2)
#lg = LogisticRegression(solver = 'lbfgs', multi_class = 'auto', class_weight = 'balanced')
#rf_score = cross_validate(rf, x, dose_clusters, cv = 5)
#lg_score = cross_validate(lg, x, dose_clusters, cv = 5)
#
#print('random_forest', np.mean(rf_score['test_score']), rf_score['test_score'])
#print('logistic_regression', np.mean(lg_score['test_score']), lg_score['test_score'])
#
##asymetric_lymph_similarity = get_lymph_similarity(db)
#percent_diff = lambda x,y: 1 - np.abs(x-y)/np.max([x,y])
#symmetric_lymph_similarity = get_sim(db, lambda d,x,y: Rankings.jaccard_distance(db.lymph_nodes[x], db.lymph_nodes[y]))
#age_similarity = get_sim(db, lambda d,x,y: np.abs(d.ages[x] - d.ages[y])/d.ages.max())
#
#subsite_similarity = get_sim(db, lambda d,x,y: 1 if d.subsites[x] == d.subsites[y] else 0)
#
#total_dose_similarity = get_sim(db, lambda d,x,y: percent_diff(d.prescribed_doses[x], d.prescribed_doses[y]))
#
#gender_similarity = get_sim(db, lambda d,x,y: 1 if d.genders[x] == d.genders[y] else 0)
#
#n_category_similarity = get_sim(db, n_category_sim)
#t_category_similarity = get_sim(db, t_category_sim)
#
#gtv_volume_similarity = get_sim(db, gtv_volume_sim)
#gtv_count_similarity = get_sim(db, gtv_count_sim)
#
#    
db.change_classes()
organ_similarity = get_sim(db, gtv_organ_sim)
print(KnnEstimator(match_type='clusters').evaluate(organ_similarity, db).mean())

#class_similarity = get_sim(db, lambda d,x,y: 1 if db.classes[x] == db.classes[y] else 0)
#distance_similarity = TsimModel().get_similarity(db)
#
#nca_tumor_similarity = get_nca_similarity(db, min_nca_components = 15)
#nca_distance_similarity = get_nca_similarity(db, 'distances', min_nca_components = 15)
##nca_lymph_similarity = get_nca_similarity(db, 'lymph')
##nca_organ_similarity = get_nca_similarity(db, 'organs')
#best_score = 1000
#best_k = 0
#best_min_matches = 0
#estimator = SimilarityFuser()
#similarity = estimator.get_similarity(db, [nca_tumor_similarity, distance_similarity])
##similarity = nca_distance_similarity*class_similarity
#print(similarity)
#print('similarity finished')
#for k in np.linspace(.5, 1, 20):
#    for min_matches in range(1, 30):
#        result = KnnEstimator(match_threshold = k, min_matches = min_matches).evaluate(similarity, db)
#        if result.mean() < best_score:
#            best_score = copy.copy(result.mean())
#            best_k = copy.copy(k)
#            best_min_matches = copy.copy(min_matches)
#            
#export(db, similarity = similarity, estimator = KnnEstimator(match_threshold = best_k, min_matches = best_min_matches))
#print(best_k, best_min_matches, best_score)
#print(KnnEstimator(match_type = 'clusters').evaluate(similarity, db).mean())
#print(KnnEstimator(match_type = 'clusters').evaluate(nca_tumor_similarity, db).mean())
#print(KnnEstimator(match_type = 'clusters').evaluate(distance_similarity, db).mean())
#print(KnnEstimator(match_type = 'clusters').evaluate(nca_distance_similarity, db).mean())
#print('\n')
#print(KnnEstimator(match_type = 'clusters').evaluate(similarity*class_similarity, db).mean())
#print(KnnEstimator(match_type = 'clusters').evaluate(nca_tumor_similarity*class_similarity, db).mean())
#print(KnnEstimator(match_type = 'clusters').evaluate(distance_similarity*class_similarity, db).mean())
#print(KnnEstimator(match_type = 'clusters').evaluate(nca_distance_similarity*class_similarity, db).mean())


#from sklearn.mixture import GaussianMixture
#kmeans = KMeans(n_clusters = 5)
#dbscan = DBSCAN(eps = 80, min_samples = 2)
#gaussian_clusterer = GaussianMixture(n_components = 3, random_state = 4)
#gaussian_clusterer.fit(db.doses)
#clusters = dbscan.fit_predict(db.doses)
#good_clusters = np.where(clusters >= 0)[0]
#clusters = KMeans(n_clusters = 5).fit(db.doses).predict(db.doses)
#features = get_input_features(db)
#
#nca = NeighborhoodComponentsAnalysis(n_components = min([10, features.shape[1]]),
#                                     max_iter = 200,
#                                     init = 'pca')
#features = (features - features.mean(axis=0))/(features.std(axis=0))
#
#sampler = SMOTE()
#clean_features, clean_clusters = sampler.fit_resample(features[good_clusters], clusters[good_clusters])
#nca.fit(clean_features, clean_clusters)
#features = nca.transform(features)
#similarity = np.zeros((db.get_num_patients(), db.get_num_patients()))
#for p in range(db.get_num_patients()):
#    x1 = features[p,:]
#    for p2 in range(p+1, db.get_num_patients()):
#        x2 = features[p2, :]
#        similarity[p,p2] = 1/np.linalg.norm(x1 - x2)
#similarity /= similarity.max()
#similarity += similarity.transpose()
#
#db.change_classes()
#result = KnnEstimator(match_type = 'clusters').evaluate(similarity, db)
#print(result.mean())
#kmeans.fit(features[good_clusters,:])
#db.classes = kmeans.predict(features).astype('int') + 1
#export(db, similarity = similarity)

#from sklearn.svm import LinearSVC
#from sklearn.neural_network import MLPClassifier
#from sklearn.ensemble import RandomForestClassifier
#input_model = RandomForestClassifier(n_estimators = 10, min_samples_split = 4)
#model = ClassifierSimilarity(input_model)
##export(db, model = model)
#similarity = model.get_similarity(db)
#result = KnnEstimator().evaluate(similarity, db)
#print(result.mean())

#result = TreeEstimator(num_pca_components = 6,
#                       n_estimators = 45,
#                       min_samples_split = 4,
#                       max_depth = 45).evaluate(db)
#print(result.mean())

#gtvs = db.gtvs
#scores = np.zeros((db.get_num_patients(), db,get_num_patients()))
#for p1 in range(db.get_num_patients()):
#    for p2 in range(p1 + 1, db.get_num_patients()):
#        gtvs1 = gtvs[p1]
#        gtvs2 = gtvs[p2]
#