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
from collections import OrderedDict
import matplotlib.pyplot as plt

def export(data_set, patient_data_file = 'data\\patient_dataset.json', model = None, estimator = None, similarity = None):
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
        entry['total_Dose'] = data_set.prescribed_doses[x]
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
        def default(o):
            if isinstance(o, np.int32):
                return int(o)
        with open(patient_data_file, 'w+') as f:  # generate JSON
            json.dump( export_data, f, indent=4, default = default)
        print('successfully save patient data to ', patient_data_file)
        #save a labeled matrix of similarity scores for other people
    except:
        print('error exporting patient data to json')
#        try:
#            raw_scores = self.gen_score_matrix(1, classes = False)
#            score_df = pd.DataFrame(raw_scores, index = self.ids, columns = self.ids)
#            score_df.to_csv(score_file)
#            print('successfully saved similarity score matrix to ', score_file)
#        except:
#            print('error saving ssim score matrix')
    return

def get_input_features(data, num_pca_components = 6):
    num_patients = data.get_num_patients()
    pca = lambda x: Rankings.pca(x, num_pca_components)
    distances = pca(data.tumor_distances)
    lymph_nodes = pca(data.lymph_nodes)
    tumor_volumes = np.zeros((num_patients, 2))
    for i in range(num_patients):
        gtvs = data.gtvs[i]
        gtvp_volume = gtvs[0].volume
        gtvn_volume = 0
        for gtvn in gtvs[1:]:
            gtvn_volume += gtvn.volume
        tumor_volumes[i, :] = (gtvp_volume, gtvn_volume)
    laterality = data.lateralities.reshape(num_patients, 1)
    laterality = np.vectorize(TreeSimilarity.laterality_map.__getitem__)(laterality)
    subsites = data.subsites.reshape(num_patients, 1)
    subsites = np.vectorize(TreeSimilarity.subsite_map.__getitem__)(subsites)
    total_doses = data.prescribed_doses.reshape(num_patients, 1)
    clusters = data.classes.reshape(num_patients, 1)
    features = np.hstack([distances, lymph_nodes, tumor_volumes, total_doses, subsites])
    return features

def get_true_similarity(data):
    n_patients = data.get_num_patients()
    doses = data.doses
    error_matrix = np.zeros((n_patients, n_patients))
    for p1 in range(n_patients):
        for p2 in range(p1 + 1, n_patients):
            dose_difference = np.abs(doses[p1,:] - doses[p2, :])
            error_matrix[p1, p2] = np.mean(dose_difference)
    similarity_matrix = 1/error_matrix
    similarity_matrix = similarity_matrix/similarity_matrix.max()
    similarity_matrix += similarity_matrix.transpose()
    return similarity_matrix

#model = TsimModel()
#model = NodeSimilarityModel()

from sklearn.cluster import KMeans
from fastNCA import NCA

#db = PatientSet(root = 'data\\patients_v*\\',
#                class_name = None,
#                use_distances = False)
#

#distance_similarity = TsimModel().get_similarity(db)

num_patients = db.get_num_patients()
subsite_similarity = np.zeros((num_patients, num_patients))
tumor_volume_similarity = np.zeros((num_patients, num_patients))
total_dose_similarity = np.zeros((num_patients, num_patients))
num_tumors_similarity = np.zeros((num_patients, num_patients))

age_similarity = np.zeros((num_patients, num_patients))
gender_similarity = np.zeros((num_patients, num_patients))
n_category_similarity = np.zeros((num_patients, num_patients))
pathological_grade_similarity = np.zeros((num_patients, num_patients))

for p1 in range(num_patients):
    for p2 in range(num_patients):
        subsite_similarity[p1, p2] = 1 if db.subsites[p1] == db.subsites[p2] else 0
        total_dose_similarity[p1, p2] = 1 if db.prescribed_doses[p1] == db.prescribed_doses[p2] else 0
        age_similarity[p1, p2] = np.abs(db.ages[p1] - db.ages[p2])/db.ages.max()
        gender_similarity[p1, p2] = 1 if db.genders[p1] == db.genders[p2] else 0
        pathological_grade_similarity[p1,p2] = 1 if db.pathological_grades[p1] == db.pathological_grades[p2] else 0
        n_category_similarity[p1,p2] = 1 if db.n_categories[p1] == db.n_categories[p2] else 0
        
        gtvs1 = db.gtvs[p1]
        gtvs2 = db.gtvs[p2]
        num_tumors = max([len(gtvs1), len(gtvs2)])
        tumor_volumes_1 = np.zeros((num_tumors,))
        tumor_volumes_2 = np.zeros((num_tumors,))
        for tumor in range(num_tumors):
            try:
                tumor_volumes_1[tumor] = gtvs1[tumor].volume
            except: 
                pass
            try:
                tumor_volumes_2[tumor] = gtvs2[tumor].volume
            except: 
                pass
        num_tumors_similarity[p1, p2] = min([len(gtvs1), len(gtvs2)])/num_tumors
        tumor_volume_similarity[p1, p2] = np.sum(tumor_volumes_1.sum() - tumor_volumes_2.sum())
estimator = SimilarityFuser()
similarity = estimator.get_similarity(db, [distance_similarity, subsite_similarity, 
                            tumor_volume_similarity, total_dose_similarity, 
                            num_tumors_similarity, age_similarity])

result = KnnEstimator().evaluate(similarity, db)
print(result.mean())
export(db, similarity = similarity)

#clusterer = KMeans(n_clusters = 7)
#clusterer.fit(db.doses)
#clusters = clusterer.predict(db.doses)
#
#nca =NCA(dim=None)
#features = get_input_features(db)
#features = (features - features.mean(axis=0))/(features.std(axis=0))
#features = nca.fit_transform(features, clusters)
#similarity = np.zeros((db.get_num_patients(), db.get_num_patients()))
#for p in range(db.get_num_patients()):
#    x1 = features[p,:]
#    for p2 in range(p+1, db.get_num_patients()):
#        x2 = features[p2, :]
#        similarity[p,p2] = 1/np.linalg.norm(x1 - x2)
#similarity /= similarity.max()
#similarity += similarity.transpose()
#
#clusterer.fit(features)
#db.classes = clusterer.predict(features) + 1
##similarity = .99 * (similarity - similarity.max())/(similarity.max() - similarity.min())
#result = KnnEstimator().evaluate(similarity, db)
#print(result.mean())
#export(db, similarity = similarity)
#print(clusters)

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