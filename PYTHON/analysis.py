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



#model = TsimModel()
#model = NodeSimilarityModel()

from sklearn.cluster import KMeans
from sklearn.neighbors import NeighborhoodComponentsAnalysis

db = PatientSet(root = 'data\\patients_v*\\',
                class_name = None,
                use_distances = False)

clusterer = KMeans(n_clusters = 6)
clusterer.fit(db.doses)
clusters = clusterer.predict(db.doses)
db.classes = clusters+1

nca = NeighborhoodComponentsAnalysis()
new_features = nca.fit_transform(ClassifierSimilarity().get_input_features(db), clusters)
similarity = np.zeros((db.get_num_patients(), db.get_num_patients()))
for p in range(db.get_num_patients()):
    x1 = new_features[p,:]
    for p2 in range(p+1, db.get_num_patients()):
        x2 = new_features[p2, :]
        similarity[p,p2] = np.linalg.norm(x1 - x2)
result = KnnEstimators().evaluate(similarity, db)
print(result.mean())
export(db, similarity)
print(clusters)

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