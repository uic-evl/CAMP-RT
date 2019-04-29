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
    
#model = TsimModel()

#db = PatientSet(root = 'data\\patients_v*\\',
#                class_name = None,
#                use_distances = False)
#
#distance_similarity = TsimModel().get_similarity(db)

#asymetric_lymph_similarity = get_lymph_similarity(db)
percent_diff = lambda x,y: 1 - np.abs(x-y)/np.max([x,y])
symmetric_lymph_similarity = get_sim(db, lambda d,x,y: Rankings.jaccard_distance(db.lymph_nodes[x], db.lymph_nodes[y]))
age_similarity = get_sim(db, lambda d,x,y: np.abs(d.ages[x] - d.ages[y])/d.ages.max())

subsite_similarity = get_sim(db, lambda d,x,y: 1 if d.subsites[x] == d.subsites[y] else 0)

total_dose_similarity = get_sim(db, lambda d,x,y: percent_diff(d.prescribed_doses[x], d.prescribed_doses[y]))

class_similarity = get_sim(db, lambda d,x,y: 1 if db.classes[x] == db.classes[y] else 0)

gender_similarity = get_sim(db, lambda d,x,y: 1 if d.genders[x] == d.genders[y] else 0)

n_category_similarity = get_sim(db, n_category_sim)
t_category_similarity = get_sim(db, t_category_sim)

gtv_volume_similarity = get_sim(db, gtv_volume_sim)
gtv_count_similarity = get_sim(db, gtv_count_sim)
best_score = 1000
best_k = 0
best_min_matches = 0
from sklearn.svm import SVC
estimator = SimilarityFuser(model = SVC(kernel = 'linear', probability = True))
similarity = estimator.get_similarity(db, [distance_similarity,
                                           total_dose_similarity,
                                           n_category_similarity,
                                           t_category_similarity,
                                           subsite_similarity,
                                           gtv_volume_similarity,
                                           gtv_count_similarity])
import copy
print('similarity finished')
for k in np.linspace(.3, 1, 20):
    for min_matches in range(1, 30):
        result = KnnEstimator(match_threshold = k, min_matches = min_matches).evaluate(similarity, db)
        if result.mean() < best_score:
            best_score = copy.copy(result.mean())
            best_k = copy.copy(k)
            best_min_matches = copy.copy(min_matches)
print(best_k, best_min_matches, best_score)
print(KnnEstimator(match_type = 'clusters').evaluate(similarity, db).mean())
print(KnnEstimator(match_type = 'clusters').evaluate(base_similarity, db).mean())
print(KnnEstimator(match_type = 'clusters').evaluate(distance_similarity, db).mean())
print(KnnEstimator(match_type = 'clusters').evaluate(distance_similarity*class_similarity, db).mean())



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