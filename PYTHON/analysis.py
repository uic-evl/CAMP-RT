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
import cv2
from scipy.spatial.distance import directed_hausdorff
from scipy.spatial import ConvexHull, procrustes
from sklearn.cluster import KMeans, DBSCAN
from NCA import NeighborhoodComponentsAnalysis
from imblearn.over_sampling import RandomOverSampler, SMOTE, ADASYN



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

def dist_to_sim(distance):
    distance = copy.copy(distance)
    distance -= distance.min()
    distance /= distance.max()
    diagonals = ( np.arange(distance.shape[0]), np.arange(distance.shape[0]) )
    sim = 1-distance
    sim[diagonals] = 0
    return sim

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
    return Rankings.jaccard_distance(v1,v2) #1 if np.linalg.norm(v1 - v2) == 0 else 0 

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
    laterality = np.vectorize(TreeEstimator.laterality_map.__getitem__)(laterality)
    subsites = data.subsites.reshape(num_patients, 1)
    subsites = np.vectorize(TreeEstimator.subsite_map.__getitem__)(subsites)
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
            
def get_gtv_vectors(db):
    vectors = np.zeros((db.get_num_patients(), 6))
    for p in range(db.get_num_patients()):
        gtvs = db.gtvs[p]
        center = np.mean([x.position*x.volume for x in gtvs], axis = 0)/np.sum([x.volume for x in gtvs], axis = 0)
        secondary_points = np.zeros((3,))
        secondary_tumor_volume = np.sum([tumor.volume for tumor in gtvs])
        if secondary_tumor_volume > 0:
            for t in range(len(gtvs)):
                weight = gtvs[t].volume/secondary_tumor_volume
                secondary_points = secondary_points + weight*(gtvs[t].position - center)
            slope = secondary_points/np.linalg.norm(secondary_points)
        else:
            slope = np.zeros((3,))
        vectors[p] = np.hstack([center, slope])
    return vectors

def get_vector_sim(db, p1, p2):
    vectors = get_gtv_vectors(db)
    return np.dot(vectors[p1, 3:], vectors[p2, 3:])

def threshold_grid_search(db, similarity, start_k = .4, max_matches = 20, print_out = True, n_itters = 20):
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
    return((best_score, best_threshold, best_min_matches))

def get_autoencoderish_model(features):
    input_x = Input(shape=(features.shape[1],))
    encoder = Sequential([
            Dense(45, input_dim=features.shape[1], activation = 'relu'),
            Dense(100, activation = 'relu'),
            Dense(100, activation = 'relu'),
            Dense(100, activation = 'relu'),
            ])(input_x)
        
    decoder = Sequential([
            Dense(100,input_dim = 4, activation = 'relu',
                  activity_regularizer = regularizers.l2(.01)),
            Dense(45, activation = 'relu'),
            ])(encoder)
    model = Model(input_x, decoder)
    encoder_model= Model(input_x, encoder)
#    optimizer = optimizers.SGD(lr = .01, decay = 1e-12, momentum = .1)
    optimizer = optimizers.Adam()
    model.compile(loss = losses.mean_absolute_error, 
                  optimizer = optimizer)
    return(model, encoder_model)
    
def get_regression_model(features, activation = 'relu', lr = .01):
    model = Sequential([
            Dense(45, input_dim=features.shape[1], activation = activation),
            Dense(100, activation = activation),
            Dense(200, activation = activation),
            Dense(45, activation = activation)
            ])
    optimizer = optimizers.SGD(lr = lr, decay = 1e-4, momentum = 0.05)
    model.compile(loss = losses.mean_absolute_error, 
                  optimizer = optimizer)
    return(model)

def run_autoencoder(db):
    from keras.models import Sequential, Model
    from keras.layers import Dense, Activation, Input
    from keras import losses, optimizers,regularizers
    from sklearn.model_selection import LeaveOneOut
    features = get_input_distance_features(db)
    features = (features - features.mean(axis = 0))/features.std(axis = 0)
    clusters = db.classes.astype('int32')
    doses = db.doses
    
    loo = LeaveOneOut()
    loo.get_n_splits(features)
    regression_errors = []
    nn_sim = np.zeros((db.get_num_patients(),db.get_num_patients()))
    p1 = 0
    for train,test in loo.split(features, doses):
        model, encoder_model = get_autoencoderish_model(features)
        x_train = features[train]
        y_train = doses[train]
        x_test = features[test]
        y_test = doses[test]
        model.fit(x_train, y_train, epochs = 3, batch_size = 4, shuffle = True, verbose = 0)
        regression_error = model.evaluate(x_test, y_test)
        regression_errors.append(regression_error)
        print(regression_error)
        
        x_embedding = encoder_model.predict(features)
        for p2 in range(db.get_num_patients()):
            if p1 == p2:
                continue
            nn_sim[p1, p2] = 1/np.linalg.norm(x_embedding[p1] - x_embedding[p2])
        p1 += 1
    print(np.mean(regression_errors))
    nn_sim = (nn_sim - nn_sim.min(axis = 0))/(nn_sim.max(axis = 0) - nn_sim.min(axis = 0))
    
    threshold_grid_search(db, nn_sim)
    
def get_features(db, holdout = set([]) ):
    n_patients = db.get_num_patients()
    x = get_input_distance_features(db)
    normalizer = Normalizer()
    normalizer.fit( x[[ i for i in np.arange(n_patients) if i not in holdout] ] )
    x = normalizer.transform(x)
    a = []
    b = []
    y = []
    a_validate = []
    b_validate = []
    y_validate = []
    val_pairs = []
    true_error = SimilarityFuser(min_matches = 7, max_error = .20).get_true_matches(db)
    for p1 in range(n_patients):
        for p2 in range(p1 + 1, n_patients):
            loss = 0 if true_error[p1,p2] > 0 else 1
            if p1 in holdout or p2 in holdout:
                a_validate.append(x[p1])
                b_validate.append(x[p2])
                y_validate.append( loss )
                val_pairs.append((p1,p2))
            else:
                a.append(x[p1])
                b.append(x[p2])
                y.append( loss )
    a = np.array(a)
    b = np.array(b)
    y = np.array(y).ravel()
    a_validate = np.array(a_validate)
    b_validate = np.array(b_validate)
    y_validate = np.array(y_validate)
    
    args = np.arange(len(y))
    np.random.shuffle(args)
    
    return (a[args], b[args], y[args], a_validate, b_validate, y_validate, val_pairs)

def get_similarity_model(n_features, encoding_size = 25, reg = .000001):
    patient_a = layers.Input(shape = (n_features,))
    patient_b = layers.Input(shape = (n_features,))
    activation = 'selu'
    encoder = Sequential([
                Dense(50, input_dim = n_features, activation = activation,
                      activity_regularizer = regularizers.l2( reg )),
                Dense(100, activation = activation,
                      activity_regularizer = regularizers.l2( reg )),
                layers.Dropout(.5, seed = 0),
                Dense(encoding_size, activation = 'relu'),
                ])
    encoded_a = encoder(patient_a)
    encoded_b = encoder(patient_b)
    distance_layer = layers.dot([encoded_a, encoded_b], axes = 1, normalize = True)
    model = Model([patient_a, patient_b], distance_layer)
#    optimizer = optimizers.SGD(lr = .001, decay = 1e-8, momentum = .01)
    optimizer = optimizers.Adam()
    model.compile(optimizer = optimizer, loss = losses.mean_absolute_error)
    return(model)
    
def get_distance_model(n_features, encoding_size = 25, reg = 0.000001):
    patient_a = layers.Input(shape = (n_features,))
    patient_b = layers.Input(shape = (n_features,))
    activation = 'relu'
    encoder = Sequential([
                Dense(500, input_dim = n_features, activation = activation,
                      activity_regularizer = regularizers.l2( reg )),
                layers.Dropout(.1, seed = 0),
                Dense(encoding_size, activation = 'relu'),
                layers.BatchNormalization(),
                ])
    encoded_a = encoder(patient_a)
    encoded_b = encoder(patient_b)
    distance_layer = layers.Lambda(lambda x: K.expand_dims(K.mean(K.square(x[0] - x[1]),axis=-1),1),
                                   output_shape=(1,))([encoded_a, encoded_b])
    distance_activation = Activation('sigmoid')(distance_layer)
    model = Model(inputs=[patient_a, patient_b], outputs = distance_activation)
    distance_model = Model(inputs=[patient_a, patient_b], outputs = distance_layer)
#    optimizer = optimizers.SGD(lr = .01, decay = 1e-8, momentum = .01, nesterov = True)
    optimizer = optimizers.Adam(lr = .001, decay = 1e-8)
#    optimizer = optimizers.Adadelta()
#    optimizer = optimizers.RMSprop()
    model.compile(optimizer = optimizer, 
                  loss = losses.mean_squared_error)
    return(model, distance_model)

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
    
def get_transformed_tumor_centroids(gtvs, transform):
    t_centroids = np.ones((len(gtvs), 4))
    pos = 0
    for gtv in gtvs:
        t_centroids[pos, 0:3] = gtv.position
        pos += 1
    new_centroids = np.dot(transform, t_centroids.T).T
    return new_centroids

def get_tumor_emd_sim(db, p1, p2, centroids):
    def weighted_point(p):
        centroid = centroids[p]
        volume = np.array([g.volume for g in db.gtvs[p]]).reshape(-1,1)
        max_volume = volume.max()
        for v in range(len(volume)):
            if volume[v] == max_volume:
                volume[v] = 3
            else:
                volume[v] = np.sign(volume[v])
        return np.hstack([np.sqrt(volume), centroid]).astype('float32')
    point1 = weighted_point(p1)
    point2 = weighted_point(p2)
    return cv2.EMD(point1, point2, cv2.DIST_C)[0]

def get_organ_emd_sim(db, p1, p2, centroids):
    def weighted_point(p):
        centroid = centroids[p,:,:]
        volumes = db.volumes[p,:].reshape(-1,1)
        return np.hstack([np.sqrt(volumes), centroid]).astype('float32')
    point1 = weighted_point(p1)
    point2 = weighted_point(p2)
    return cv2.EMD(point1, point2, cv2.DIST_C)[0]

def undirected_hausdorff_distance(db,p1,p2,centroids):
    h1 = directed_hausdorff(centroids[p1], centroids[p2], 1)
    h2 = directed_hausdorff(centroids[p2], centroids[p1], 1)
    return max([h1[0], h2[0]])

def procrustes_distance(db,p1,p2,centroids, max_points = 1000):
    #max points is in case I want to use a different # later
    #centroids should be a list of n x 3 np arrays
    c1 = centroids[p1]
    c2 = centroids[p2]
    n_points = max([len(c1), len(c2), max_points])
    def scale_centroid(c):
        if len(c) < n_points:
            n_dummy_points = n_points - len(c)
            dummy_values = np.zeros((n_dummy_points, c.shape[1]))
            c = np.vstack([c, dummy_values])
        elif len(c) > n_points:
            c = c[:n_points, :]
        return c
    c1 = scale_centroid(c1)
    c2 = scale_centroid(c2)
    return procrustes(c1, c2)[2]

def single_convex_hull_projection(point_cloud, centroids, cuttoff = 15):
    hull = ConvexHull(point_cloud)
    def project(point, plane):
        distance = np.dot(plane, np.append(point, 1))
        displacement = distance*plane[0:3]/np.linalg.norm(plane[0:3])
        return (point + displacement, distance)
    projections = np.zeros(centroids.shape)
    for idx in range(centroids.shape[0]):
        point = centroids[idx]
        min_dist = np.inf
        for hull_idx in range(hull.equations.shape[0]):
            plane = hull.equations[hull_idx]
            projection, distance = project(point,plane)
            if np.abs(distance) < min_dist and projection[1] < cuttoff: #check if point is roughly in anerior of the head
                min_dist = distance
                projections[idx,:] = projection
    return projections
    
def convex_hull_projection(point_cloud, centroids):
    #looks at every plane in the convex hull and finds the projection
    #of the closest tumor?
    hull = ConvexHull(point_cloud)
    def project(point, plane):
        distance = np.dot(plane, np.append(point, 1))
        displacement = distance*plane[0:3]/np.linalg.norm(plane[0:3])
        return (point + displacement, distance)
    projections = []
    for hull_idx in range(hull.equations.shape[0]):
        plane = hull.equations[hull_idx]
        min_dist = np.inf
        curr_projection = np.zeros((3,))
        for idx in range(centroids.shape[0]):
            point = centroids[idx]
            projection, distance = project(point,plane)
            if distance < min_dist:
                curr_projection = projection
                min_dist = distance
        projections.append(curr_projection)
    return np.vstack(projections)

#db = PatientSet(root = 'data\\patients_v*\\',
#                use_distances = False)
#class_similarity = get_sim(db, lambda d,x,y: 1 if db.classes[x] == db.classes[y] else 0)    

reference_centroids = db.centroids.mean(axis = 0)
new_tumor_centroids = []
centered_tumor_centroids = []
new_organ_centroids = np.empty(db.centroids.shape)
tumor_projections = []
all_tumor_projections = []
for p in range(db.get_num_patients()):
    centroids = db.centroids[p,:,:]
    transform = cv2.estimateAffine3D(centroids, reference_centroids)[1]
    ref = np.hstack([centroids, np.ones((Constants.num_organs, 1))]).T
    transformed_centroids = np.dot(transform, ref)
    new_organ_centroids[p,:,:] = transformed_centroids.T
    patient_tumor_centroids = get_transformed_tumor_centroids(db.gtvs[p], transform)
    new_tumor_centroids.append(patient_tumor_centroids)
    tumor_projections.append(single_convex_hull_projection(transformed_centroids.T, patient_tumor_centroids))
    all_tumor_projections.append(convex_hull_projection(transformed_centroids.T, patient_tumor_centroids))
    max_tumor_pos = np.argmax([g.volume for g in db.gtvs[p]])
    centered_tumor_centroids.append( patient_tumor_centroids - patient_tumor_centroids[max_tumor_pos]  )

#organ_emd = lambda d,p1,p2: get_organ_emd_sim(d, p1, p2, new_organ_centroids)
#tumor_emd = lambda d,p1,p2: get_tumor_emd_sim(d,p1,p2, new_tumor_centroids)
#tumor_haus = lambda d,p1,p2: undirected_hausdorff_distance(d,p1,p2, new_tumor_centroids)
#procrustes_d = lambda d,p1,p2: procrustes_distance(d,p1,p2, new_tumor_centroids)
#organ_emd_sim = dist_to_sim(get_sim(db, organ_emd))
#tumor_emd_sim = dist_to_sim(get_sim(db, tumor_emd))
#hausdorf_sim = dist_to_sim(get_sim(db, tumor_haus))
#procrustes_sim = dist_to_sim(get_sim(db, procrustes_d))
#jaccard_tsim = TsimModel(organs = [Constants.organ_list.index(o) for o in Constants.optimal_jaccard_organs], similarity_function = Rankings.jaccard_distance).get_similarity(db)
#
#
#threshold_grid_search(db, similarity = hausdorf_sim, start_k = .6, n_itters = 8)
#threshold_grid_search(db, similarity = organ_emd_sim, start_k = .6, n_itters = 8)
#threshold_grid_search(db, similarity = tumor_emd_sim, start_k = .6, n_itters = 8)
#threshold_grid_search(db, similarity = procrustes_sim, start_k = .6, n_itters = 8)
#threshold_grid_search(db, similarity = jaccard_tsim, start_k = .6, n_itters = 8)

#p = np.array([0,1,2])
#nn_sim = np.zeros((db.get_num_patients(), db.get_num_patients()))
#nn_sim2 = np.zeros(nn_sim.shape)
#while p.min() < db.get_num_patients():
#    while p.max() >= db.get_num_patients():
#        p = p[:-1]
#    (x1, x2, y, x1_val, x2_val, y_val, val_ids) = get_features(db, holdout = set(p))
#    p = p + len(p)
#    model, distance_model = get_distance_model(x1.shape[1])
#    model.fit([x1, x2], y, 
#              epochs = 100, 
#              batch_size = 90*4, 
#              shuffle = True, 
#              verbose = 0,
#              validation_data = ([x1_val, x2_val], y_val))
#    y_pred = distance_model.predict([x1_val, x2_val])
#    output = model.predict([x1_val, x2_val])
#    print(p.max(), model.evaluate([x1_val, x2_val], y_val))
#    print((y_pred[0:5]).ravel(), y_pred.mean())
#    print(output[0:5].ravel(), output.mean()) 
#    for idx in range(len(y_pred)):
#        score = y_pred[idx]
#        (p1, p2) = val_ids[idx]
#        nn_sim[p1, p2] = score
#        nn_sim2[p1, p2] = output[idx]
#nn_sim += nn_sim.transpose()
#nn_sim2 += nn_sim2.transpose()
#nn_sim = (nn_sim - nn_sim.min())/(nn_sim.max() - nn_sim.min())
#threshold_grid_search(db, nn_sim)

#asymetric_lymph_similarity = get_lymph_similarity(db)
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
#vector_similarity = get_sim(db, get_vector_sim)
#organ_similarity = get_sim(db, gtv_organ_sim)

#ajcc_similarity = get_sim(db, 
#                          lambda d,x,y: percent_diff(d.ajcc8[x] + 1, d.ajcc8[y] + 1))

#class_similarity = get_sim(db, lambda d,x,y: 1 if db.classes[x] == db.classes[y] else 0)
#distance_similarity = TsimModel(
#        organs = [Constants.organ_list.index(o) for o in Constants.optimal_organs]).get_similarity(db)
#best_score, best_k, best_min_matches = threshold_grid_search(db, distance_similarity)
#best_score, best_k, best_min_matches = threshold_grid_search(db, distance_similarity*class_similarity)

#nca_tumor_similarity = get_nca_similarity(db, min_nca_components = 15, lmnn = True)
#nca_distance_similarity = get_nca_similarity(db, 'distances', min_nca_components = 15)
#lmnn_distance_similarity = get_nca_similarity(db, 'distances', min_nca_components = 15, lmnn=True)
#nca_lymph_similarity = get_nca_similarity(db, 'lymph')
#nca_organ_similarity = get_nca_similarity(db, 'organs')

#estimator = SimilarityFuser()
#similarity = estimator.get_similarity(db, [distance_similarity*total_dose_similarity])

#print('similarity finished')

#best_score, best_k, best_min_matches = threshold_grid_search(db, similarity)
#best_score, best_k, best_min_matches = threshold_grid_search(db, nca_distance_similarity)
#best_score, best_k, best_min_matches = threshold_grid_search(db, lmnn_distance_similarity)
           
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


