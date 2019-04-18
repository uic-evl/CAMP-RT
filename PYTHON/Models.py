# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 10:15:21 2019

@author: Andrew
"""
import numpy as np
from Constants import Constants

class Rankings():
    #ranking functions that generate a score, takes in pateint objects
    def pca(points, n_components = 2):
        points = points - np.mean(points, axis = 0)
        cov = np.cov(points, rowvar = False)
        ev, eig = np.linalg.eigh(cov)
        args = np.argsort(ev)[::-1]
        ev = ev[args[0:n_components]]
        eig = eig[:, args[0:n_components]]
        principle_components = np.dot(points, eig)
        return(principle_components)

class KnnEstimator():
    def __init__(self):
        return

    def predict_doses(self, similarity_matrix, data):
        dose_matrix = data.doses
        clusters = data.classes
        predicted_doses = np.zeros(dose_matrix.shape)
        for p in range( dose_matrix.shape[0] ):
            num_matches = self.get_num_matches(p, similarity_matrix, clusters)
            if num_matches == 0:
                print(p, 'no matches')
            scores = similarity_matrix[p, :]
            args = np.argsort(-scores)
            args = args[0 : num_matches]
            predicted_doses[p, :] = self.get_prediction(data, scores, args, p)
        return(predicted_doses)
    
    def get_prediction(self, data, scores, args, patient):
        dose_matrix = data.doses
        matched_scores = scores[args].reshape(len(args), 1)
        matched_doses = dose_matrix[args, :]
        if matched_scores.mean() > 0:
                predicted_doses = np.mean(matched_scores*matched_doses, axis = 0)/matched_scores.mean()
        else:
            print(patient, 'doesnt have any matches')
            predicted_doses = dose_matrix.mean(axis=0) #return average if bad
        return predicted_doses

    def get_matches(self, similarity_matrix, data):
        dose_matrix = data.doses
        clusters = data.classes
        #should return a list of matched patients
        matches = []
        for p in range( dose_matrix.shape[0] ):
            num_matches = self.get_num_matches(p, similarity_matrix, clusters)
            scores = similarity_matrix[p, :]
            args = np.argsort(-scores)
            args = args[0 : num_matches] + 1
            matches.append(args)
        return(matches)

    def get_num_matches(self, p, similarity, clusters):
        #for later better use probs
        num_cluster_values = len(np.where(clusters == clusters[p])[0])
        num_matches = np.max([np.sqrt(num_cluster_values), 3])
#        max_val = similarity[p].max()
#        num_matches = len( np.where(similarity[p] == max_val)[0] )
        return int(num_matches + 1)

    def get_error(self, predicted_doses, dose_matrix):
        differences = np.abs(predicted_doses - dose_matrix)
        percent_error = np.sum(differences, axis = 1)/np.sum(dose_matrix, axis = 1)
        return percent_error

    def evaluate(self, similarity_matrix, data):
        predicted_doses = self.predict_doses(similarity_matrix, data)
        percent_error = self.get_error(predicted_doses, data.doses)
        return(percent_error)
    
class KnnTreeEstimator(KnnEstimator):
    
    def __init__(self, num_pca_components = 10, n_estimators = 50, min_samples_split = 4, max_depth = None):
        from sklearn.ensemble import RandomForestRegressor
        self.model = RandomForestRegressor(n_estimators = n_estimators)
        self.num_pca_components = num_pca_components
        
    
    def get_prediction(self, data, scores, args, p):
        dose_matrix = data.doses
#        matched_scores = scores[args].reshape(len(args), 1
        matched_doses = dose_matrix[args, :]
        features = self.get_features(data)
        self.model.fit(features[args], matched_doses)
        return self.model.predict(features[p, :].reshape(1, features.shape[1]))
        
    def get_features(self, data):
        num_patients = data.get_num_patients()
        pca = lambda x: Rankings.pca(x, self.num_pca_components)
        lymph_nodes = pca(data.lymph_nodes)
        distances = pca(data.tumor_distances)
        tumor_volumes = np.zeros((num_patients, 2))
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
        features = np.hstack([distances, lymph_nodes, tumor_volumes, total_doses, laterality, subsites])
        return features
    
    def get_num_matches(self, p, similarity, clusters):
        #for later better use probs
        num_cluster_values = len(np.where(clusters == clusters[p])[0])
        num_matches = np.max(np.sqrt([2*num_cluster_values, 3]))
        return int(num_matches)
        
class TreeSimilarity():
    
    subsite_map = {'BOT': 0, 'GPS': 1, 'NOS': 2, 'Soft palate': 3, 'Tonsil': 4}
    laterality_map = {'Bilateral': 0, 'L': 1, 'R': 2}
    
    def __init__(self, num_pca_components = 10, n_estimators =10, min_samples_split = 4, max_depth = None):
        from sklearn.ensemble import RandomForestRegressor
        self.model = RandomForestRegressor(min_samples_split = min_samples_split, 
                                           n_estimators = n_estimators,
                                           max_depth = max_depth)
        self.num_pca_components = num_pca_components
        
    def get_similarity(self, data):
        true_similarity = self.get_true_similarity(data)
        similarity = np.zeros(true_similarity.shape)
        x = self.get_input_features(data)
        for p in range(data.get_num_patients()):
            y_train = np.delete(true_similarity, p, axis = 0)
            x_train = np.delete(x, p, axis = 0)
            self.model.fit(x_train, y_train)
            similarity[p, :] = self.model.predict(x[p])
            similarity[p, p] = 0
        return similarity
            
    
    def get_true_similarity(self, data):
        n_patients = data.get_num_patients()
        doses = data.doses
        error_matrix = np.zeros((n_patients, n_patients))
        for p1 in range(n_patients):
            for p2 in range(p1 + 1, n_patients):
                dose_difference = doses[p1,:] - doses[p2, :]
                error_matrix[p1, p2] = np.mean(dose_difference)
        similarity_matrix = 1 - (error_matrix - error_matrix.max())/(error_matrix.max() - error_matrix.min())
        similarity_matrix += similarity_matrix.transpose()
        return similarity_matrix
    
    def get_input_features(self, data):
        num_patients = data.get_num_patients()
        pca = lambda x: Rankings.pca(x, self.num_pca_components)
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
        features = np.hstack([distances, lymph_nodes, tumor_volumes, total_doses, laterality, subsites])
        return features
        
class TsimModel():
    def __init__(self, max_distance = 80, patients = None, organs = None):
        self.max_distance = max_distance
        #for subsetting the data later
        self.patients = patients
        self.organs = organs

    def get_adjacency_lists(self, organ_distance_matrix, organs):
        #this code is much simpler than expected
        organ_distances = organ_distance_matrix[organs][:, organs]
        adjacency_lists = []
        for row in organ_distances:
            adjacent_args = np.argwhere(row < self.max_distance)
            adjacency_lists.append(adjacent_args.ravel())
        return adjacency_lists

    def get_similarity(self, data):
        #data is assumed to be a patientset object for now
        if self.patients is None:
            patients = range(len(data.ids))
        else:
            patients = self.patients
        if self.organs is None:
            organs = range(Constants.num_organs)
        else:
            organs = self.organs
        index = np.ix_(patients, organs)
        adjacency = self.get_adjacency_lists(data.organ_distances, organs)
        distances = data.tumor_distances[index]
        volumes = data.volumes[index]
        clusters = data.classes[patients]
        scores = self.similarity(adjacency, distances, volumes, clusters)
        return scores

    def similarity(self, adjacency, distances, volumes, clusters, similarity_function = None):
        if similarity_function is None:
            similarity_function = self.local_ssim
        num_patients, num_organs = distances.shape
        score_matrix = np.zeros((num_patients, num_patients))
        for patient1 in range(0, num_patients - 1):
            for patient2 in range(patient1 + 1, num_patients):
#                if clusters[patient1] != clusters[patient2]:
#                    continue
                scores = []
                for organ in range(num_organs):
                    adjacent_args = adjacency[organ]
                    if len(adjacent_args) < 1:
                        continue
                    d1 = distances[patient1, adjacent_args]
                    d2 = distances[patient2, adjacent_args]
                    v1 = volumes[patient1, adjacent_args]
                    v2 = volumes[patient2, adjacent_args]
                    scores.append( similarity_function(d1,d2,v1,v2) )
                score_matrix[patient1, patient2] = np.mean(scores)
        score_matrix += np.transpose(score_matrix)
        #scale to between 0 and .99
        #score_matrix = .99*(score_matrix - score_matrix.min())/(score_matrix.max() - score_matrix.min())
        return score_matrix

    def local_ssim(self, x,y,v = None, w = None):
        c1 = .000001
        c2  = .000001
        mean_x = np.mean(x)
        mean_y = np.mean(y)
        covariance = np.cov(x,y)
        numerator = (2*mean_x*mean_y + c1) * (covariance[0,1] + covariance[1,0] + c2)
        denominator = (mean_x**2 + mean_y**2 + c1)*(np.var(x) + np.var(y) + c2)
        if v is not None and 2 is not None:
            mean_v = np.mean(v)
            mean_w = np.mean(w)
            numerator *= (2*mean_v*mean_w + c1)
            denominator *= (mean_v**2 + mean_w**2 + c1)
        if denominator > 0:
            return numerator/denominator
        else:
            print('error, zero denomiator in ssim function')
            return 0

class NodeSimilarityModel():

    def __init__(self):
        pass

    def get_similarity(self, db):
        node_matrix = db.lymph_nodes
        subsites = db.subsites
        lateralities = db.lateralities
        num_patients = node_matrix.shape[0]
        similarity = np.zeros((num_patients, num_patients))
        for i1 in range(num_patients):
            for i2 in range(num_patients):
                if i1 == i2:
                    continue
                p1 = node_matrix[i1, :]
                p2 = node_matrix[i2, :]
                subsite1 = subsites[i1]
                subsite2 = subsites[i2]
                laterality1 = lateralities[i1]
                laterality2 = lateralities[i2]
                same_laterality = 1 if (laterality1 == laterality2) else 0
                same_subsite = 1 if subsite1 == subsite2 else 0
                similarity[i1, i2] = self.similarity(p1, p2, same_laterality, same_subsite)
        return similarity

    def similarity(self, x, y, j, k):
        numerator = x.dot(y)
        denominator = x.dot(x) + y.dot(y) - x.dot(y)
        if numerator == 0 or denominator == 0:
            return 0
        return numerator/denominator
    
class TreeEstimator():
    subsite_map = {'BOT': 0, 'GPS': 1, 'NOS': 2, 'Soft palate': 3, 'Tonsil': 4}
    laterality_map = {'Bilateral': 0, 'L': 1, 'R': 2}
    def __init__(self, num_pca_components = 10, n_estimators =10, min_samples_split = 4, max_depth = None):
        from sklearn.ensemble import RandomForestRegressor
        self.model = RandomForestRegressor(min_samples_split = min_samples_split, 
                                           n_estimators = n_estimators,
                                           max_depth = max_depth)
        self.num_pca_components = num_pca_components

    def evaluate(self, data):
        x = self.get_input_features(data)
        y = data.doses
        i = 0
        errors = []
        while i < (data.get_num_patients() - 10):
            values = np.arange(i, i+10, 1)
            x_train = np.delete(x, values, axis = 0)
            y_train = np.delete(y, values, axis = 0)
            x_test = x[values, :]
            y_test = y[values, :]
            i += 10
            y_predict = self.predict_doses( x_train, y_train, x_test, y_test)
            error = self.get_error(y_predict, y_test)
            errors.append(error)
        return np.mean(errors, axis = 0)
        
    def predict_doses(self, x_train, y_train, x_test, y_test):
        self.model.fit(x_train, y_train)
        y_predict = self.model.predict(x_test)
        return y_predict
        
    def get_error(self, predicted_doses, dose_matrix):
        differences = np.abs(predicted_doses - dose_matrix)
        percent_error = np.sum(differences, axis = 1)/np.sum(dose_matrix, axis = 1)
        return percent_error
    
    def get_input_features(self, data):
        num_patients = data.get_num_patients()
        pca = lambda x: Rankings.pca(x, self.num_pca_components)
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
        laterality = np.vectorize(TreeEstimator.laterality_map.__getitem__)(laterality)
        subsites = data.subsites.reshape(num_patients, 1)
        subsites = np.vectorize(TreeEstimator.subsite_map.__getitem__)(subsites)
        total_doses = data.prescribed_doses.reshape(num_patients, 1)
        clusters = data.classes.reshape(num_patients, 1)
        features = np.hstack([distances, lymph_nodes, tumor_volumes, total_doses, laterality, subsites])
        return features
            