# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 10:15:21 2019

@author: Andrew
"""
import numpy as np
from Constants import Constants
from ErrorChecker import ErrorChecker
import Metrics

class KnnEstimator():
    #class that works as a modified knn predictor for doses
    #mainn function is evaluate, which gives you the precent prediction error given a patientset and a similarity matrix
    def __init__(self, match_threshold = .65, match_type = 'threshold', min_matches = 8):
        #match type is how the number of matches are selected
        #give clusters to make it based on the class.  
        #threshold uses similarity score, uses max(min_matches, patients with score > match threshold)
        
        self.min_matches = min_matches
        #similarity needed to be a match
        self.match_threshold = match_threshold
        self.match_type = match_type
        return

    def predict_doses(self, similarity_matrix, data):
        dose_matrix = data.doses
        clusters = data.classes
        predicted_doses = np.zeros(dose_matrix.shape)
        outliers = ErrorChecker().get_data_outliers(data.doses)
        similarity = np.copy(similarity_matrix)
        for p in range( dose_matrix.shape[0] ):
            num_matches = self.get_num_matches(p, similarity, clusters)
            scores = similarity[p, :]
            if p not in outliers:
                scores[list(outliers)] = 0
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
#            print(patient, 'doesnt have any matches')
            predicted_doses = dose_matrix.mean(axis=0) #return average if bad
        return predicted_doses

    def get_matches(self, similarity_matrix, data):
        dose_matrix = data.doses
        clusters = data.classes
        #should return a list of matched patients
        matches = []
        outliers = ErrorChecker().get_data_outliers(data.doses)
        similarity = np.copy(similarity_matrix)
        for p in range( dose_matrix.shape[0] ):
            num_matches = self.get_num_matches(p, similarity, clusters)
            scores = similarity[p, :]
            if p not in outliers:
                scores[list(outliers)] = 0
            args = np.argsort(-scores)
            
            args = args[0 : num_matches] + 1
            matches.append(args)
        return(matches)

    def get_num_matches(self, p, similarity, clusters):
        #for later better use probs
        if self.match_type == 'threshold':
            good_matches= len(np.where(similarity[p,:] > self.match_threshold)[0])
            matches = max([self.min_matches, good_matches])
        elif self.match_type == 'clusters':
            num_cluster_values = len(np.where(clusters == clusters[p])[0])
            matches = max([int(np.sqrt(num_cluster_values) + 1), self.min_matches])
        else:
            matches = self.min_matches
        return matches

    def get_error(self, predicted_doses, dose_matrix):
        differences = np.abs(predicted_doses - dose_matrix)
        percent_error = np.sum(differences, axis = 1)/np.sum(dose_matrix, axis = 1)
        return percent_error

    def evaluate(self, similarity_matrix, data):
        predicted_doses = self.predict_doses(similarity_matrix, data)
        percent_error = self.get_error(predicted_doses, data.doses)
        return(percent_error)
        
class TsimModel():
    #orginal-ish similarity model that gives a similarity matrix from get_simirity based on a spatial ssim 
    def __init__(self, max_distance = 50, patients = None, organs = None, 
                 similarity_function = None, use_classes = False):
        self.max_distance = max_distance
        #for subsetting the data later
        self.patients = patients
        self.organs = organs
        self.use_classes = use_classes
        if similarity_function is None:
            self.similarity_function = self.local_ssim
        else:
            self.similarity_function = similarity_function

    def get_adjacency_lists(self, organ_distance_matrix, organs = None):
        if organs is None:
            organs = self.organs
        #this code is much simpler than expected
        organ_distances = organ_distance_matrix[organs][:, organs]
        adjacency_lists = []
        for row in organ_distances:
            adjacent_args = np.argwhere(row < self.max_distance)
            adjacency_lists.append(adjacent_args.ravel())
        return adjacency_lists

    def get_similarity(self, data, denoise = True):
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
        num_patients, num_organs = distances.shape
        score_matrix = np.zeros((num_patients, num_patients))
        for patient1 in range(0, num_patients - 1):
            for patient2 in range(patient1 + 1, num_patients):
                scores = self.pairwise_similarity(patient1, patient2, 
                                                     distances, volumes,
                                                     clusters, adjacency)
                score_matrix[patient1, patient2] = scores
        score_matrix += np.transpose(score_matrix)
#        scale to between 0 and .99
        score_matrix = (score_matrix - score_matrix.min())
        score_matrix = .99*score_matrix/score_matrix.max()
        for p in range(0, num_patients):
            score_matrix[p,p] = 0
        return score_matrix

    def pairwise_similarity(self, patient1, patient2, distances, volumes, clusters, adjacency = None):
        if self.use_classes and clusters[patient1] != clusters[patient2]:
            return 0
        scores = []
        for organ in range(distances.shape[1]):
            adjacent_args = adjacency[organ]
            if len(adjacent_args) < 1:
                continue
            d1 = distances[patient1, adjacent_args]
            d2 = distances[patient2, adjacent_args]
            v1 = volumes[patient1, adjacent_args]
            v2 = volumes[patient2, adjacent_args]
            scores.append( self.similarity_function(d1,d2,v1,v2) )
        return np.mean(scores)


    def local_ssim(self, x,y,v = None, w = None):
        c1 = .000001
        c2  = .000001
        mean_x = np.mean(x)
        mean_y = np.mean(y)
        covariance = np.cov(x,y)
        numerator = (2*mean_x*mean_y + c1) * (covariance[0,1] + covariance[1,0] + c2)
        denominator = (mean_x**2 + mean_y**2 + c1)*(np.var(x) + np.var(y) + c2)
        if v is not None and w is not None:
            mean_v = np.mean(v)
            mean_w = np.mean(w)
            numerator *= (2*mean_v*mean_w + c1)
            denominator *= (mean_v**2 + mean_w**2 + c1)
        if denominator > 0:
            return numerator/denominator
        else:
            print('error, zero denomiator in ssim function')
            return 0
        
class TJaccardModel(TsimModel):
    #smalle variant of the tsim model that uses jaccard by default.  saves like two lines of work later
    def __init__(self, max_distance = 50, patients = None, organs = None, 
                 similarity_function = None, use_classes = False):
        super(TJaccardModel, self).__init__(max_distance, patients, organs, similarity_function, use_classes)
        self.similarity_function = Metrics.jaccard_distance
        if organs is None:
            self.organs = [Constants.organ_list.index(o) for o in Constants.optimal_jaccard_organs]

class OsimModel(TsimModel):
    #variant that calculates the tsim similarity using organ-organ distances, rather than tumor-organ distances
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
        adjacency = self.get_adjacency_lists(data.organ_distances, organs)
        distances = data.all_organ_distances
        volumes = data.volumes
        clusters = data.classes
        scores = self.similarity(adjacency, distances, volumes, clusters)
        return scores
    
    def similarity(self, adjacency, distances, volumes, clusters, similarity_function = None):
        if similarity_function is None:
            similarity_function = self.local_ssim
        num_patients, num_organs = (distances.shape[2], distances.shape[0])
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
                    d1 = distances[adjacent_args, :, patient1][:, adjacent_args]
                    d2 = distances[adjacent_args, :, patient2][:, adjacent_args]
                    similarity_score = similarity_function(d1,d2)
                    scores.append( similarity_score )
                score_matrix[patient1, patient2] = np.mean(scores)
        score_matrix += np.transpose(score_matrix)
        #scale to between 0 and .99
        score_matrix = .99*(score_matrix - score_matrix.min())/(score_matrix.max() - score_matrix.min())
        return score_matrix

class SimilarityFuser():
    #class that uses (logistic regression) to map a list of similarity scores to a single score
    #attempts to classify each vector as a neighbor (dose error < some number) 
    #and the match probability is used as the similarity
    def __init__(self, model = None, min_matches = 4, max_error = .1):
        self.min_matches = min_matches
        self.max_error = max_error
        if model is None:
            from sklearn.linear_model import LogisticRegression
            model = LogisticRegression(class_weight = 'balanced',
                                       solver = 'lbfgs',
                                       max_iter=500,
                                       random_state = 0)
        self.model = model
        
    def get_similarity(self, db, similarity_matrices, classes = None):
        [x,y, positions] = self.extract_features(similarity_matrices, db.get_num_patients(), db, classes)
        final_similarity = np.zeros(similarity_matrices[0].shape)
        x = (x-x.min(axis=0))/(x.max(axis=0) - x.min(axis=0))
        for p in range(db.get_num_patients() - 1):
            test_pairs = np.where(positions == p)[0] 
            x_train = np.delete(x, test_pairs, axis = 0)
            y_train = np.delete(y, test_pairs, axis = 0)
            predict_pairs  = np.where(positions[:,0] == p)[0]
            x_test = x[predict_pairs,:]
            self.model.fit(x_train, y_train)
            y_pred = self.model.predict_proba(x_test)[:, 1]
            final_similarity[p, :] = np.insert(y_pred, p, 0)
#        final_similarity += final_similarity.transpose()
        return(final_similarity)
        
        
    def extract_features(self, similarities, num_patients, data, classes):
        true_similarity = self.get_true_matches(data)
        if classes is not None:
            classes = (classes - np.min(classes))/(np.max(classes) - np.min(classes))
        x = []
        y = []
        positions = []
        for p1 in range(num_patients):
            for p2 in range(num_patients):
                if p1 == p2:
                    continue
                x_row = []
                for similarity in similarities:
                    x_row.append(similarity[p1, p2])
                if classes is not None:
                    target_class = classes[p2]
                    x_row.append(target_class)
                x.append(x_row)
                y.append(true_similarity[p1, p2])
                positions.append([p1,p2])
        x = np.array(x)
        y = np.array(y)
        positions = np.array(positions)
        return [x, y, positions]
    
    def get_true_matches(self, data, negative_class = 0):
        min_matches = self.min_matches
        dose_error = self.get_match_error(data)
        match_matrix = np.zeros(dose_error.shape) + negative_class
        n_patients = data.get_num_patients()
        for p in range(n_patients):
            max_error = self.max_error
            errors = dose_error[p, :]
            matches = []
            max_error = max_error
            while len(matches) < min_matches:
                matches = np.where(errors < max_error)[0]
                max_error = max_error + .01
            match_matrix[p, matches] = 1
        return match_matrix.astype('int32')
        
    def get_match_error(self, data):
        n_patients = data.get_num_patients()
        doses = data.doses
        error_matrix = np.zeros((n_patients, n_patients))
        for p1 in range(n_patients):
            for p2 in range(p1 + 1, n_patients):
                dose_difference = np.abs(doses[p1,:] - doses[p2, :])
                error_matrix[p1, p2] = np.sum(dose_difference)/np.sum(doses[p1, :])
        error_matrix += error_matrix.transpose()
        return error_matrix