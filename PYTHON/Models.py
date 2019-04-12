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
        
    def predict_doses(self, similarity_matrix, dose_matrix, clusters):
        predicted_doses = np.zeros(dose_matrix.shape)
        for p in range( dose_matrix.shape[0] ):
            num_matches = self.get_num_matches(p)
            scores = similarity_matrix[p, :]
            args = np.argsort(-scores)
            args = args[0 : num_matches]
            matched_scores = scores[args].reshape(len(args), 1)
            matched_doses = dose_matrix[args, :]
            predicted_doses[p,:] = np.mean(matched_scores*matched_doses, axis = 0)/matched_scores.mean()
        return(predicted_doses)

    def get_num_matches(self, row):
        #for later better use probs
        return 10
    
    def evaluate(self, similarity_matrix, dose_matrix, clusters = None):
        predicted_doses = self.predict_doses(similarity_matrix, dose_matrix, clusters)
        differences = np.abs(predicted_doses - dose_matrix)
#        total_mean_error = np.mean(differences)
        percent_error = np.sum(differences, axis = 1)/np.sum(dose_matrix, axis = 1)
#        total_rmse = np.sqrt(np.mean(differences**2))    
        return(percent_error)
        
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
        scores = self.similarity(adjacency, distances, volumes)
        return scores
    
    def similarity(self, adjacency, distances, volumes, similarity_function = None):
        if similarity_function is None:
            similarity_function = self.local_ssim
        num_patients, num_organs = distances.shape
        score_matrix = np.zeros((num_patients, num_patients))
        for patient1 in range(0, num_patients - 1):
            for patient2 in range(patient1 + 1, num_patients):
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
        print(score_matrix)
        score_matrix = .99*(score_matrix - score_matrix.min())/(score_matrix.max() - score_matrix.min())
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