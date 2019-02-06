# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 17:22:31 2019
@author: Andrew
"""
from glob import glob
from re import findall
import numpy as np
import pandas as pd
from collections import OrderedDict
from skimage.measure import compare_ssim, compare_mse
import random
import matplotlib.pyplot as plt
import pickle
from Constants import Constants
from Patient import Patient

class Rankings():
    #ranking functions that generate a score, takes in pateint objects

    def vector_ssim(p1, p2):
        upper_triangle = np.triu_indices(len((p1.distances)))
        d1 = p1.distances[upper_triangle].ravel()
        d2 = p2.distances[upper_triangle].ravel()
        return(compare_ssim(d1, d2, win_size = 5))

    def volume_weights_ssim(p1, p2):
        volume_matrix1 = p1.volumes.reshape(45,1)*p1.volumes.reshape(1,45)
        volume_matrix2 = p2.volumes.reshape(45,1)*p2.volumes.reshape(1,45)
        d1 = p1.distances*volume_matrix1
        d2 = p2.distances*volume_matrix2
        return(compare_ssim(d1, d2, win_size = 3))

    def ssim(p1, p2):
        return(compare_ssim(p1.distances, p2.distances, win_size = 3))

    def ssim_with_laterality(p1, p2, weights = np.array([1,2,.05])):
        scores = np.zeros((3,))
        scores[0] = Rankings.ssim(p1,p2)
        scores[1] = 1 if p1.laterality == p2.laterality else 0
        scores[2] = 1 - np.abs(float(p1.gtvp_volume) - float(p2.gtvp_volume))/(
                float(p1.gtvp_volume) + float(p2.gtvp_volume) + .000001)
        final_score = np.sum(scores*weights)/np.mean(weights)
        return(final_score)

    def geometric_distance(p1,p2):
        dist = np.sqrt(np.sum((p1 - p2)**2))
        return(1/(dist + .000001))


    def mse(p1,p2):
        return( 1/( compare_mse(p1.distances, p2.distances) + .000001) )

    def volume_mse(p1, p2):
        return( 1/(compare_mse(p1.volumes, p2.volumes) + .000001) )

    def emd(patient_1, patient_2):
        #simplified earth movers distance - here it's just work done to move organs in the less massive
        #patient into the same positions of the more massive patient
        volumes = np.abs(patient_1.volumes - patient_2.volumes)
        dists = np.sqrt(np.sum((patient_1.centroids - patient_2.centroids)**2, axis = 1))
        work = np.sum(dists*volumes)
        #we this to work with the other functions so >0 and high scores = closer
        return(1/(work + .000001))

    def min_dose_error(p1, p2):
        error = np.mean(np.abs(p1.doses - p2.doses))
        return(1/(error + .000001))

    def experimental(p1, p2, weights):
        #this is basically just an ensemble of different distances metrics at this point
        scores = np.zeros((len(weights),))
        make_matrix = lambda x: x.reshape(len(x),1)*np.absolute(x.reshape(1,len(x)))
        percent_different = lambda x,y: 1- np.abs(x - y)/(x + y + .0000001)
        if(weights[0] > 0):
            #ssim seems to do better than other things?
            scores[0] = Rankings.ssim(p1, p2)
        if(weights[1]) > 0:
            #this one is the most important
            scores[1] = compare_ssim( make_matrix(p1.tumor_distances),
                  make_matrix(p2.tumor_distances), win_size = 3)
        if(weights[2] > 0):
            #difference in prescribed dose 'total dose' for the tumor
            scores[2] = percent_different(p1.prescribed_dose, p2.prescribed_dose)
        if(weights[3] > 0):
            #differences in tumor volume?
            scores[3] = percent_different(p1.tumor_volume, p2.tumor_volume)
        if(weights[4] > 0):
            scores[4] = Rankings.volume_mse(p1,p2)
        if(weights[5] > 0):
            scores[5] = Rankings.emd(p1,p2)
        #I was normalizing here, but moved it so I can rescale the data
        #final_score = np.sum(scores*weights)/np.mean(weights)
        return(scores)

class PatientSet():

    def __init__(self, patient_set = None, outliers = [], root = 'patient_files\\'):
        if patient_set is None:
            outliers = outliers
            self.total_dose_predictor = None
            (self.patients, self.doses, self.total_doses, self.num_patients, self.id_map) = self.read_patient_data(root, outliers)
        else:
            self.total_dose_predictor = patient_set.total_dose_predictor
            self.patients = patient_set.patients
            self.doses = patient_set.doses
            self.total_doses = patient_set.total_doses
            self.num_patients = patient_set.num_patients
            self.id_map = patient_set.id_map
        print('\npatient data loaded...\n')

    def read_patient_data(self, root, outliers):
        #sorts by size of largest integer string, which is the id for our files
        file_sort = lambda x: sorted(x, key =
                                     lambda file:
                                         max([int(x) for x in findall("[0-9]+", file)])
                                )
        distance_files = file_sort(glob(root + '**/*distances.csv'))
        dose_files = file_sort(glob(root + '**/*centroid*.csv'))
        #maps ids to position so I can do that too?
        id_map = {max([int(x) for x in findall('[0-9]+', file)]): distance_files.index(file)  for file in distance_files}
        ids = sorted(list(id_map.keys()))
        #delete patient files with an id in the outliers
        for outlier_id in sorted(outliers, reverse = True):
            if outlier_id in id_map:
                pos = id_map[outlier_id]
                del distance_files[pos]
                del dose_files[pos]
                del ids[pos]
        metadata_file = 'data\\patient_info.csv'
        assert(len(distance_files) == len(dose_files))
        #maps a position 0-len(files) to the dummy id for a patient
        num_patients = len(ids)
        metadata = pd.read_csv(metadata_file,
                               index_col = 0, #index is the "Dummy ID"
                               usecols = [0,1,2,3,4,5,6,7,8,9,10,11,31]
                               ).loc[ids]
        patients = OrderedDict()
        dose_matrix = np.zeros((num_patients, len(Constants.organ_list)))
        total_dose_vector = np.zeros((num_patients,))
        #putting all the data into a patient object for further objectification
        for patient_index in range(0, num_patients):
            #these are indexed by name of organ
            #we only use 3 rows but half of them have a comma missing in the header between the last two rows
            distances = pd.read_csv(distance_files[patient_index],
                                    usecols = [0,1,2]).dropna()
            #renames anything that is equivalent to GTVp/GTVn to the correct format
            distances.replace(Constants.tumor_aliases, inplace = True)
            doses = pd.read_csv(dose_files[patient_index],
                                usecols = [0,1,2,3,4,5,6,7]).dropna()
            doses.columns = Constants.centroid_file_names
            doses.replace(Constants.tumor_aliases, inplace = True)
    #            doses = doses.astype({'volume': np.int32, 'mean_dose': np.float32,
    #                                         'x': np.float32, 'y': np.float32, 'z': np.float32})
            info = metadata.loc[ids[patient_index]]
            new_patient = Patient(distances, doses,
                                  ids[patient_index], patient_index, info)
            patients[patient_index] = new_patient
            dose_matrix[patient_index, :] = new_patient.doses
            total_dose_vector[patient_index] = new_patient.total_dose
        return((patients, dose_matrix, total_dose_vector, num_patients, id_map))

    def get_by_id(self, p_id):
        if p_id in self.id_map:
            pos = self.id_map[p_id]
            return(self.patients[pos])
        else:
            print('invalid id')

    def get_patients(self):
        return list(self.patients.values())

    def gen_score_matrix(self, weights, rank_function = 'ssim'):
        #generates a score matrix based on a rank function
        #function should rank more similar people with a higher number
        scores = np.zeros((self.num_patients, self.num_patients,  len(weights)))
        for row in range(0, self.num_patients):
            for col in range(row + 1, self.num_patients):
                scores[row, col, :] = self.compare_traits(self.patients[row], self.patients[col],
                      rank_function = rank_function, weights = weights)
        rescale = lambda x: (x - np.min(x))/(np.max(x) - np.min(x) + .000001)
        for score_idx in range(0, len(weights)):
            scores[:, :, score_idx] = rescale(scores[:,:, score_idx])
        scores *= weights
        scores = np.sum(scores, axis = 2)/np.mean(weights)
        #formats it as a symetric matrix with a zero diagonal
        scores += scores.transpose()
        #basically normalize the score so the max is 1?
        scores = scores/scores.max()
        return(scores)

    def estimate_patient_doses(self, ranks, num_matches):
        sorted_matches = np.argsort(-ranks)
        top_matches = sorted_matches[0:num_matches]
        scores = ranks[top_matches]
        matched_dosages = self.doses[tuple(top_matches), :]
        #weight things by their scores
        for match_idx in range(0, num_matches):
            matched_dosages[match_idx, :] = scores[match_idx]*matched_dosages[match_idx, :]
        patient_estimates = np.mean(matched_dosages, axis = 0)/np.mean(scores)
        return(patient_estimates)

    def predict_doses(self, rank_function = 'ssim', weights = 1, num_matches = 5,
                      td_rank_function = None, td_weights = None):
        #generates an ndarray of dose estimates based on algorithm parameters
        #rank_function and weights are the functions used to match dose
        estimates = np.zeros(self.doses.shape)
        ranks = self.gen_score_matrix(rank_function = rank_function, weights = weights)
        for patient_idx in range(0, self.num_patients):
            rank_row = ranks[patient_idx, :]
            estimates[patient_idx, :] = self.estimate_patient_doses(rank_row, num_matches)
        #if a seprate prediction is set for total dose, use that
        if self.total_dose_predictor is not None:
            #normalize the estimates by patient
            estimates /= np.sum(estimates, axis = 1).reshape((self.num_patients,1))
            x = self.gen_patient_feature_matrix()
            total_dose_prediction = self.total_dose_predictor.predict(x)
            estimates *= total_dose_prediction.reshape((self.num_patients,1))
        elif td_rank_function is not None and td_weights is not None:
            estimates /= np.sum(estimates, axis = 1).reshape((self.num_patients,1))
            total_dose_estimates = np.zeros((self.num_patients,))
            td_ranks = self.gen_score_matrix(rank_function = td_rank_function, weights = td_weights)
            for p_idx in range(0,self.num_patients):
                td_rank_row = td_ranks[p_idx, :]
                estimates_for_total_dose = self.estimate_patient_doses(td_rank_row, num_matches)
                total_dose_estimates[p_idx] = np.sum(estimates_for_total_dose)
            estimates *= total_dose_estimates.reshape((self.num_patients,1))
        return(estimates)

    def get_average_patient_data(self, key = 'all'):
        #generates a dictionary with *some* (positions, distances, and volumes) of the
        #average data accross all patients
        avg_centroids = np.zeros((45,3))
        avg_volumes = np.zeros((45,))
        avg_distances = np.zeros((45,45))
        avg_tumor_distances = np.zeros((45,))
        avg_tumor_volume = 0.0
        for patient in self.get_patients():
            avg_centroids += patient.centroids
            avg_volumes += patient.volumes
            avg_distances += patient.distances
            avg_tumor_distances += patient.tumor_distances
            avg_tumor_volume += patient.tumor_volume
        p_avg = {}
        p_avg['centroids'] = avg_centroids / self.num_patients
        p_avg['volumes'] = avg_volumes / self.num_patients
        p_avg['distances'] = avg_distances / self.num_patients
        p_avg['tumor_distances'] = avg_tumor_distances / self.num_patients
        p_avg['tumor_volume'] = avg_tumor_volume / self.num_patients
        #defaults to a dict, adding in a parameter to only look at one thing
        if key == 'all':
            return(p_avg)
        else:
            return(p_avg[key])

    def set_total_dose_prediction(self, predictor):
        #should be like a skikit model, so takes numppy array with .predict(x)
        #and returns a 1d numpy array of y
        self.total_dose_predictor = predictor

    def gen_patient_feature_matrix(self):
        #function to get a matrix I can try some dose prediction on?
        features = np.zeros((self.num_patients, 14))
        laterality_map = {'L': -1, 'R': 1, 'Bilateral': 0}
#        norm = lambda x: np.sqrt(np.sum(x*x))
        brainstem_pos = Constants.organ_list.index('Brainstem')
        lt_thyroid_pos = Constants.organ_list.index('Lt_thyroid_lobe')
        rt_thyroid_pos = Constants.organ_list.index('Rt_thyroid_lobe')
        for patient_idx in range(0, self.num_patients):
            patient = self.patients[patient_idx]
            features[patient_idx, 0] = patient.gtvp_volume
            features[patient_idx, 1] = patient.gtvn_volume
            features[patient_idx, 2] = patient.prescribed_dose
            features[patient_idx, 3] = laterality_map[patient.laterality]
            features[patient_idx, 4:7] = patient.gtvp_position[:]
            features[patient_idx, 7:10] = patient.gtvn_position[:]
            features[patient_idx, 10] = np.sum(patient.volumes)
            features[patient_idx, 11] = np.sum(patient.tumor_distances[brainstem_pos])
            features[patient_idx, 12] = np.sum(patient.tumor_distances[lt_thyroid_pos])
            features[patient_idx, 13] = np.sum(patient.tumor_distances[rt_thyroid_pos])
        #standarization, not needed for binary trees though
        features = (features - np.mean(features, axis = 0))/np.std(features, axis = 0)
        return(features)

    def evaluate(self, rank_function = 'ssim', weights = 1, num_matches = 5,
                 td_rank_function = None, td_weights = None):
        #gives a bunch of different metrics for evaluating a given metric
        estimates = self.predict_doses(rank_function, weights, num_matches,
                                       td_rank_function, td_weights)
        differences = self.doses - estimates
        patient_mean_error = np.mean(np.abs(differences), axis = 1)
        total_mean_error = np.mean(patient_mean_error)
        total_rmse = np.sqrt(np.mean(differences**2))
        result_dict = {'prediction': estimates,
                       'patient_mean_error': patient_mean_error,
                       'mean_error': total_mean_error,
                       'rmse': total_rmse,
                       'differences': differences}
        return(result_dict)

    def compare_traits(self, p1, p2, weights, rank_function = 'ssim'):
        #calculates an overall scores
        #currently: comparison function, if laterality is equal, difference in tumor volume, tumor distances
        if rank_function == 'ssim':
            score = Rankings.ssim_with_laterality(p1, p2, weights = weights)
        elif rank_function == 'experimental':
            score = Rankings.experimental(p1,p2,weights = weights)
        elif rank_function == 'min_dose_error':
            score = Rankings.min_dose_error(p1,p2)
        elif rank_function == 'random':
            score = random.random()
        elif rank_function == 'mse':
            score = Rankings.mse(p1,p2)
        else:
            print('error, invalid rank method: ', rank_function)
        return(score)

    def run_study(self, max_matches = 20, rank_function = 'ssim', weights = 1,
                  td_rank_function = None, td_weights = None):
        #tests out a metric for a range of difference potential matches and gives a minimum
        error_hist = []
        for num_matches in range(2, max_matches):
            estimates = self.predict_doses(rank_function, weights, num_matches,
                                           td_rank_function, td_weights)
            error = np.mean(np.abs(self.doses - estimates))
            error_hist.append(error)
        print(rank_function, ': error of', min(error_hist), ' at ', np.argmin(error_hist) + 2)
        return(error_hist)

def cluster_organs(db):
    avg = db.get_average_patient_data()
    centroids = avg['centroids']
    named_centroids = zip(Constants.organ_list, centroids)
    from sklearn.cluster import AffinityPropagation
    estimator = AffinityPropagation()
    estimator.fit(centroids)
    organs = [(estimator.labels_[x], centroids[x,:], Constants.organ_list[x]) for x in range(0, Constants.num_organs)]
    organs = sorted(organs, key = lambda x: x[0])
    return([x[2] for x in organs])

def train_total_dose_tree(db):
    x = db.gen_patient_feature_matrix()
    y = db.total_doses[:].reshape((db.num_patients, 1))
    partition = db.num_patients//3
    x_test = x[:partition, :]
    x_train = x[partition:, :]
    y_test = y[:partition]
    y_train = y[partition:]
    from sklearn.tree import DecisionTreeRegressor

    tree = DecisionTreeRegressor(max_depth = 10, criterion = 'mae')
    tree.fit(x_train,y_train)
    return(tree)

#db = pickle.load( open('data\\patient_data_v2_only.p', 'rb' ))
#db = PatientSet(patient_set = None, root = 'data\\patients_v2\\',
#                outliers = Constants.v2_bad_entries + Constants.v2_half_dosed)
#pickle.dump(db, open('data\\patient_data_v2_only.p', 'wb'))


db.set_total_dose_prediction(None)
weights = np.array([1, 1, 0, 0, 1, 0])
td_weights = np.array([1, 1, 0, 0, 1, 0])
num_matches = 5
prediction = db.predict_doses(rank_function = 'experimental',
                              weights = weights,
                              num_matches = num_matches,
                              td_weights = td_weights,
                              td_rank_function = 'experimental')
total_predicted_doses = np.sum(prediction, axis = 1)
base_rmse = np.sqrt(np.sum((db.total_doses - total_predicted_doses)**2)/db.num_patients)
print('not decision tree', base_rmse)
result = db.evaluate(rank_function = 'experimental',
                              weights = weights,
                              num_matches = num_matches,
                              td_weights = td_weights,
                              td_rank_function = 'experimental')
print(' mean error: ',result['mean_error'],' rmse: ', result['rmse'])
print(len(np.where(result['patient_mean_error'] > 8)[0]))


tree = train_total_dose_tree(db)
db.set_total_dose_prediction(tree)
prediction = db.predict_doses(rank_function = 'experimental', weights = weights, num_matches = num_matches)
total_predicted_doses = np.sum(prediction, axis = 1)
base_rmse = np.sqrt(np.sum((db.total_doses - total_predicted_doses)**2)/db.num_patients)
print('decision tree', base_rmse)


result = db.evaluate(rank_function = 'experimental', weights = weights, num_matches = num_matches)
print(result['rmse'], ' ',result['mean_error'])
print(len(np.where(result['patient_mean_error'] > 8)[0]))

