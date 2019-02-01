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

    def experimental(p1, p2, weights = np.array([1,.3, 1])):
        #this is basically just an ensemble of different distances metrics at this point
        scores = np.zeros((2,))
        #ssim seems to do better than other things?
        scores[0] = Rankings.ssim(p1, p2)
        #this one is the most important
        make_matrix = lambda x: x.reshape(len(x),1)*x.reshape(1,len(x))
        scores[1] = compare_ssim( make_matrix(p1.tumor_distances),
              make_matrix(p2.tumor_distances))
#        scores[2] = 1 if p1.laterality == p2.laterality else 0
#        percent_different = lambda x,y: 1- np.abs(x - y)/(x + y + .0000001)
#        scores[3] = percent_different(p1.tumor_volume, p2.tumor_volume)
#        scores[4] = Rankings.volume_mse(p1, p2)
#        scores[5] = percent_different(p1.prescribed_dose, p2.prescribed_dose)
        final_score = np.sum(scores*weights)/np.mean(weights)
        return(final_score)


class Patient():

    def __init__(self, distances, doses, p_id, position, info):
        #patient ID number
        self.id = p_id
        #basically ordinality of the id, so where it will be in an index
        self.pos = position
        self.laterality = info['Tm Laterality (R/L)']
        self.age = info['Age at Diagnosis (Calculated)']
        self.prescribed_dose = info['Total dose']
        centroid_data = self.get_doses_file_info(doses)
        self.doses = centroid_data[:, 4]
        #####normalize to total dose and then dose proportions
        self.total_dose = np.sum(self.doses)
        ######################
        self.volumes = centroid_data[:, 3]
        self.centroids = centroid_data[:, 0:3]
        self.distances = self.gen_distance_matrix(distances)
        (self.gtvp_dists, self.gtvn_dists) = self.get_tumor_distances(distances)
        #store the entries without gtvp for future study
        (self.tumor_volume, self.tumor_distances) = self.get_main_tumor()
        self.check_missing_organs(distances, doses)
        #report if there is no primary tumor
        if self.tumor_volume == 0 or np.sum(self.tumor_distances) == 0:
            Constants.no_tumor.append(self.id)

    def check_missing_organs(self, distances, doses):
        #check if any organs are missing using the dose file, and store them
        organs = set(Constants.organ_list[:])
        dose_organs = set(doses['ROI'].unique())
        diff = organs - dose_organs
        #for missing_organ in diff:
            #print('patient ', self.id, ' at index ', self.pos, ' is missing organ ', missing_organ)
        if len(diff) > 0:
            Constants.missing_organs[self.pos] = {'id': self.id, 'organs': diff}
        return

    def get_doses_file_info(self, doses):
        #rename the columns so they're consistent
        doses.columns = Constants.centroid_file_names
        #move centroids so the center of the cloud is at zero?
        centroids = self.center_centroids(doses)
        centroids = centroids.set_index('ROI')
        #extract the primary tumor info.
        try:
            gtvp = centroids.loc['GTVp']
            self.gtvp_volume = gtvp.volume
        except:
            self.gtvp_volume = float(0)
        try:
            self.gtvp_position = gtvp[['x','y','z']].values
        except:
            self.gtvp_position = np.array([0,0,0])
        #extract a secondary tumor (only gets the first one?)
        #several patients have no gtvp but a gtvn
        try:
            gtvn = centroids.loc['GTVn']
            if not isinstance(gtvn.volume, float):
                Constants.multiple_gtvn.append(self.id)
                gtvn = gtvn.iloc[0]
            self.gtvn_volume = gtvn.volume
        except:
            self.gtvn_volume = float(0)
        try:
            self.gtvn_position = gtvn[['x','y','z']].values
        except:
            self.gtvn_position = np.array([0,0,0])
        #get the info the centers, volumes, nad doses for all the things
        centroid_matrix = np.zeros((Constants.num_organs,5)) #row = x,y,z,volume,dose
        for idx in range(0, Constants.num_organs):
            organ = Constants.organ_list[idx]
            try:
                organ_entry = centroids.loc[organ]
                centroid_matrix[idx, 0:3] = organ_entry[['x','y','z']].values
                centroid_matrix[idx, 3] = organ_entry.volume
                centroid_matrix[idx, 4] = organ_entry.mean_dose
            except:
                pass
                #print('patient ', self.id, ' is missing organ ', organ, ' centroid data')
        return(centroid_matrix)

    def center_centroids(self, centroids):
        #subtract off the mean so the pointcloud is centered at 0
        #should I just use a reference organ instead?  or rotate?
        centroids.x -= centroids.x.mean()
        centroids.y -= centroids.y.mean()
        centroids.z -= centroids.z.mean()
        return(centroids)

    def gen_distance_matrix(self, dists):
        #generates a symetric 45x45 matrix of organ-organ distances
        dist_matrix = np.zeros(( Constants.num_organs, Constants.num_organs))
        dists = dists.set_index(['Reference ROI', 'Target ROI']).sort_index()
        for row in range(0, Constants.num_organs):
            for col in range(row + 1, Constants.num_organs):
                organ1 = Constants.organ_list[row]
                organ2 = Constants.organ_list[col]
                try:
                    dist_matrix[row, col] = (dists.loc[organ1, organ2])['Eucledian Distance (mm)']
                except:
                    dist_matrix[row, col] = 0
        dist_matrix += np.transpose(dist_matrix)
        return(dist_matrix)

    def get_tumor_distances(self, dists):
        #gets the tumor-organ distances
        gtvp_dists = np.zeros((Constants.num_organs,))
        gtvn_dists = np.zeros((Constants.num_organs,))
        dists = dists.set_index(['Reference ROI', 'Target ROI']).sort_index()
        for idx in range(0, Constants.num_organs):
            organ = Constants.organ_list[idx]
            try:
                tumor_row = dists.loc['GTVp', organ]
                gtvp_dists[idx] = tumor_row['Eucledian Distance (mm)']
            except:
                try:
                    tumor_row = dists.loc[organ, 'GTVp']
                    gtvp_dists[idx] = tumor_row['Eucledian Distance (mm)']
                except:
                    gtvp_dists[idx] = float(0)
            try:
                tumor_row = dists.loc['GTVn', organ]
                gtvn_dists[idx] = tumor_row['Eucledian Distance (mm)']
            except:
                try:
                    tumor_row = dists.loc[Constants.organ_list[idx], 'GTVn']
                    gtvn_dists[idx] = tumor_row['Eucledian Distance (mm)']
                except:
                    gtvn_dists[idx] = float(0)
        return((gtvp_dists, gtvn_dists))

    def get_main_tumor(self):
        #basically gives a proxy so we use only the most important tumor?
        if self.gtvn_volume > 0 and self.gtvp_volume > 0:
            tumor_volume = (self.gtvn_volume + self.gtvp_volume)
            tumor_distances = (self.gtvn_dists*self.gtvn_volume + self.gtvp_dists*self.gtvp_volume)/(self.gtvn_volume + self.gtvp_volume)
        elif self.gtvn_volume > 0 and self.gtvp_volume == 0:
            tumor_volume = self.gtvn_volume
            tumor_distances = self.gtvn_dists
        else:
            tumor_volume = self.gtvp_volume
            tumor_distances = self.gtvp_dists
        if (self.gtvn_volume != 0 and np.sum(self.gtvn_dists) == 0) or (
                self.gtvn_volume == 0 and np.sum(self.gtvn_dists) != 0) or (
                        self.gtvp_volume != 0 and np.sum(self.gtvp_dists) == 0) or (
                            self.gtvp_volume == 0 and np.sum(self.gtvp_dists) != 0):
            print('patient ', self.id, 'is having some issues with tumor consistency')
        return(tumor_volume, tumor_distances)

class PatientSet():

    def __init__(self, outliers = [], root = 'patient_files\\'):
        outliers = outliers
        (self.patients, self.doses, self.total_doses, self.num_patients) = self.read_patient_data(root, outliers)
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
        self.id_map = {max([int(x) for x in findall('[0-9]+', file)]): distance_files.index(file)  for file in distance_files}
        ids = sorted(list(self.id_map.keys()))
        for outlier_id in sorted(outliers, reverse = True):
            if outlier_id in self.id_map:
                pos = self.id_map[outlier_id]
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
            distances = pd.read_csv(distance_files[patient_index], usecols = [0,1,2])
            #renames anything that is equivalent to GTVp/GTVn to the correct format
            distances.replace(Constants.tumor_aliases, inplace = True)
            doses = pd.read_csv(dose_files[patient_index])
            doses.replace(Constants.tumor_aliases, inplace = True)
            info = metadata.loc[ids[patient_index]]
            new_patient = Patient(distances, doses,
                                  ids[patient_index], patient_index, info)
            patients[patient_index] = new_patient
            dose_matrix[patient_index, :] = new_patient.doses
            total_dose_vector[patient_index] = new_patient.total_dose
        return((patients, dose_matrix, total_dose_vector, num_patients))

    def get_by_id(self, p_id):
        if p_id in self.id_map:
            pos = self.id_map[p_id]
            return(self.patients[pos])
        else:
            print('invalid id')

    def get_patients(self):
        return list(self.patients.values())

    def gen_score_matrix(self, rank_function = 'ssim', weights = 1):
        #generates a score matrix based on a rank function
        #function should rank more similar people with a higher number
        scores = np.zeros((self.num_patients, self.num_patients))
        for row in range(0, self.num_patients):
            for col in range(row + 1, self.num_patients):
                scores[row, col] = self.compare_traits(self.patients[row], self.patients[col],
                      rank_function = rank_function, weights = weights)
        #formats it as a symetric matrix with a zero diagonal
        scores += scores.transpose()
        #basically normalize the score so the max is 1?
        scores = scores/scores.max()
        return(scores)

    def predict_doses(self, rank_function = 'ssim', weights = 1, num_matches = 5):
        #generates an ndarray of dose estimates based on algorithm parameters
        estimates = np.zeros(self.doses.shape)
        ranks = self.gen_score_matrix(rank_function = rank_function, weights = weights)
        for patient_idx in range(0, self.num_patients):
            sorted_matches = np.argsort(-ranks[patient_idx,:])
            top_matches = sorted_matches[0:num_matches]
            scores = ranks[patient_idx, tuple(top_matches)]
            matched_dosages = self.doses[tuple(top_matches), :]
            #weight things by their scores
            for match_idx in range(0, num_matches):
                matched_dosages[match_idx, :] = scores[match_idx]*matched_dosages[match_idx, :]
            estimates[patient_idx, :] = np.mean(matched_dosages, axis = 0)/np.mean(scores)
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

    def evaluate(self, rank_function = 'ssim', weights = 1, num_matches = 5):
        #gives a bunch of different metrics for evaluating a given metric
        estimates = self.predict_doses(rank_function, weights, num_matches)
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

    def run_study(self, max_matches = 20, rank_function = 'ssim', weights = 1):
        #tests out a metric for a range of difference potential matches and gives a minimum
        error_hist = []
        for num_matches in range(2, max_matches):
            estimates = self.predict_doses(rank_function, weights, num_matches)
            error = np.mean(np.abs(self.doses - estimates))
            error_hist.append(error)
        print(rank_function, ': error of', min(error_hist), ' at ', np.argmin(error_hist) + 2)
        return(error_hist)

#db = pickle.load(open('data\\patient_data_v2_only.p', 'rb'))

#db = PatientSet(root = 'data\\patients_v2\\', outliers = [239,2009, 10034, 10164])
#pickle.dump(db, open('data\\patient_data_v2_only.p', 'wb'))


weights = np.array([1, .3]) #ssim, tumor_distance ssim, laterality,tumor volume percent similarity, volume vector mse,
result = db.evaluate(rank_function = 'experimental', weights = weights, num_matches = 5)
print(result['mean_error'])
print(len(np.where(result['patient_mean_error'] > 8)[0]))

