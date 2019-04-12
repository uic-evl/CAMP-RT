# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 17:22:31 2019
@author: Andrew
"""
from glob import glob
from re import findall
import json
import numpy as np
import pandas as pd
from collections import OrderedDict
from Constants import Constants
from Patient import Patient
from ErrorChecker import ErrorChecker 

class PatientSet():

    def __init__(self, outliers = [], root = 'data\\patients_v2*\\', class_name = None, 
                 use_distances = False, use_clean_subset = True):
        if class_name is not None: #default signifies you want to overwrite it
            classes = pd.read_csv('data//rt_plan_clusters.csv',
                                   index_col = 1)
            classes = classes.drop(labels = ['Unnamed: 0'], axis = 1)
            classes.columns = classes.columns.str.strip()
            self.classes = classes[class_name]
            self.num_classes = len(self.classes.unique())
        else:
            self.classes = None
            self.num_classes = 0
        self.read_patient_data(root, outliers, use_distances)
        if use_clean_subset:
            self.clean_values()
        print('\npatient data loaded...\n')

    def read_patient_data(self, root, outliers, use_distances):

        #sorts by size of largest integer string, which is the id for our files
        file_sort = lambda x: sorted(x, key =
                                     lambda file:
                                         max([int(x) for x in findall("[0-9]+", file)])
                                )
        distance_files = file_sort(glob(root + '**/*distances.csv'))
        dose_files = file_sort(glob(root + '**/*centroid*.csv'))
        #maps ids to position so I can do that too?
        ids = self.delete_outliers(outliers, distance_files, dose_files)
        metadata_file = 'data\\patient_info.csv'
        assert(len(distance_files) == len(dose_files))
        #maps a position 0-len(files) to the dummy id for a patient
        num_patients = len(ids)
        metadata = pd.read_csv(metadata_file,
                               index_col = 0, #index is the "Dummy ID"
                               usecols = [0,1,2,3,4,5,6,7,8,9,10,11,18,31]
                               ).loc[ids]

        patients = OrderedDict()
        dose_matrix = np.zeros((num_patients, Constants.num_organs))
        max_dose_matrix = np.zeros((num_patients, Constants.num_organs))
        min_dose_matrix = np.zeros((num_patients, Constants.num_organs))
        centroid_matrix = np.zeros((num_patients, Constants.num_organs, 3))
        total_dose_vector = np.zeros((num_patients,))
        prescribed_dose_vector = np.zeros((num_patients,))
        tumor_distance_matrix = np.zeros((num_patients, Constants.num_organs))
        organ_distance_matrix = np.zeros((Constants.num_organs, Constants.num_organs, num_patients))
        volume_matrix = np.zeros((num_patients,Constants.num_organs))
        laterality_list = []
        subsite_list = []
        gtv_list = []
        classes = np.zeros((num_patients,))
        #putting all the data into a patient object for further objectification

        for patient_index in range(0, num_patients):
            dataset_version = int(findall('patients_v([0-9])', distance_files[patient_index])[0])
            assert(dataset_version in [2,3])
            #these are indexed by name of organ
            #we only use 3 rows but half of them have a comma missing in the header between the last two rows
            distances = pd.read_csv(distance_files[patient_index],
                                    usecols = [0,1,2]).dropna()
            #renames anything that is equivalent to GTVp/GTVn to the correct format
            distances = self.fix_tumor_names(distances)
            doses = pd.read_csv(dose_files[patient_index],
                                usecols = [0,1,2,3,4,5,6,7]).dropna()
            #pateints_v3 dataset has a different way of ording the columns (and different spelling)
            if dataset_version == 2:
                doses.columns = Constants.centroid_file_names_v2
            elif dataset_version == 3:
                doses.columns = Constants.centroid_file_names_v3
            doses = self.fix_tumor_names(doses)
            #misc patient info - laterality, subsite, total dose, etc
            info = metadata.loc[ids[patient_index]]
            group = self.get_patient_class(ids[patient_index], doses.set_index('ROI').mean_dose)
            new_patient = Patient(distances, doses,
                                  ids[patient_index], group,
                                  info, use_distances = use_distances)
            patients[patient_index] = new_patient
            classes[patient_index] = group
            laterality_list.append(new_patient.laterality)
            subsite_list.append(new_patient.tumor_subsite)
            gtv_list.append(new_patient.gtvs)
            dose_matrix[patient_index, :] = new_patient.doses
            max_dose_matrix[patient_index, :] = new_patient.max_doses
            min_dose_matrix[patient_index, :] = new_patient.min_doses
            tumor_distance_matrix[patient_index, :] = new_patient.tumor_distances
            total_dose_vector[patient_index] = new_patient.total_dose
            prescribed_dose_vector[patient_index] = new_patient.prescribed_dose
            volume_matrix[patient_index, :] = new_patient.volumes
            centroid_matrix[patient_index, :, :] = new_patient.centroids
            if use_distances:
                organ_distance_matrix[:, :, patient_index] = np.nan_to_num(new_patient.distances)
        self.doses = np.nan_to_num(dose_matrix)
        self.max_doses = np.nan_to_num(max_dose_matrix)
        self.min_doses = np.nan_to_num(min_dose_matrix)
        self.tumor_distances = np.nan_to_num(tumor_distance_matrix)
        self.volumes = np.nan_to_num(volume_matrix)
        self.classes = np.nan_to_num(classes)
        if use_distances:
            self.organ_distances = np.nan_to_num(organ_distance_matrix).mean(axis = 2)
        else:
            self.organ_distances = self.load_saved_distances()
        self.prescribed_doses = np.nan_to_num(prescribed_dose_vector)
        self.centroids = np.nan_to_num(centroid_matrix)
        self.lateralities = np.array(laterality_list)
        self.subsites = np.array(subsite_list)
        self.ids = np.array(ids)
        self.gtvs = gtv_list
        
    def clean_values(self):
        error_checker = ErrorChecker()
        p = error_checker.get_clean_subset(self)
        p = sorted(p)
        self.doses = self.doses[p]
        self.max_doses = self.max_doses[p]
        self.min_doses = self.min_doses[p]
        self.tumor_distances = self.tumor_distances[p]
        self.volumes = self.volumes[p]
        self.classes = self.classes[p]
        self.prescribed_doses = self.prescribed_doses[p]
        self.centroids = self.centroids[p]
        self.lateralities = self.lateralities[p]
        self.subsites = self.subsites[p]
        self.ids = self.ids[p]
        new_gtvs = []
        for patient in p:
            new_gtvs.append(self.gtvs[patient])
        self.gtvs = new_gtvs
    
    def get_num_patients(self):
        return( self.doses.shape[0] )
    
    def load_saved_distances(self, file = 'data/mean_organ_distances.csv'):
        try:
            distances = pd.read_csv(file, index_col = 0)
            distances = distances.values
        except:
            print('error, no mean-organ distance file found')
            distances = np.zeros((Constants.num_organs, Constants.num_organs))
        return distances
    
    def get_patient_class(self, patient_id, doses):
        #if a vector of classes is used
        group = self.get_default_class(patient_id, doses)
        if self.classes is not None:
            try:
                subclass = self.classes[patient_id]
                group = (group-1)*self.num_classes + subclass
            except:
                print('patient ', patient_id, 'not in class list, defaulting to 0')
        return int(group)

    def get_default_class(self, patient_id, dose_vector):
        full_dose, left_biased = self.check_if_full_dose(dose_vector)
        if not full_dose:
            group = (3 if left_biased == True else 4)
        elif patient_id in Constants.v2_high_throat_dose:
            group = 2
        else:
            group = 1
        return group

    def check_if_full_dose(self, dose_vector):
        #checks difference in sternoceldomastoids to seperate out unilaterally dosed patients?
        #may be used for getting classes eventually?
        try:
            if isinstance(dose_vector, pd.core.series.Series):
                ls = dose_vector.loc['Lt_Sternocleidomastoid_M']
                rs = dose_vector.loc['Rt_Sternocleidomastoid_M']
            else:
                ls_pos = Constants.organ_list.index('Lt_Sternocleidomastoid_M')
                rs_pos = Constants.organ_list.index('Rt_Sternocleidomastoid_M')
                ls = dose_vector[ls_pos]
                rs = dose_vector[rs_pos]
        except:
            print('error in getting dose?')
            ls = 1
            rs = 1
        if np.abs(ls - rs)/max([ls, rs]) < .6:
            full_dose = True
        else:
            full_dose = False
        return(full_dose, (ls > rs))

    def delete_outliers(self, outliers, distance_files, dose_files):
        id_map = {max([int(x) for x in findall('[0-9]+', file)]): distance_files.index(file)  for file in distance_files}
        ids = sorted(list(id_map.keys()))
        #delete patient files with an id in the outliers
        for outlier_id in sorted(outliers, reverse = True):
            if outlier_id in ids:
                pos = id_map[outlier_id]
                del distance_files[pos]
                del dose_files[pos]
                del ids[pos]
        return(ids)

    def fix_tumor_names(self, dataframe):
        #this should probably not need to return anything, but does.
        #replace the aliases for GTVp or GTVn(1) with a consistent name
        dataframe.replace(Constants.tumor_aliases, inplace = True)
        dataframe.replace(to_replace = r'GTV.*N', value = 'GTVn', regex = True, inplace = True)
        return dataframe

    def change_classes(self, class_name):
        if class_name is not None:
            classes = pd.read_csv('data//rt_plan_clusters.csv',
                                       index_col = 1)
            classes = classes.drop(labels = ['Unnamed: 0'], axis = 1)
            classes.columns = classes.columns.str.strip()
            self.classes = classes[class_name]
            self.num_classes = len(self.classes.unique())
        else:
            self.classes = None
            self.num_classes = 0
        for p in self.get_patients():
            p.group = self.get_patient_class(p.id, p.doses)

    def save_organ_distances(self, file = 'data/mean_organ_distances.csv'):
        mean_dists = self.organ_distances.mean(axis = 2)
        if np.sum(mean_dists) == 0:
            print('error, trying to save emtpy organ list')
            return
        else:
            organ_dist_df = pd.DataFrame(mean_dists, index = Constants.organ_list, columns = Constants.organ_list)
        organ_dist_df.to_csv(file)
        
#    def export(self, weights = np.array([0,1]) ,
#               rank_function = 'tumor_organ_ssim',
#               num_matches = 9,
#               patient_data_file = 'data\\patient_dataset.json',
#               score_file = 'data\\all_ssim_scores.csv'):
#        #exports the dataset into the json format peter is using for the frontend
#        data = []
#        scores = self.gen_score_matrix(weights)
#        dose_estimates = self.predict_doses(weights, num_matches)
#        patient_mean_error = np.mean(np.absolute(self.doses - dose_estimates), axis = 1)
#        dose_pca = Rankings.pca(self.doses)
#        distance_pca = Rankings.pca( self.gen_tumor_distance_matrix() )
#        print(distance_pca)
#        for p_idx in range(0, self.num_patients):
#            patient = self.patients[p_idx]
#            entry = patient.to_ordered_dict(dose_estimates[p_idx, :])
#            ssim_scores = scores[p_idx,:]
#            ssim_scores[p_idx] = 1
#            zipped_scores = sorted(zip(ssim_scores, np.arange(1, len(ssim_scores) + 1)),
#                                   key = lambda x: -x[0])
#            patient_scores, internal_ids = zip(*zipped_scores)
#            entry['ID_internal'] = p_idx + 1
#            num_matches = self.get_num_matches(patient) + 1
#            while patient_scores[num_matches - 1] <= 0:
#                num_matches -= 1
#            matches = internal_ids[:num_matches]
#            match_scores = patient_scores[:num_matches]
#            entry['similarity_ssim'] = matches
#            entry['scores_ssim'] = match_scores
#            entry['dose_pca'] = dose_pca[p_idx,:].tolist()
#            entry['distance_pca'] = distance_pca[p_idx, :].tolist()
#            entry['mean_error'] = round(patient_mean_error[p_idx], 4)
#            data.append(entry)
#        #save the vast dictionary of data for the front-end
#        try:
#            def default(o):
#                if isinstance(o, np.int32):
#                    return int(o)
#            with open(patient_data_file, 'w+') as f:  # generate JSON
#                json.dump( data, f, indent=4, default = default)
#            print('successfully save patient data to ', patient_data_file)
#            #save a labeled matrix of similarity scores for other people
#        except:
#            print('error exporting patient data to json')
#        try:
#            raw_scores = self.gen_score_matrix(1, classes = False)
#            score_df = pd.DataFrame(raw_scores, index = self.ids, columns = self.ids)
#            score_df.to_csv(score_file)
#            print('successfully saved similarity score matrix to ', score_file)
#        except:
#            print('error saving ssim score matrix')



#    def get_patients(self):
#        return list(self.patients.values())

#    def get_score_dataframe(self, weights, rank_function):
#        scores = self.gen_score_matrix(weights)
#        score_df = pd.DataFrame(scores, index = self.ids, columns = self.ids)
#        return(score_df)
#
#    def get_average_patient_data(self, key = 'all'):
#        #generates a dictionary with *some* (positions, distances, and volumes) of the
#        #average data accross all patients
#        avg_centroids = np.zeros((45,3))
#        avg_volumes = np.zeros((45,))
#        avg_distances = np.zeros((45,45))
#        avg_tumor_distances = np.zeros((45,))
#        avg_doses = np.zeros((45,))
#        avg_tumor_volume = 0.0
#        for patient in self.get_patients():
#            avg_centroids += patient.centroids
#            avg_volumes += patient.volumes
#            avg_distances += patient.distances
#            avg_tumor_distances += patient.tumor_distances
#            avg_tumor_volume += patient.tumor_volume
#            avg_doses += patient.doses
#        p_avg = {}
#        p_avg['centroids'] = avg_centroids / self.num_patients
#        p_avg['volumes'] = avg_volumes / self.num_patients
#        p_avg['distances'] = avg_distances / self.num_patients
#        p_avg['tumor_distances'] = avg_tumor_distances / self.num_patients
#        p_avg['tumor_volume'] = avg_tumor_volume / self.num_patients
#        p_avg['doses'] = avg_doses / self.num_patients
#        #defaults to a dict, adding in a parameter to only look at one thing
#        if key == 'all':
#            return(p_avg)
#        else:
#            return(p_avg[key])
#
#    def evaluate(self, rank_function = 'tumor_organ_ssim', key = None, weights = np.array([0,1]), num_matches = 10):
#        #gives a bunch of different metrics for evaluating a given metric
#        estimates = self.predict_doses(weights, num_matches)
#        differences = self.doses - estimates
#        patient_mean_error = self.labeled_mean_error(differences, axis = 1)
#        organ_mean_error = self.labeled_mean_error(differences, axis = 0)
#        total_mean_error = np.mean(np.abs(differences))
#        percent_error = np.sum(np.abs(differences), axis = 1)/np.sum(self.doses, axis = 1)
#        total_rmse = np.sqrt(np.mean(differences**2))
#        result_dict = {'prediction': estimates,
#                       'patient_mean_error': patient_mean_error,
#                       'mean_error': total_mean_error,
#                       'rmse': total_rmse,
#                       'differences': differences,
#                       'organ_mean_error': organ_mean_error,
#                       'percent_error': percent_error}
#        if key is not None:
#            try:
#                return(result_dict[key])
#            except:
#                print('invalid key passed to PatientSet.evaluate()')
#                return(result_dict)
#        return(result_dict)
#
#    def labeled_mean_error(self, differences, axis):
#        #gives us a nice sorted list organ or patient total mean error as a labeled tuple
#        error = np.mean(np.abs(differences), axis = axis)
#        if axis == 0: #axis 0 is features, so organs here
#            labels = Constants.organ_list
#        else:
#            labels = self.ids #ids for the patients in sorted order?
#        name_error_tuples = [ (labels[x], error[x] ) for x in range(0, len(error))]
#        name_error_tuples = sorted(name_error_tuples, key = lambda x: x[1])
#        return(name_error_tuples)
#
#    def run_study(self, max_matches = 20,  weights = np.array([0,1])):
#        #tests out a metric for a range of difference potential matches and gives a minimum
#        error_hist = []
#        for num_matches in range(2, max_matches):
#            estimates = self.predict_doses(weights, num_matches)
#            error = np.mean(np.abs(self.doses - estimates))
#            error_hist.append(error)
#        print('min error of', min(error_hist), ' at ', np.argmin(error_hist) + 2)
#        return(error_hist)
#
#    def get_organ_clusters(self, max_dist):
#        #creates an index-based dictionary of k-nearest organs for each organ
#        #used in clusters
#        mean_dists = self.get_average_patient_data('distances')
#        organ_kmeans = {}
#        for organ_row in range(0, Constants.num_organs):
#            organ_dists = mean_dists[organ_row, :]
#            k = len(np.where( organ_dists < max_dist)[0])
#            organ_args = np.argsort(organ_dists)[0: k]
#            organ_kmeans[organ_row] = organ_args
#        return(organ_kmeans)
#
#    def gen_score_matrix(self, weights, classes = True):
#        #weights basically scales the importants of the knn cluster (one centered on each organ)
#        #currently not in use?
#        try:
#            if weights.shape != (Constants.num_organs,):
#                weights = np.ones((Constants.num_organs,))
#        except:
#            weights = np.ones((Constants.num_organs,))
#        scores = np.zeros((self.num_patients, self.num_patients))
#        #go compare each patient to every other in a nice for loop
#        for row in range(0, self.num_patients):
#            #skip comparing the patient to itself - will default to 0 so it doesn't appear in the comparison
#            for col in range(row + 1, self.num_patients):
#                p1 = self.patients[row]
#                p2 = self.patients[col]
#                #divide gorups up
#                if classes == True:
#                    if p1.group != p2.group:
#                        scores[row, col] = 0
#                        continue
#                score = []
#                #run the ssim or whatever for each knn cluster of organ-tumor distances
#                for key, value in self.organ_kmeans.items():
#                    d1 = p1.tumor_distances[value]
#                    d2 = p2.tumor_distances[value]
#                    v1 = p1.volumes[value]
#                    v2 = p2.volumes[value]
#                    #changing weights so now it's based on number of organs in common with tumor-organ overlap?
#                    weight = 1
#                    for x in range(0,len(value)):
#                        if d1[x] <= 0 and d2[x] <= 0:
#                            weight += 1/len(value)
#                    similarity = Rankings.local_ssim(d1,d2,v1,v2)
#                    if np.isnan(similarity):
#                        similarity = 0
#                    score.append(similarity*weight)
#                scores[row, col] = np.nanmean(score)
#        #scores is a semetric matrix, normalize from 0-1 in case I use something that's not the ssim
#        scores += np.transpose(scores)
#        scores = .999*(scores - scores.min()) / (scores.max() - scores.min())
#        return(scores)
#
#    def predict_doses(self,weights, num_matches = 5):
#        #generates an ndarray of dose estimates based on algorithm parameters
#        #rank_function and weights are the functions used to match dose
#        estimates = np.zeros(self.doses.shape)
#        ranks = self.gen_score_matrix(weights = weights)
#        for patient_idx in range(0, self.num_patients):
#            rank_row = ranks[patient_idx, :]
#            p = self.patients[patient_idx]
#            matches = self.get_num_matches(p)
#            estimates[patient_idx, :] = self.estimate_patient_doses(rank_row, matches)
#        return(estimates)
#
#    def get_num_matches(self, patient):
#        #function for determining the number of matches to use, so this can be changed easily
#        matches = 4
#        if patient.group == 1:
#            matches = 11
#        return matches
#
#    def get_matches(self, ranks, num_matches):
#        sorted_matches = np.argsort(-ranks)
#        top_matches = sorted_matches[0:num_matches]
#        return(top_matches)
#
#    def estimate_patient_doses(self, ranks, num_matches):
#            top_matches = self.get_matches(ranks, num_matches)
#            scores = ranks[top_matches]
#            matched_dosages = self.doses[tuple(top_matches), :]
#            #weight things by their scores
#            for match_idx in range(0, num_matches):
#                matched_dosages[match_idx, :] = scores[match_idx]*matched_dosages[match_idx, :]
#            if np.mean(scores) > 0:
#                patient_estimates = np.mean(matched_dosages, axis = 0)/np.mean(scores)
#            else:
#                patient_estimates = self.get_average_patient_data('doses')
#            return(patient_estimates)