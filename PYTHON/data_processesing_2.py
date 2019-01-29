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

class Constants():
    #all non-tumor organs being included (some datasets have more)
#    organ_list = [
#        'Brainstem','Cricoid_cartilage','Cricopharyngeal_Muscle',
#        'Esophagus','Extended_Oral_Cavity','Genioglossus_M',
#        #'Glottic_Area',
#        'Hard_Palate','Hyoid_bone',
#        'IPC','Larynx','Lower_Lip',
#        'Lt_Ant_Digastric_M','Lt_Anterior_Seg_Eyeball',
#        'Lt_Brachial_Plexus','Lt_Lateral_Pterygoid_M',
#        'Lt_Masseter_M','Lt_Mastoid',
#        'Lt_Medial_Pterygoid_M','Lt_Parotid_Gland',
#        'Lt_Posterior_Seg_Eyeball','Lt_Sternocleidomastoid_M',
#        'Lt_Submandibular_Gland','Lt_thyroid_lobe',
#        'Mandible','MPC','Mylogeniohyoid_M',
#        'Rt_Ant_Digastric_M','Rt_Anterior_Seg_Eyeball',
#        'Rt_Brachial_Plexus','Rt_Lateral_Pterygoid_M',
#        'Rt_Masseter_M','Rt_Mastoid',
#        'Rt_Medial_Pterygoid_M','Rt_Parotid_Gland',
#        'Rt_Posterior_Seg_Eyeball','Rt_Sternocleidomastoid_M',
#        'Rt_Submandibular_Gland','Rt_thyroid_lobe',
#        'Soft_Palate','SPC','Spinal_Cord',
#        'Supraglottic_Larynx','Thyroid_cartilage','Tongue',
#        'Upper_Lip'
#    ]
    organ_list = ['Extended_Oral_Cavity',
                  'Genioglossus_M',
                  'Hard_Palate',
                  'Lower_Lip',
                  'Lt_Ant_Digastric_M',
                  'Lt_Lateral_Pterygoid_M',
                  'Lt_Masseter_M',
                  'Lt_Medial_Pterygoid_M',
                  'Mandible',
                  'Mylogeniohyoid_M',
                  'Rt_Ant_Digastric_M',
                  'Rt_Lateral_Pterygoid_M',
                  'Rt_Masseter_M',
                  'Rt_Medial_Pterygoid_M',
                  'Soft_Palate',
                  'Tongue',
                  'Upper_Lip',
                  'Cricoid_cartilage',
                  'Cricopharyngeal_Muscle',
                  'Esophagus',
                  'Hyoid_bone',
                  'IPC',
                  'Larynx',
                  'Lt_Sternocleidomastoid_M',
                  'Lt_thyroid_lobe',
                  'MPC',
                  'Rt_Sternocleidomastoid_M',
                  'Rt_thyroid_lobe',
                  'SPC',
                  'Supraglottic_Larynx',
                  'Thyroid_cartilage',
                  'Lt_Parotid_Gland',
                  'Lt_Submandibular_Gland',
                  'Rt_Parotid_Gland',
                  'Rt_Submandibular_Gland',
                  'Lt_Anterior_Seg_Eyeball',
                  'Lt_Posterior_Seg_Eyeball',
                  'Rt_Anterior_Seg_Eyeball',
                  'Rt_Posterior_Seg_Eyeball',
                  'Brainstem',
                  'Lt_Brachial_Plexus',
                  'Rt_Brachial_Plexus',
                  'Spinal_Cord',
                  'Lt_Mastoid',
                  'Rt_Mastoid'
            ]
    num_organs = len(organ_list)
    #GTVn analgoues: 'GTV node', 'GTV-N', 'GTV_n', 'GTVn2', 'GTVn1'
    #GTVp analogues: 'GTV primary', 'GTV-P', 'GTV_p'
    tumor_aliases = {'GTV node': 'GTVn', 'GTV-N': 'GTVn',
                     'GTV_n': 'GTVn', 'GTVn1': 'GTVn',
                     'GTVn2': 'GTVn','GTV primary': 'GTVp',
                     'GTV-P': 'GTVp', 'GTV_p': 'GTVp',
                     'GTV nodes': 'GTVn', 'GTV tongue': 'GTVp',
                     'GTV-N1': 'GTVn', 'GTV-N2': 'GTVn'}
    #names to use for the dataframe when I read a centroid file.  it's cleaner
    centroid_file_names = ['ROI', 'x', 'y', 'z', 'mean_dose', 'volume', 'min_dose', 'max_dose']
    #patients without organs, using their sorted order (12th patient in the sorted list ect)
    missing_organs = {34: {'id': 131, 'organs': {'Lt_Parotid_Gland'}},
         55: {'id': 174, 'organs': {'Lt_Mastoid'}},
         56: {'id': 175, 'organs': {'Rt_Ant_Digastric_M'}},
         74: {'id': 205, 'organs': {'Lt_thyroid_lobe'}},
         100: {'id': 249, 'organs': {'Spinal_Cord','Supraglottic_Larynx','Thyroid_cartilage','Tongue','Upper_Lip'}},
         137: {'id': 5001, 'organs': {'Lt_Posterior_Seg_Eyeball'}},
         138: {'id': 5003, 'organs': {'Rt_Ant_Digastric_M'}},
         149: {'id': 5042,'organs': {'Lt_Posterior_Seg_Eyeball', 'Rt_Posterior_Seg_Eyeball'}},
         151: {'id': 5053,'organs': {'Lt_Posterior_Seg_Eyeball','Rt_Anterior_Seg_Eyeball','Rt_Posterior_Seg_Eyeball'}},
         168: {'id': 10009, 'organs': {'Rt_Lateral_Pterygoid_M'}},
         169: {'id': 10011, 'organs': {'Cricoid_cartilage'}},
         170: {'id': 10013, 'organs': {'Rt_Parotid_Gland'}},
         172: {'id': 10015, 'organs': {'Rt_Brachial_Plexus'}},
         173: {'id': 10018,'organs': {
            'Cricopharyngeal_Muscle',
           'Extended_Oral_Cavity',
           'Lt_Anterior_Seg_Eyeball',
           'Lt_Brachial_Plexus',
           'Lt_Posterior_Seg_Eyeball',
           'Mylogeniohyoid_M',
           'Rt_Anterior_Seg_Eyeball',
           'Rt_Brachial_Plexus',
           'Rt_Posterior_Seg_Eyeball',
           'Rt_Sternocleidomastoid_M',
           'Spinal_Cord',
           'Supraglottic_Larynx',
           'Thyroid_cartilage'}},
         175: {'id': 10020,'organs': {
            'Cricopharyngeal_Muscle',
           'Extended_Oral_Cavity',
           'Lt_Anterior_Seg_Eyeball',
           'Lt_Brachial_Plexus',
           'Lt_Posterior_Seg_Eyeball',
           'Mylogeniohyoid_M',
           'Rt_Anterior_Seg_Eyeball',
           'Rt_Brachial_Plexus',
           'Rt_Posterior_Seg_Eyeball',
           'Rt_Sternocleidomastoid_M',
           'Spinal_Cord',
           'Supraglottic_Larynx'}},
         177: {'id': 10022,'organs': {
            'Cricopharyngeal_Muscle',
           'Extended_Oral_Cavity',
           'Lt_Anterior_Seg_Eyeball',
           'Lt_Brachial_Plexus',
           'Lt_Posterior_Seg_Eyeball',
           'Mylogeniohyoid_M',
           'Rt_Anterior_Seg_Eyeball',
           'Rt_Brachial_Plexus',
           'Rt_Posterior_Seg_Eyeball',
           'Rt_Sternocleidomastoid_M',
           'Spinal_Cord',
           'Supraglottic_Larynx'}},
         178: {'id': 10029, 'organs': {'Rt_Ant_Digastric_M'}},
         184: {'id': 10044,'organs': {
            'Cricopharyngeal_Muscle',
           'Extended_Oral_Cavity',
           'Lt_Anterior_Seg_Eyeball',
           'Lt_Brachial_Plexus',
           'Lt_Posterior_Seg_Eyeball',
           'Mylogeniohyoid_M',
           'Rt_Anterior_Seg_Eyeball',
           'Rt_Brachial_Plexus',
           'Rt_Posterior_Seg_Eyeball',
           'Rt_Sternocleidomastoid_M',
           'Spinal_Cord',
           'Supraglottic_Larynx'}},
         190: {'id': 10065,'organs': {'Lt_Anterior_Seg_Eyeball', 'Rt_Anterior_Seg_Eyeball'}},
         199: {'id': 10085, 'organs': {'Lt_Parotid_Gland'}},
         203: {'id': 10094, 'organs': {'Hyoid_bone'}},
         205: {'id': 10103, 'organs': {'Hyoid_bone'}},
         211: {'id': 10130, 'organs': {'Brainstem', 'Cricoid_cartilage'}},
         213: {'id': 10132, 'organs': {'Thyroid_cartilage'}},
         219: {'id': 10140, 'organs': {'Rt_Anterior_Seg_Eyeball'}},
         230: {'id': 10157, 'organs': {'Rt_Posterior_Seg_Eyeball'}},
         236: {'id': 10176, 'organs': {'Thyroid_cartilage'}},
         241: {'id': 10191, 'organs': {'Rt_Brachial_Plexus'}}}
    no_gtvp = []
    no_gtvn = []

class Rankings():
    #ranking functions that generate a score, takes in pateint objects
    old_weights = np.array([3432,32423,1])

    def vector_ssim(p1, p2):
        upper_triangle = np.triu_indices(len((p1.distances)))
        d1 = p1.distances[upper_triangle].ravel()
        d2 = p2.distances[upper_triangle].ravel()
        return(compare_ssim(d1, d2, win_size = 5))

    def ssim(p1, p2):
        return(compare_ssim(p1.distances, p2.distances, win_size = 3))

    def ssim_with_laterality(p1, p2, weights = np.array([1,2,.05])):
        scores = np.zeros((3,))
        scores[0] = Rankings.ssim(p1,p2)
        scores[1] = 1 if p1.laterality == p2.laterality else 0
        scores[2] = 1 - np.abs(p1.gtvp_volume - p2.gtvp_volume)/(
                np.max([p1.gtvp_volume, p2.gtvp_volume]) + .000001)
        final_score = np.sum(scores*weights)/np.mean(weights)
        return(final_score)

    def mse(p1,p2):
        return( 1/( compare_mse(p1.distances, p2.distances) + .000001) )

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

    def experimental(p1, p2, weights = np.array([1,.4,1,.05,.1])):
        if not np.array_equal(Rankings.old_weights, weights):
            Rankings.old_weights = weights
            print('weights ', weights)
        scores = np.zeros((5,))
        scores[0] = compare_ssim(p1.gtvp_dists, p2.gtvp_dists, win_size = 3)
        scores[1] = Rankings.vector_ssim(p1, p2)
        scores[2] = 1 if p1.laterality == p2.laterality else 0
        scores[3] = 1 - np.abs(p1.gtvp_volume - p2.gtvp_volume)/(
                p1.gtvp_volume + p2.gtvp_volume + .000001)
        scores[4] = 1/(compare_mse(p1.volumes, p2.volumes) + .000001)
        final_score = np.sum(scores*weights)/np.mean(weights)
        return(final_score)


class Patient():

    def __init__(self, distances, doses, p_id, position, info):
        #patient ID number
        self.id = p_id
        #basically ordinality of the id, so where it will be in an index
        self.pos = position
        self.check_missing_organs(distances, doses)
        self.laterality = info['Tm Laterality (R/L)']
        self.age = info['Age at Diagnosis (Calculated)']
        print(len(distances.columns), ' ', self.id)
        centroid_data = self.get_doses_file_info(doses)
        self.doses = centroid_data[:, 4]
        self.volumes = centroid_data[:, 3]
        self.centroids = centroid_data[:, 0:3]
        self.distances = self.gen_distance_matrix(distances)
        (self.gtvp_dists, self.gtvn_dists) = self.get_tumor_distances(distances)
        #store the entries without gtvp for future study
        if(self.gtvp_volume == 0):
            Constants.no_gtvp.append(self.id)
        if(self.gtvn_volume == 0):
            Constants.no_gtvn.append(self.id)

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
            self.gtvp_position = gtvp[['x','y','z']].values
        except:
            self.gtvp_volume = 0
            self.gtvp_position = np.array([0,0,0])
        #extract a secondary tumor (only gets the first one?)
        #several patients have no gtvp but a gtvn
        try:
            gtvn = centroids.loc['GTVn']
            if(len(gtvn.volume) > 1): #there are often multiple gtvns, how do I handle that
                gtvn = gtvn.iloc[0]
            self.gtvn_volume = gtvn.volume
            self.gtvn_position = gtvn[['x','y','z']].values
        except:
            self.gtvn_volume = 0
            self.gtvn_position = np.array([0,0,0])
        #get the info the centers, volumes, nad doses for all the things
        centroid_matrix = np.zeros((Constants.num_organs,5)) #row = x,y,z,volume,dose
        for idx in range(0, Constants.num_organs):
            organ = Constants.organ_list[idx]
            try:
                organ_entry = centroids.loc[organ]
                centroid_matrix[idx, 0:3] = organ_entry[['x','y','z']].values[0]
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
            try:
                tumor_row = dists.loc['GTVp', Constants.organ_list[idx]]
                gtvp_dists[idx] = tumor_row['Eucledian Distance (mm)']
            except:
                gtvp_dists[idx] = -1
            try:
                tumor_row = dists.loc['GTVn', Constants.organ_list[idx]]
                gtvn_dists[idx] = tumor_row['Eucledian Distance (mm)']
            except:
                gtvn_dists[idx] = -1
        return((gtvp_dists, gtvn_dists))

class PatientSet():

    def __init__(self, outliers = []):
        self.outliers = outliers
        (self.patients, self.doses, self.num_patients) = self.read_patient_data()
        print('\npatient data loaded...\n')

    def read_patient_data(self):
        #sorts by size of largest integer string, which is the id for our files
        file_sort = lambda x: sorted(x, key =
                                     lambda file:
                                         max([int(x) for x in findall("[0-9]+", file)])
                                )
        distance_files = file_sort(glob('patient_files\\' + '**/*distances.csv'))
        dose_files = file_sort(glob('patient_files\\' + '**/*centroid*.csv'))
        for file_pos in sorted(self.outliers, reverse = True):
            del distance_files[file_pos]
            del dose_files[file_pos]
        metadata_file = 'data\\patient_info.csv'
        assert(len(distance_files) == len(dose_files))
        #maps a position 0-len(files) to the dummy id for a patient
        ids = [max([int(x) for x in findall('[0-9]+', file)]) for file in distance_files]
        num_patients = len(ids)
        metadata = pd.read_csv(metadata_file,
                               index_col = 0, #index is the "Dummy ID"
                               usecols = [0,1,2,3,4,5,6,7,8,9,10,11]
                               ).loc[ids]
        patients = OrderedDict()
        dose_matrix = np.zeros((num_patients, len(Constants.organ_list)))
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
        return((patients, dose_matrix, num_patients))

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
            for match_idx in range(0, num_matches):
                matched_dosages[match_idx, :] = scores[match_idx]*matched_dosages[match_idx, :]
            estimates[patient_idx, :] = np.mean(matched_dosages, axis = 0)/np.mean(scores)
        return(estimates)

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

db = PatientSet()
pickle.dump(db, open('data\\patient_data.p', 'wb'))
#db = pickle.load(open('data\\patient_data.p', 'rb'))
#weights = np.array([1,2,.05])
#test_weights = np.array([2,0.5,2,.05,1])
#max_count = 20
#
#from scipy.optimize import minimize
#func = lambda x: max(db.run_study(rank_function = 'experimental', weights = x, max_matches = 8))
#result = minimize(func, test_weights, method = 'CG', options ={'disp': True, 'eps': 1})
#
#rand_hist = db.run_study(rank_function = 'random', max_matches = max_count, weights = 1)
#ssim_hist = db.run_study(rank_function = 'ssim', weights = weights, max_matches = max_count)
#mse_hist = db.run_study(rank_function = 'experimental', weights = test_weights, max_matches = max_count)
#min_error_hist = db.run_study(rank_function = 'min_dose_error', max_matches = max_count, weights = 1)
#
#x = list(range(2, max_count))
#plt.plot(x, rand_hist, x, ssim_hist, x, mse_hist, x, min_error_hist)
##plt.plot(x, rand_hist, x, mse_hist, x, min_error_hist)
#plt.legend(['random','ssim','test function','min error'])