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

class Constants():
    #all non-tumor organs being included (some datasets have more)
    organ_list = [
        'Brainstem','Cricoid_cartilage','Cricopharyngeal_Muscle',
        'Esophagus','Extended_Oral_Cavity','Genioglossus_M',
        'Glottic_Area','Hard_Palate','Hyoid_bone',
        'IPC','Larynx','Lower_Lip',
        'Lt_Ant_Digastric_M','Lt_Anterior_Seg_Eyeball',
        'Lt_Brachial_Plexus','Lt_Lateral_Pterygoid_M',
        'Lt_Masseter_M','Lt_Mastoid',
        'Lt_Medial_Pterygoid_M','Lt_Parotid_Gland',
        'Lt_Posterior_Seg_Eyeball','Lt_Sternocleidomastoid_M',
        'Lt_Submandibular_Gland','Lt_thyroid_lobe',
        'Mandible','MPC','Mylogeniohyoid_M',
        'Rt_Ant_Digastric_M','Rt_Anterior_Seg_Eyeball',
        'Rt_Brachial_Plexus','Rt_Lateral_Pterygoid_M',
        'Rt_Masseter_M','Rt_Mastoid',
        'Rt_Medial_Pterygoid_M','Rt_Parotid_Gland',
        'Rt_Posterior_Seg_Eyeball','Rt_Sternocleidomastoid_M',
        'Rt_Submandibular_Gland','Rt_thyroid_lobe',
        'Soft_Palate','SPC','Spinal_Cord',
        'Supraglottic_Larynx','Thyroid_cartilage','Tongue',
        'Upper_Lip'
    ]
    #GTVn analgoues: 'GTV node', 'GTV-N', 'GTV_n', 'GTVn2', 'GTVn1'
    #GTVp analogues: 'GTV primary', 'GTV-P', 'GTV_p'
    tumor_aliases = {'GTV node': 'GTVn', 'GTV-N': 'GTVn',
                     'GTV_n': 'GTVn', 'GTVn1': 'GTVn',
                     'GTVn2': 'GTVn','GTV primary': 'GTVp',
                     'GTV-P': 'GTVp', 'GTV_p': 'GTVp'}
    #names to use for the dataframe when I read a centroid file.  it's cleaner
    centroid_file_names = ['ROI', 'x', 'y', 'z', 'mean_dose', 'volume', 'min_dose', 'max_dose']
    #patients without organs, using their sorted order (12th patient in the sorted list ect)
    patients_without_organs = [12,13,30,55,69]
    
class Rankings():
    #ranking functions that generate a score, takes in pateint objects
    def ssim(p1, p2):
        return(compare_ssim(p1.distances, p2.distances))
    
    def emd(patient_1, patient_2):
        #simplified earth movers distance - here it's just work done to move organs in the less massive
        #patient into the same positions of the more massive patient
        p_ref = patient_1
        p_new = patient_2
        if(np.sum(patient_1.volumes) > np.sum(patient_2.volumes)):
            p_ref = patient_2
            p_new = patient_1
        volumes = p_ref.volumes
        dists = np.sqrt(np.sum(np.abs(p_ref.centroids**2 - p_new.centroids**2), axis = 1))
        work = np.sum(dists*volumes)
        #we this to work with the other functions so >0 and high scores = closer
        return(1/(work + .000001))
    
class Patient():

    def __init__(self, distances, doses, p_id, position, info):
        #patient ID number
        self.id = p_id
        #basically ordinality of the id, so where it will be in an index
        self.pos = position
        self.laterality = info['Tm Laterality (R/L)']
        self.age = info['Age at Diagnosis (Calculated)']
        self.distances = self.gen_distance_matrix(distances)
        centroid_data = self.get_doses_file_info(doses)
        self.doses = centroid_data.mean_dose.values
        self.volumes = centroid_data.volume.values
        self.centroids = centroid_data[['x','y','z']].values
        
        
    def get_doses_file_info(self, doses):
        doses.columns = Constants.centroid_file_names
        doses = doses.sort_values(by = ['ROI'])
        centroids = self.center_doses(doses)
        self.gtvp = doses.loc[doses['ROI'] == 'GTVp']
        self.gtvn = doses.loc[doses['ROI'] == 'GTVn']
        if self.gtvp.empty:
            print('patient ', self.id, 'is missing gtvp')
            self.gtvp.volume = 0
            print(self.gtvp.volume)
        if self.gtvn.empty:
            self.gtvn.volume = 0
        centroids = centroids.loc[doses['ROI'].isin(Constants.organ_list)]
        #extracts a vector of the mean dose values with zero for missing organs added in
        missing_organs = set(Constants.organ_list) - set(centroids['ROI'].unique())
        #in the future I should get a mean value to replace here since 0 is a bad number
        for organ in missing_organs:
            #print(self.id, 'is getting ',organ, 'added')
            centroids = centroids.append({'ROI': organ, 'mean_dose': 0, 'x': 0, 'y': 0, 'z': 0, 'volume': 0}, ignore_index = True)
        centroids = centroids.sort_values(by = 'ROI', kind = 'mergesort')
        return(centroids)

    def center_doses(self, doses):
        #subtract off the mean so the pointcloud is centered at 0
        #should I just use a reference organ instead?  or rotate?
        doses.x -= doses.x.mean()
        doses.y -= doses.y.mean()
        doses.z -= doses.z.mean()
        return(doses)

    def gen_distance_matrix(self, distances):
        #extracts a 45x45 matrix of the distancess between organs
        #should fill missing values with 0? doesn't include tumors since those aren't always there
        distances = distances.loc[(
                distances['Reference ROI'].isin(Constants.organ_list)) & (
                distances['Target ROI'].isin(Constants.organ_list))]
        distances = self.fill_missing_organs(distances)
        distances = distances.sort_values(by = ['Reference ROI', 'Target ROI'], kind = 'mergesort')
        distances = distances.set_index(['Reference ROI', 'Target ROI'])
        #extracts the matrix of values.  set_index sorts by default
        dist_matrix = distances['Eucledian Distance (mm)'].unstack(fill_value = 0).values
        #dist_matrix += np.transpose(dist_matrix)
        if(len(dist_matrix) != 45):
            print('Patient with id ', self.pos, ' is missing organs')
        #so this result is a 45x45 upper triangular.
        #values are more like the upper triangle of a 46x46 matrix with trace = 0
        good_indices = np.triu_indices(len(dist_matrix))
        #unraveled vector for now?
        dist_matrix = dist_matrix[good_indices]
        return(dist_matrix)

    def fill_missing_organs(self, distance_df):
        #beautiful, understandable code that gets the organs neither in 'Reference ROI'
        #nor 'Target ROI' for the organ master list
        missing_organs = set(Constants.organ_list) - set(distance_df['Reference ROI'].append(distance_df['Target ROI']).unique())
        #I feel like this might not work if the first or last organ is missing, but that hasn't happened yet
        for organ in missing_organs:
            distance_df = distance_df.append({'Reference ROI': organ,
                                              'Target ROI': organ,
                                              'Eucledian Distance (mm)': 0},
                                             ignore_index = True)
        return(distance_df)
    
class PatientSet():
    
    def __init__(self):
        (self.patients, self.doses, self.num_patients) = self.read_patient_data()
        
    def read_patient_data(self):
        #sorts by size of largest integer string, which is the id for our files
        file_sort = lambda x: sorted(x, key =
                                     lambda file:
                                         max([int(x) for x in findall("[0-9]+", file)])
                                )
        distance_files = file_sort(glob('patients_v2\\' + '**/*distances.csv'))
        dose_files = file_sort(glob('patients_v2\\' + '**/*centroid*.csv'))
        metadata_file = 'data\\patient_info.csv'
        assert(len(distance_files) == len(dose_files))
        num_patients = len(distance_files)
        #maps a position 0-len(files) to the dummy id for a patient
        ids = [max([int(x) for x in findall('[0-9]+', file)]) for file in distance_files]
        metadata = pd.read_csv(metadata_file,
                               index_col = 0, #index is the "Dummy ID"
                               usecols = [0,1,2,3,4,5,6,7,8,9,10,11]
                               ).loc[ids]
        patients = OrderedDict()
        dose_matrix = np.zeros((num_patients, len(Constants.organ_list)))
        #putting all the data into a patient object for further objectification
        for patient_index in range(0, num_patients):
            #these are indexed by name of organ
            distances = pd.read_csv(distance_files[patient_index])
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
    
    def gen_score_matrix(self, score_function = 'ssim'):
        #generates a score matrix based on a rank function
        #function should rank more similar people with a higher number
        scores = np.zeros((self.num_patients, self.num_patients))
        for row in range(0, self.num_patients):
            for col in range(row + 1, self.num_patients):
                scores[row, col] = self.compare_traits(self.patients[row], self.patients[col],
                      rank_function = score_function, weights = np.array([1,.1,1]))
        #formats it as a symetric matrix with a zero diagonal
        scores += scores.transpose()
        #basically normalize the score so the max is 1?
        scores = scores/scores.max()
        return(scores)
    
    def compare_traits(self, p1, p2, rank_function = 'ssim', weights = 1):
        score = np.zeros((3,))
        if rank_function == 'ssim':
            score[0] = Rankings.ssim(p1, p1)
        elif rank_function == 'emd':
            score[0] = Rankings.emd(p1, p2)
        else:
            print('error, invalid rank method: ', rank_function)
        score[1] = 1 if p1.laterality == p2.laterality else 0
        #try making this also transformation distance?
        score[2] = 1 - np.abs(p1.gtvp.get('volume') - p2.gtvp.get('volume'))/(
                np.max([p1.gtvp.get('volume'), p2.gtvp.get('volume')]) + .000001)
        #gtvn for later?
        #score[3] = 1 - np.abs(p1.gtvn.volume - p2.gtvn.volume)/(max([p1.gtvn.volume, p2.gtvn.volume]) + .000001)
        final_score = np.sum(score*weights)/np.mean(weights)
        return(final_score)
    

db = PatientSet()
scores = db.gen_score_matrix(Rankings.emd)