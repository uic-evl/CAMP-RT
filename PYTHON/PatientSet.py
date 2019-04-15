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
        