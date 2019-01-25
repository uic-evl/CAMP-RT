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

    #patients without organs, using their sorted order (12th patient in the sorted list ect)
    patients_without_organs = [12,13,30,55,69]

class Patient():

    def __init__(self, distances, doses, p_id, position, metadata):
        #patient ID number
        self.id = p_id
        #basically ordinality of the id, so where it will be in an index
        self.pos = position
        self.distances = self.gen_distance_matrix(distances)
        self.doses = doses

    def gen_distance_matrix(self, distances):
        #extracts a 45x45 matrix of the distancess between organs
        #should fill missing values with 0? doesn't include tumors since those aren't always there
        #if an entry is gone idk what to do
        distances = distances.loc[(
                distances['Reference ROI'].isin(Constants.organ_list)) & (
                distances['Target ROI'].isin(Constants.organ_list))]
        distances = self.fill_missing_organs(distances)
        distances = distances.set_index(['Reference ROI', 'Target ROI'])
        distances.sort_values(by = ['Reference ROI', 'Target ROI'], inplace = True)
        #extracts the matrix of values.  set_index sorts by default
        dist_matrix = distances['Eucledian Distance (mm)'].unstack(fill_value = 0).values
        #dist_matrix += np.transpose(dist_matrix)
        if(len(dist_matrix) != 45):
            print('Patient with id ', self.pos, ' is missing organs')
        #so this result is a 45x45 upper triangular.
        #values are more like the upper triangle of a 46x46 matrix with trace = 0
        #I don't know what to do with that information, so this is just arbi
        good_indices = np.triu_indices(len(dist_matrix))
        dist_matrix = dist_matrix[good_indices]
        return(dist_matrix)

    def fill_missing_organs(self, distance_df):
        #beautiful, understandable code that gets the organs neither in 'Reference ROI'
        #nor 'Target ROI' for the organ master list
        missing_organs = set(Constants.organ_list) - set(distance_df['Reference ROI'].append(distance_df['Target ROI']).unique())
        #I feel like this might not work if the first or last organ is missing, but that hasn't happened yet
        for organ in missing_organs:
            print(self.pos, 'is getting ',organ, 'added')
            distance_df = distance_df.append({'Reference ROI': organ,
                                              'Target ROI': organ,
                                              'Eucledian Distance (mm)': 0},
                                             ignore_index = True)
        return(distance_df)

#sorts by size of largest integer string, which is the id for our files
file_sort = lambda x: sorted(x, key =
                             lambda file:
                                 max([int(x) for x in findall("[0-9]+", file)])
                        )
distance_files = file_sort(glob('patients_v2\\' + '**/*distances.csv'))
dose_files = file_sort(glob('patients_v2\\' + '**/*centroid*.csv'))
metadata_file = 'data\\patient_info.csv'
assert(len(distance_files) == len(dose_files))

#maps a position 0-len(files) to the dummy id for a patient
ids = [max([int(x) for x in findall('[0-9]+', file)]) for file in distance_files]
metadata = pd.read_csv(metadata_file,
                       index_col = 0, #index is the "Dummy ID"
                       usecols = [0,1,2,3,4,5,6,7,8,9,10,11]
                       ).loc[ids]
patients = OrderedDict()
#putting all the data into a patient object for further objectification
for patient_index in range(0, len(distance_files)):
    #these are indexed by name of organ
    distances = pd.read_csv(distance_files[patient_index]).replace(Constants.tumor_aliases)
    doses = pd.read_csv(dose_files[patient_index]).replace(Constants.tumor_aliases)
    info = metadata.loc[ids[patient_index]]
    new_patient = Patient(distances, doses,
                          ids[patient_index], patient_index, info)
    patients[patient_index] = new_patient