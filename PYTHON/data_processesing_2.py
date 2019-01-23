# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 17:22:31 2019

@author: Andrew
"""
from glob import glob
from re import findall
import numpy as np
import pandas as pd

ORGAN_LIST = [
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
    'Upper_Lip','GTVp','GTVn'
]

#sorts by size of largest integer string, which is the id for our files
file_sort = lambda x: sorted(x, key =
                             lambda file: max([int(x) for x in findall("[0-9]+", file)])
                        )
distance_files = file_sort(glob('patients_v2\\' + '**/*distances.csv'))
dose_files = file_sort(glob('patients_v2\\' + '**/*meandoses.csv'))
#this would be faster if I zipped them maybe?
ids = [max([int(x) for x in findall('[0-9]+', file)]) for file in distance_files]

patient_distances = pd.read_csv(distance_files[37])
patient_dose = pd.read_csv(dose_files[37])