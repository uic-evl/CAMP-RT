# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 14:47:15 2019

@author: Andrew
"""

class Constants():
    #all non-tumor organs being included (some datasets have more)
    #currently soreted based on the organ atals
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
    #organ lists sorted by botton to top or right to left incase these are better
    organ_list_ud_sorted = ['Lt_Anterior_Seg_Eyeball','Rt_Anterior_Seg_Eyeball',
     'Lt_Posterior_Seg_Eyeball','Rt_Posterior_Seg_Eyeball','Brainstem',
     'Lt_Lateral_Pterygoid_M','Rt_Lateral_Pterygoid_M','Hard_Palate','Upper_Lip','Rt_Masseter_M',
     'Lt_Mastoid','Lt_Masseter_M','Rt_Mastoid','Soft_Palate','Rt_Medial_Pterygoid_M',
     'Lt_Medial_Pterygoid_M','Rt_Parotid_Gland','Lt_Parotid_Gland',
     'Tongue','Extended_Oral_Cavity','Mandible','Lower_Lip',
     'SPC','Genioglossus_M','Mylogeniohyoid_M','Rt_Submandibular_Gland',
     'Lt_Submandibular_Gland','Lt_Ant_Digastric_M','Rt_Ant_Digastric_M',
     'Hyoid_bone','MPC','Supraglottic_Larynx','Spinal_Cord',
     'Lt_Sternocleidomastoid_M','Rt_Sternocleidomastoid_M','IPC',
     'Thyroid_cartilage','Larynx','Lt_Brachial_Plexus','Cricopharyngeal_Muscle',
     'Rt_Brachial_Plexus','Cricoid_cartilage','Rt_thyroid_lobe',
     'Lt_thyroid_lobe','Esophagus']

    organ_list_lr_sorted = ['Rt_Parotid_Gland','Rt_Mastoid','Rt_Masseter_M','Rt_Sternocleidomastoid_M',
     'Rt_Lateral_Pterygoid_M','Rt_Posterior_Seg_Eyeball','Rt_Submandibular_Gland',
     'Rt_Anterior_Seg_Eyeball','Rt_Medial_Pterygoid_M','Rt_Brachial_Plexus',
     'Rt_thyroid_lobe','Rt_Ant_Digastric_M','Brainstem','Hyoid_bone','Spinal_Cord',
     'Hard_Palate','Mandible','Thyroid_cartilage','Genioglossus_M', 'Extended_Oral_Cavity',
     'Soft_Palate','Tongue','SPC','Larynx','Cricoid_cartilage','Upper_Lip',
     'Mylogeniohyoid_M','MPC','Supraglottic_Larynx','IPC','Cricopharyngeal_Muscle',
     'Lower_Lip','Esophagus','Lt_Ant_Digastric_M','Lt_thyroid_lobe','Lt_Brachial_Plexus',
     'Lt_Medial_Pterygoid_M','Lt_Posterior_Seg_Eyeball','Lt_Submandibular_Gland',
     'Lt_Anterior_Seg_Eyeball','Lt_Lateral_Pterygoid_M','Lt_Sternocleidomastoid_M',
     'Lt_Masseter_M','Lt_Mastoid','Lt_Parotid_Gland']

    num_organs = len(organ_list)
    #GTVn analgoues: 'GTV node', 'GTV-N', 'GTV_n', 'GTVn2', 'GTVn1'
    #GTVp analogues: 'GTV primary', 'GTV-P', 'GTV_p'
    tumor_aliases = {'GTV node': 'GTVn', 
                     'GTV-N': 'GTVn',
                     'GTV_n': 'GTVn', 
                     'GTVn1': 'GTVn',
                     #'GTVn2': 'GTVn', #this one is causing issues with getting distance data
                     'GTV primary': 'GTVp',
                     'GTV-P': 'GTVp', 
                     'GTV_p': 'GTVp',
                     'GTV nodes': 'GTVn', 
                     #'GTV tongue': 'GTVp',
                     'GTV-N1': 'GTVn', 
                     'GTV-N2': 'GTVn'}
    #names to use for the dataframe when I read a centroid file.  it's cleaner
    centroid_file_names = ['ROI', 'x', 'y', 'z', 'mean_dose', 'volume', 'min_dose', 'max_dose']
    #patients without organs, using their sorted order (12th patient in the sorted list ect)
    missing_organs = {}
    no_tumor = []
    multiple_gtvn = []