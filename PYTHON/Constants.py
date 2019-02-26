# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 14:47:15 2019

@author: Andrew
"""

class Constants():
    #all non-tumor organs being included (some datasets have more)
    #currently soreted based on the organ atals
#    organ_list = ['Extended_Oral_Cavity',
#                  'Genioglossus_M',
#                  'Hard_Palate',
#                  'Lower_Lip',
#                  'Lt_Ant_Digastric_M',
#                  'Lt_Lateral_Pterygoid_M',
#                  'Lt_Masseter_M',
#                  'Lt_Medial_Pterygoid_M',
#                  'Mandible',
#                  'Mylogeniohyoid_M',
#                  'Rt_Ant_Digastric_M',
#                  'Rt_Lateral_Pterygoid_M',
#                  'Rt_Masseter_M',
#                  'Rt_Medial_Pterygoid_M',
#                  'Soft_Palate',
#                  'Tongue',
#                  'Upper_Lip',
#                  'Cricoid_cartilage',
#                  'Cricopharyngeal_Muscle',
#                  'Esophagus',
#                  'Hyoid_bone',
#                  'IPC',
#                  'Larynx',
#                  'Lt_Sternocleidomastoid_M',
#                  'Lt_thyroid_lobe',
#                  'MPC',
#                  'Rt_Sternocleidomastoid_M',
#                  'Rt_thyroid_lobe',
#                  'SPC',
#                  'Supraglottic_Larynx',
#                  'Thyroid_cartilage',
#                  'Lt_Parotid_Gland',
#                  'Lt_Submandibular_Gland',
#                  'Rt_Parotid_Gland',
#                  'Rt_Submandibular_Gland',
#                  'Lt_Anterior_Seg_Eyeball',
#                  'Lt_Posterior_Seg_Eyeball',
#                  'Rt_Anterior_Seg_Eyeball',
#                  'Rt_Posterior_Seg_Eyeball',
#                  'Brainstem',
#                  'Lt_Brachial_Plexus',
#                  'Rt_Brachial_Plexus',
#                  'Spinal_Cord',
#                  'Lt_Mastoid',
#                  'Rt_Mastoid'
#            ]
    #organ list set via affinity propogation clustering
    organ_list = ['Esophagus',
         'Spinal_Cord',
         'Lt_Brachial_Plexus',
         'Rt_Brachial_Plexus',
         'Cricopharyngeal_Muscle',
         'Lt_thyroid_lobe',
         'Rt_thyroid_lobe',
         'Cricoid_cartilage',
         'IPC',
         'MPC',
         'Brainstem',
         'Larynx',
         'Thyroid_cartilage',
         'Rt_Sternocleidomastoid_M',
         'Rt_Mastoid',
         'Rt_Parotid_Gland',
         'Rt_Medial_Pterygoid_M',
         'Rt_Lateral_Pterygoid_M',
         'Rt_Masseter_M',
         'Lt_Sternocleidomastoid_M',
         'Lt_Mastoid',
         'Lt_Parotid_Gland',
         'Lt_Submandibular_Gland',
         'Lt_Medial_Pterygoid_M',
         'Lt_Lateral_Pterygoid_M',
         'Lt_Masseter_M',
         'Supraglottic_Larynx',
         'SPC',
         'Rt_Submandibular_Gland',
         'Hyoid_bone',
         'Soft_Palate',
         'Genioglossus_M',
         'Tongue',
         'Rt_Ant_Digastric_M',
         'Lt_Ant_Digastric_M',
         'Mylogeniohyoid_M',
         'Extended_Oral_Cavity',
         'Mandible',
         'Hard_Palate',
         'Lt_Posterior_Seg_Eyeball',
         'Rt_Posterior_Seg_Eyeball',
         'Lt_Anterior_Seg_Eyeball',
         'Rt_Anterior_Seg_Eyeball',
         'Lower_Lip',
         'Upper_Lip']

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
                     'GTV_P': 'GTVp',
                     'GTV P': 'GTVp',
                     'GTV nodes': 'GTVn',
                     #'GTV tongue': 'GTVp',
                     'GTV-N1': 'GTVn',
                     'GTV N': 'GTVn',
                     'GTV-NR': 'GTVn' #there is also an NL?
                     #'GTV-N2': 'GTVn'
                     }
    v3_bad_entries = [11, 33, 37, 55, 56, 67, 69, 72, 73, 74, 76, 77, 79,
                      81, 17, 39, 101, 115, 2020, 10011, 10034,10043,
                      10074,10080, 10094, 10145, 10148, 10164, 10174, 
                      10181, 10086, 10147, #missing oragans
                      131, 5001, 5003, 5042, 5053, 10009, 10013, 10015, 
                      10018,10020, 10022, 10029, 10044, 10046, 10065, 10070, 
                      10075, 10077, 10085, 10103, 10130, 10132, 10140, 10157, 10176, 10191,
                      10086, 10108, 10131, 10143, 10159, 10181, 10193]
    v2_bad_entries = [239, 10034, 10164, 174, 175, 205, 249, 2009, 5059]
    v2_half_dosed = [160, 5077, 178, 155, 221, 251, 2007, 212, 234,154]
    #pateints with like, really high neck radation, chosen manually
    v2_high_throat_dose = [184, 276, 5084, 265, 246, 289, 283, 5068, 197]
    #names to use for the dataframe when I read a centroid file.  it's cleaner
    centroid_file_names_v2 = ['ROI', 'x', 'y', 'z', 'mean_dose', 'volume', 'min_dose', 'max_dose']
    centroid_file_names_v3 = ['ROI', 'x', 'y', 'z', 'volume', 'min_dose', 'mean_dose', 'max_dose']
    #patients without organs, using their sorted order (12th patient in the sorted list ect)
    missing_organs = {}
    no_tumor = []