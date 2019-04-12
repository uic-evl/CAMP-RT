# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 14:47:15 2019

@author: Andrew
"""

class Constants():
    #all non-tumor organs being included (some datasets have more)
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
                     'GTV primary': 'GTVp',
                     'GTV-P': 'GTVp',
                     'GTV_p': 'GTVp',
                     'GTV_P': 'GTVp',
                     'GTV P': 'GTVp',
                     'GTV nodes': 'GTVn',
                     'GTV-N1': 'GTVn',
                     'GTV N': 'GTVn',
                     'GTV-NR': 'GTVn2', #I am only aware of this for 10144, may need more robust solution later
                     'GTV-NL': 'GTVn3'
                     }
    v3_bad_entries = list(set([131, 5001, 5003, 5042, 5053, 10009, 10011, 10013, 10015, 10018, 10020, 10022, 10029, 10044, 10065, 10075, 10085, 10086, 10094, 10103, 10130, 10132, 10140, 10157, 10176,10191, 17, 39, 101, 115, 2020, 10011, 10043, 10074, 10080, 10094, 10145, 10148, 10174, 10181, 10018, 10019, 10022, 10041, 10044, 10129, 31, 37, 101, 105, 112, 115, 121, 133, 165, 2000, 2016, 5025, 10019, 10021, 10040, 10046, 10113, 10147, 10148, 10154, 3, 10, 112, 131, 145, 2000, 2016, 5041, 5043, 10021, 10040, 10041, 10044, 10046, 10071, 10080, 10103, 10129, 10130, 10138, 10143, 10154, 10155, 10181, 10184, 10191, 10193, 10197, 10197, 46, 121, 5039, 10020, 10044, 10022, 10018, 5041, 131, 10021, 10071, 10184, 10065]))
    v3_real_bad_entries = list(set([10138, 10103, 10086, 10147, 2000, 3, 10191]))
    v3_no_tumor = [17, 39, 101, 115, 2020, 10011, 10043, 10074, 10080, 10094, 10145, 10148, 10174, 10181]
    v3_bad_positions = [10018, 10022, 10020, 10044, 5041, 5309]
    v2_bad_entries = [239, 10034, 10164, 174, 175, 205, 249, 2009, 5059]
    v2_half_dosed = [160, 5077, 178, 155, 221, 251, 2007, 212, 234,154]
    #pateints with like, really high neck radation, chosen manually
    v2_high_throat_dose = [184, 276, 5084, 246, 289, 283, 5068, 197, 5007, 100]
    #names to use for the dataframe when I read a centroid file.  it's cleaner
    centroid_file_names_v2 = ['ROI', 'x', 'y', 'z', 'mean_dose', 'volume', 'min_dose', 'max_dose']
    centroid_file_names_v3 = ['ROI', 'x', 'y', 'z', 'volume', 'min_dose', 'mean_dose', 'max_dose']
    #patients without organs, using their sorted order (12th patient in the sorted list ect)
    missing_organs = {}
    no_tumor = []