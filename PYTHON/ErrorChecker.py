# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 12:07:35 2019

@author: Andrew
"""
import numpy as np
from Constants import Constants
from collections import namedtuple

OrganTuple = namedtuple('OrganTuple', ['organ', 'dose'])
class ErrorChecker():

    def __init__(self):
        self.ra_eye = Constants.organ_list.index('Rt_Anterior_Seg_Eyeball')
        self.rp_eye = Constants.organ_list.index('Rt_Posterior_Seg_Eyeball')
        self.la_eye = Constants.organ_list.index('Lt_Anterior_Seg_Eyeball')
        self.lp_eye = Constants.organ_list.index('Lt_Posterior_Seg_Eyeball')
        self.eyes = [self.ra_eye, self.rp_eye, self.la_eye, self.lp_eye]
        self.brainstem = Constants.organ_list.index('Brainstem')
        self.spinal_cord = Constants.organ_list.index('Spinal_Cord')
        self.eye_threshold = 14
        self.spine_threshold = 40
        self.brainstem_threshold = 40
        
    def get_thresholds(self, db):
        thresholds = (db.prescribed_doses.max() + 10)*np.ones((Constants.num_organs,))
        thresholds[self.eyes] = self.eye_threshold
        thresholds[self.spinal_cord] = self. spine_threshold
        thresholds[self.brainstem] = self.brainstem_threshold
        return thresholds

    def check_high_doses(self, db):
        outliers = set([])
        thresholds = self.get_thresholds(db)
        for organ_index in range(0, Constants.num_organs):
            organ_name = Constants.organ_list[organ_index]
            doses = db.doses[:, organ_index]
            bad_doses = np.where(doses > thresholds[organ_index])[0]
            for index in bad_doses:
                bad_dose = doses[index]
                id = db.ids[index]
                organ = OrganTuple(organ_name, int(bad_dose))
                outliers[id] = outliers.get(id, set([]))
                outliers[id].add(organ)
        return outliers

    def check_missing_organs(self, db):
        #checks for absence of no-eye organs and tumorless cases
        bad_patients = set([])
        no_volume = np.where(db.volumes <= 0.00001)
        for missing_organ in range(len(no_volume[0])):
            patient = no_volume[0][missing_organ]
            organ = no_volume[1][missing_organ]
            if organ not in self.eyes:
                bad_patients.add(patient)
                
        no_dose = np.where(db.doses <= 0.00001)
        for missing_organ in range(len(no_dose[0])):
            patient = no_dose[0][missing_organ]
            organ = no_dose[1][missing_organ]
            if organ not in self.eyes:
                bad_patients.add(patient)

        #check if position variance is too low (all in the center bascially)
        organ_spread = np.var(db.centroids, axis = 1).mean(axis = 1)
        bad_centroids = np.where( organ_spread < 200)[0]
        print('bad centroids', bad_centroids)
        for bad_patient in bad_centroids:
            bad_patients.add(bad_patient)
        #check if no tumor
        for patient in range(len(db.gtvs)):
            gtv_set = db.gtvs[patient]
            tumor_volume = np.sum([gtv.volume for gtv in gtv_set])
            if tumor_volume <= .00001:
                bad_patients.add(patient)
        
        bad_patients = bad_patients | self.get_data_outliers(db)
        return bad_patients
    
    def get_data_outliers(self, db, dose_match_threshold = .2, min_matches = 1):
        outliers = set([])
        n_patients = db.get_num_patients()
        for p1 in range(n_patients):
            x1 = db.doses[p1,:]
            num_matches = 0
            for p2 in range(p1+ 1, n_patients):
                x2 = db.doses[p2,:]
                dose_diff = np.sum(np.abs(x1 - x2))/np.sum(x1)
                if dose_diff < dose_match_threshold:
                    num_matches += 1
            if num_matches < min_matches:
                outliers.add(p1)
        print('outliers', outliers)
        return outliers
    
    def get_clean_subset(self, db):
        bad_patients = self.check_missing_organs(db)
        all_patients = set(range(db.get_num_patients()))
        good_patients = all_patients - bad_patients
        print(bad_patients)
        return sorted(good_patients)
        