# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 12:07:35 2019

@author: Andrew
"""
import numpy as np
from PatientSet import PatientSet
from Constants import Constants
from collections import namedtuple

OrganTuple = namedtuple('OrganTuple', ['organ', 'dose'])
class ErrorChecker():

    def __init__(self):
        self.ra_eye = Constants.organ_list.index('Rt_Anterior_Seg_Eyeball')
        self.rp_eye = Constants.organ_list.index('Rt_Posterior_Seg_Eyeball')
        self.la_eye = Constants.organ_list.index('Lt_Anterior_Seg_Eyeball')
        self.lp_eye = Constants.organ_list.index('Lt_Posterior_Seg_Eyeball')
        self.brainstem = Constants.organ_list.index('Brainstem')
        self.spinal_cord = Constants.organ_list.index('Spinal_Cord')
        self.eye_threshold = 14
        self.spine_threshold = 40
        self.brainstem_threshold = 40
        
    def get_thresholds(self, db):
        thresholds = db.prescribed_doses + 10
        thresholds[[self.ra_eye, self.rp_eye, self.la_eye, self.lp_eye]] = self.eye_threshold
        thresholds[self.spinal_cord] = self. spine_threshold
        thresholds[self.brainstem] = self.brainstem_threshold
        return thresholds

    def check_high_doses(self, db):
        outliers = {}
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
