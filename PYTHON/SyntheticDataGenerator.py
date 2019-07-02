# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 15:39:18 2019
Big function that makes a functional (as of writing) synthetic patientset
doesn't capture all metadata, but does data needed for tsim and the json file
@author: Andrew
"""

from sklearn.ensemble import RandomForestRegressor
from collections import namedtuple
from Metrics import pca
import numpy as np
import pandas as pd
from PatientSet import PatientSet
from Constants import Constants
from copy import copy

GTV = namedtuple('GTV', ['name', 'volume', 'position', 'doses', 'dists', 'organ'])
class ClassStats():
    
    def __init__(self, db, c, extra_tumors = 2):
        self.cluster = c
        self.extra_tumors = extra_tumors
        class_args = np.argwhere(db.classes == c).ravel()
        self.density = len(class_args)/len(db.classes)
        tumor_centroids = []
        tumor_volumes = []
        tumor_distances = []
        training_organ_centroids = []
        training_input_data = []
        self.n_tumors = []
        for patient in class_args:
            tumorset = db.gtvs[patient]
            n_tumors = 0
            for tumor in tumorset:
                if tumor.volume > 0:
                    tumor_centroids.append(tumor.position)
                    tumor_volumes.append(tumor.volume)
                    tumor_distances.append(tumor.dists)
                    training_organ_centroids.append(db.centroids[patient].ravel())
                    training_input_data.append(db.prescribed_doses[patient])
                    n_tumors += 1
            self.n_tumors.append(n_tumors)
        tumor_centroids = np.vstack(tumor_centroids).astype('float32')
        tumor_volumes = np.vstack(tumor_volumes).astype('float32')
        training_input_data = np.vstack(training_input_data).astype('float32')
        training_organ_centroids = np.vstack(training_organ_centroids).astype('float32')
        tumor_distances  = np.vstack(tumor_distances).astype('float32')
        
        self.tumor_volumes = (tumor_volumes.mean(), tumor_volumes.std())
        self.tumor_volume_bounds = (tumor_volumes.min(), tumor_volumes.max())
        self.tumor_centroids = (tumor_centroids.mean(axis = 0), tumor_centroids.std(axis = 0))
        self.prescribed_doses = (db.prescribed_doses[class_args].mean(), 
                                 db.prescribed_doses[class_args].std())
        self.organ_centroids = [db.centroids[arg] for arg in class_args]
        self.organ_volumes = (db.volumes[class_args].mean(axis = 0),
                              db.volumes[class_args].mean(axis = 0))
        self.organ_volume_bounds = (db.volumes[class_args].min(axis = 0),
                              db.volumes[class_args].max(axis = 0))
        self.tumor_subsites = db.subsites[class_args]
        
        training_input = np.hstack([tumor_centroids, tumor_volumes, 
                                    training_input_data,
                                    training_organ_centroids])
    
        self.distance_generator = RandomForestRegressor(n_estimators = 20)
        self.distance_generator.fit(training_input, tumor_distances)
        self.dose_generator = RandomForestRegressor(n_estimators = 20)
        
        dose_input = np.hstack([db.tumor_distances[class_args],
                                           db.prescribed_doses[class_args].reshape(-1,1)])
        class_doses = db.doses[class_args]
        noised_doses = class_doses + np.random.normal(loc = 0, 
                                                      scale = 2*class_doses.std(axis=0), 
                                                      size = class_doses.shape)
        noised_doses = np.maximum(class_doses/2, noised_doses)
        self.dose_generator.fit(dose_input, noised_doses)
        
    def generate_one(self):
        #generates a dictionary with a bunch of values for a synthetic patient   
        sample = {}
        num_tumors = np.random.choice(self.n_tumors) +np.random.randint(0, self.extra_tumors)
        dose = int(np.random.normal(self.prescribed_doses[0], 
                                    self.prescribed_doses[1]))
        t_volumes = np.random.normal(self.tumor_volumes[0], 
                                     self.tumor_volumes[1], 
                                     num_tumors)
        t_volumes = np.minimum(t_volumes, self.tumor_volume_bounds[1])
        t_volumes = np.maximum(t_volumes, self.tumor_volume_bounds[0])
        t_volumes = sorted(t_volumes, key = lambda x: -x)
        o_centroids = self.organ_centroids[np.random.randint(0,len(self.organ_centroids))] * np.random.normal(1,.05, 3)
        gtvs= []
        min_dists = np.inf*np.ones((1,Constants.num_organs))
        sides = set([])
        for t in range(max([num_tumors,2])):
            if t == 0:
                name = 'GTVp'
            elif t == 1:
                name = 'GTVn'
            else:
                name = 'GTVn' + str(t)
                
            if t >= (num_tumors):
                vol = 0
                dists = np.zeros(min_dists.shape)
                center = np.zeros((3,))
                dose_series = pd.Series([0,0,0], 
                                        ['min_dose', 'mean_dose', 'max_dose'])
                organ = 'NA'
            else:
                vol = t_volumes[t]/(t+1)
                center = np.random.normal(self.tumor_centroids[0], 
                           self.tumor_centroids[1])
                center = np.maximum(o_centroids.min(axis = 0), center)
                center = np.minimum(o_centroids.max(axis = 0), center)
                if center[0] > 0:
                    sides.add('L')
                else:
                    sides.add('R')
                x = np.hstack([center, 
                               vol, 
                               dose,
                               o_centroids.ravel()])
                dists = self.distance_generator.predict(x.reshape(1,-1))
                min_dists = np.minimum(min_dists, dists)
                temp_doses = np.array([dose - 3*np.random.rand(), dose, dose + 3*np.random.rand()])
                dose_series = pd.Series(temp_doses, ['min_dose', 'mean_dose', 'max_dose'])
                organ = Constants.organ_list[np.argmin(dists)]
            new_gtv = GTV(name = name, position = center, doses = dose_series,
                          organ = organ, dists = dists, volume = vol)
            gtvs.append(new_gtv)
        sample['min_distances'] = min_dists
        sample['gtvs'] = gtvs 
        sample['prescribed_dose'] = dose
        sample['organ_centroids'] = o_centroids
        sample_volumes = np.random.normal(self.organ_volumes[0],
              np.sqrt(self.organ_volumes[1]))
        sample['organ_volumes'] = np.maximum(sample_volumes, self.organ_volume_bounds[0])
        dose_input = np.append(np.copy(min_dists), dose)
        sample['mean_doses'] = self.dose_generator.predict(dose_input.reshape(1,-1))
        sample['subsite'] = np.random.choice(self.tumor_subsites)
        if 'L' in sides:
            if 'R' in sides:
                sample['laterality'] = 'B'
            else:
                sample['laterality'] = 'L'
        else:
            sample['laterality'] = 'R'
        sample['cluster'] = self.cluster
        return(sample)
        
        
def generate_synthetic_dataset(db, patients_to_generate = 200, extra_tumors = 1):
    #takes a patientset and returns a copy with key values changed to generated values
    #some currently unused data won't be changed tho.  
    #patients to generate is number of patients
    #extra tumors add between 0 and x extra tumors to each patient when generating data, for testing
    #out matching for extra tumor patients
    class_stats = []
    for c in np.unique(db.classes):
        class_stats.append(ClassStats(db, c, extra_tumors))

    densities = [c.density for c in class_stats]
    patients =[]
    generated_tumor_distances = np.zeros((patients_to_generate, Constants.num_organs))
    generated_organ_volumes = np.zeros((patients_to_generate, Constants.num_organs))
    generated_doses = np.zeros((patients_to_generate, Constants.num_organs))
    generated_clusters = np.zeros((patients_to_generate,))
    generated_organ_centroids = np.zeros((patients_to_generate, Constants.num_organs, 3))
    generated_ids = []
    generated_lateralities = []
    generated_subsites = []
    generated_prescribed_doses = np.zeros((patients_to_generate,))
    generated_gtvs = []
    curr_id = 0
    for p in range(patients_to_generate):
        class_generator = np.random.choice(class_stats, p = densities)
        fake_patient = class_generator.generate_one()
        patients.append(fake_patient)
        generated_tumor_distances[p,:] = fake_patient['min_distances']
        generated_organ_volumes[p,:] = fake_patient['organ_volumes']
        generated_doses[p,:] = fake_patient['mean_doses']
        generated_clusters[p] = fake_patient['cluster']
        generated_organ_centroids[p] = fake_patient['organ_centroids']
        curr_id += np.random.randint(1,4)
        generated_ids.append(curr_id)
        generated_gtvs.append(fake_patient['gtvs'])
        generated_lateralities.append(fake_patient['laterality'])
        generated_subsites.append(fake_patient['subsite'])
        generated_prescribed_doses[p] = fake_patient['prescribed_dose']
    fake_db = copy(db)
    fake_db.tumor_distances = generated_tumor_distances
    fake_db.doses = generated_doses
    fake_db.classes = generated_clusters
    fake_db.volumes = generated_organ_volumes
    fake_db.centroids = generated_organ_centroids
    fake_db.max_doses = generated_doses + 3*np.random.random()
    fake_db.min_doses = generated_doses - 3*np.random.random()
    fake_db.ids = generated_ids
    fake_db.lateralities = generated_lateralities
    fake_db.subsites = generated_subsites
    fake_db.prescribed_doses = generated_prescribed_doses
    fake_db.gtvs = generated_gtvs
    return fake_db