# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 15:18:10 2019

@author: Andrew
"""
import numpy as np
from Constants import Constants

class Patient():
    ##class holds information for each patient.  
    ##is pased a series of dataframes (distacnes, doses, info) and extracts info
    ##p_id is a dummy id, position is the position of the patient in the whole dataset
    def __init__(self, distances, doses, p_id, position, info):
        #patient ID number
        self.id = p_id
        #basically ordinality of the id, so where it will be in an index
        self.pos = position
        self.laterality = info['Tm Laterality (R/L)']
        self.prescribed_dose = info['Total dose']
        centroid_data = self.get_doses_file_info(doses)
        self.doses = centroid_data[:, 4]
        #####normalize to total dose and then dose proportions
        self.total_dose = np.sum(self.doses)
        ######################
        self.volumes = centroid_data[:, 3]
        self.centroids = centroid_data[:, 0:3]
        self.distances = self.gen_distance_matrix(distances)
        (self.gtvp_dists, self.gtvn_dists) = self.get_tumor_distances(distances)
        #store the entries without gtvp for future study
        (self.tumor_volume, self.tumor_distances, self.tumor_position) = self.get_main_tumor()
        self.check_missing_organs(distances, doses)
        #report if there is no primary tumor
        if self.tumor_volume == 0 or np.sum(self.tumor_distances) == 0:
            Constants.no_tumor.append(self.id)

    def check_missing_organs(self, distances, doses):
        #check if any organs are missing using the dose file, and store them
        organs = set(Constants.organ_list[:])
        dose_organs = set(doses['ROI'].unique())
        diff = organs - dose_organs
        #for missing_organ in diff:
            #print('patient ', self.id, ' at index ', self.pos, ' is missing organ ', missing_organ)
        if len(diff) > 0:
            Constants.missing_organs[self.pos] = {'id': self.id, 'organs': diff}
        return

    def get_doses_file_info(self, doses):
        #rename the columns so they're consistent
        doses.columns = Constants.centroid_file_names
        #move centroids so the center of the cloud is at zero?
        centroids = self.center_centroids(doses)
        centroids = centroids.set_index('ROI')
        #extract the primary tumor info.
        try:
            gtvp = centroids.loc['GTVp']
            self.gtvp_volume = gtvp.volume
        except:
            self.gtvp_volume = float(0)
        try:
            self.gtvp_position = gtvp[['x','y','z']].values
        except:
            self.gtvp_position = np.array([0,0,0])
        #extract a secondary tumor (only gets the first one?)
        #several patients have no gtvp but a gtvn
        try:
            gtvn = centroids.loc['GTVn']
            if not isinstance(gtvn.volume, float):
                Constants.multiple_gtvn.append(self.id)
                gtvn = gtvn.iloc[0]
            self.gtvn_volume = gtvn.volume
        except:
            self.gtvn_volume = float(0)
        try:
            self.gtvn_position = gtvn[['x','y','z']].values
        except:
            self.gtvn_position = np.array([0,0,0])
        #get the info the centers, volumes, nad doses for all the things
        centroid_matrix = np.zeros((Constants.num_organs,5)) #row = x,y,z,volume,dose
        for idx in range(0, Constants.num_organs):
            organ = Constants.organ_list[idx]
            try:
                organ_entry = centroids.loc[organ]
                centroid_matrix[idx, 0:3] = organ_entry[['x','y','z']].values
                centroid_matrix[idx, 3] = organ_entry.volume
                centroid_matrix[idx, 4] = organ_entry.mean_dose
            except:
                pass
                #print('patient ', self.id, ' is missing organ ', organ, ' centroid data')
        return(centroid_matrix)

    def center_centroids(self, centroids):
        #subtract off the mean so the pointcloud is centered at 0
        #should I just use a reference organ instead?  or rotate?
        centroids.x -= centroids.x.mean()
        centroids.y -= centroids.y.mean()
        centroids.z -= centroids.z.mean()
        return(centroids)

    def gen_distance_matrix(self, dists):
        #generates a symetric 45x45 matrix of organ-organ distances
        dist_matrix = np.zeros(( Constants.num_organs, Constants.num_organs))
        dists = dists.set_index(['Reference ROI', 'Target ROI']).sort_index()
        for row in range(0, Constants.num_organs):
            for col in range(row + 1, Constants.num_organs):
                organ1 = Constants.organ_list[row]
                organ2 = Constants.organ_list[col]
                try:
                    dist_matrix[row, col] = (dists.loc[organ1, organ2])['Eucledian Distance (mm)']
                except:
                    dist_matrix[row, col] = 0
        dist_matrix += np.transpose(dist_matrix)
        return(dist_matrix)

    def get_tumor_distances(self, dists):
        #gets the tumor-organ distances
        gtvp_dists = np.zeros((Constants.num_organs,))
        gtvn_dists = np.zeros((Constants.num_organs,))
        dists = dists.set_index(['Reference ROI', 'Target ROI']).sort_index()
        for idx in range(0, Constants.num_organs):
            organ = Constants.organ_list[idx]
            try:
                tumor_row = dists.loc['GTVp', organ]
                gtvp_dists[idx] = tumor_row['Eucledian Distance (mm)']
            except:
                try:
                    tumor_row = dists.loc[organ, 'GTVp']
                    gtvp_dists[idx] = tumor_row['Eucledian Distance (mm)']
                except:
                    gtvp_dists[idx] = float(0)
            try:
                tumor_row = dists.loc['GTVn', organ]
                gtvn_dists[idx] = tumor_row['Eucledian Distance (mm)']
            except:
                try:
                    tumor_row = dists.loc[Constants.organ_list[idx], 'GTVn']
                    gtvn_dists[idx] = tumor_row['Eucledian Distance (mm)']
                except:
                    gtvn_dists[idx] = float(0)
        return((gtvp_dists, gtvn_dists))

    def get_main_tumor(self):
        #basically gives a proxy so we use only the most important tumor?
        tumor_volume = self.gtvn_volume + self.gtvp_volume
        tumor_distances = (self.gtvn_dists*self.gtvn_volume + self.gtvp_dists*self.gtvp_volume)/(
                self.gtvn_volume + self.gtvp_volume)
        tumor_position = (self.gtvn_position*self.gtvn_volume + self.gtvp_position*self.gtvp_volume)/(
                self.gtvn_volume + self.gtvp_volume)
        if (self.gtvn_volume != 0 and np.sum(self.gtvn_dists) == 0) or (
                self.gtvn_volume == 0 and np.sum(self.gtvn_dists) != 0) or (
                        self.gtvp_volume != 0 and np.sum(self.gtvp_dists) == 0) or (
                            self.gtvp_volume == 0 and np.sum(self.gtvp_dists) != 0):
            print('patient ', self.id, 'is having some issues with tumor consistency')
        return(tumor_volume, tumor_distances, tumor_position)