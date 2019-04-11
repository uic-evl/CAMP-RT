# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 15:18:10 2019

@author: Andrew
"""
import numpy as np
from Constants import Constants
from collections import OrderedDict, namedtuple

GTV = namedtuple('GTV', ['volume', 'position', 'doses', 'dists'])

class Patient():
    ##class holds information for each patient.
    ##is pased a series of dataframes (distacnes, doses, info) and extracts info
    ##p_id is a dummy id, position is the position of the patient in the whole dataset
    def __init__(self, distances, doses, p_id, group, info, use_distances = False):
        #patient ID number
        self.id = p_id
        self.group = group
        self.neck_boost = (info['Neck boost (Y/N)'] == 'Y')
        ##self.high_throat_dose = 1 if p_id in Constants.v2_high_throat_dose else 0
        #basically ordinality of the id, so where it will be in an index
        self.laterality = info['Tm Laterality (R/L)']
        self.prescribed_dose = info['Total dose']
        try:
            n_stage = info['N-category']
            self.n_stage = int(n_stage[1])
        except:
            print('error reading nstage for ', p_id)
            self.n_stage = -1
        self.tumor_subsite = info['Tumor subsite (BOT/Tonsil/Soft Palate/Pharyngeal wall/GPS/NOS)']
        centroid_data = self.get_doses_file_info(doses, distances)
        #I make a new matrix, so the order of centroid data isn't the same as the orginal csv
        self.doses = centroid_data[:, 4]
        self.min_doses = centroid_data[:, 5]
        self.max_doses = centroid_data[:, 6]
        #this is the actual total dose for all organs, 'total dose' from the data csv (66 or 70) is prescribed dose
        self.total_dose = np.sum(self.doses)
        self.volumes = centroid_data[:, 3]
        self.centroids = centroid_data[:, 0:3]
        #distances is a symetric matrix sorted by the Constants.organ_list
        if use_distances:
            self.distances = self.gen_distance_matrix(distances)
        #store the entries without gtvp for future study
        (self.tumor_volume, self.tumor_distances, self.tumor_position) = self.get_main_tumor()
        self.check_missing_organs(doses)
        ##self.check_if_full_dose()
        #report if there is no primary tumor
        if self.tumor_volume == 0 or np.sum(self.tumor_distances) == 0:
            Constants.no_tumor.append(self.id)

    def to_ordered_dict(self, dose_estimates):
        #exports local information into a dictionary
        entry = OrderedDict() #why is it ordered?
        entry['ID'] = str(self.id)
        entry['ID_int'] = int(self.id)
        entry['name'] = "Patient " + str(self.id)
        entry['tumorVolume'] = max([gtv.volume for gtv in self.gtvs])
        entry['organData'] = self.get_organ_data_dict(dose_estimates)
        entry['hasGTVp'] = str((self.gtvs[0].volume > 0)).lower()
        entry['hasGTVn'] = str((self.gtvs[1].volume > 0)).lower()
        #placeholders, will need to be populated by the Patientset class
        #there are non-ssim version in the original data but I think thats depricated (was for pearson?)
        entry['similarity_ssim'] = [0]
        entry['scores_ssim'] = [0]
        entry['laterality'] = self.laterality
        #skipping laterality int
        entry['tumorSubsite'] = self.tumor_subsite
        entry['total_Dose'] = self.prescribed_dose #this is confusing
        entry['cluster'] = self.group
        entry['gtvp_volume'] = self.gtvs[0].volume
        entry['gtvn_volume'] = self.gtvs[1].volume
        return(entry)

    def get_organ_data_dict(self, dose_estimates):
        #subset of the information for json export - is ordering important
        #should be a dictionary key = organ string, values = x,y,z,meanDose,maxDose
        data = OrderedDict()
        for x in range(0, Constants.num_organs):
            organ = Constants.organ_list[x]
            organ_dict = OrderedDict()
            organ_dict['x'] = self.centroids[x, 0]
            organ_dict['y'] = self.centroids[x, 1]
            organ_dict['z'] = self.centroids[x, 2]
            organ_dict['volume'] = self.volumes[x]
            organ_dict['meanDose'] = self.doses[x]
            organ_dict['minDose'] = self.min_doses[x]
            organ_dict['maxDose'] = self.max_doses[x]
            organ_dict['estimatedDose'] = 0 if np.isnan(dose_estimates[x]) else round(dose_estimates[x], 2)
            data[organ] = organ_dict
        if self.gtvs[0].volume > 0:
            gtvp = self.gtvs[0]
            gtvp_dict = OrderedDict()
            gtvp_dict['x'] = gtvp.position[0]
            gtvp_dict['y'] = gtvp.position[1]
            gtvp_dict['z'] = gtvp.position[2]
            gtvp_dict['volume'] = gtvp.volume
            gtvp_dict['meanDose'] = gtvp.doses[1]
            gtvp_dict['maxDose'] = gtvp.doses[2]
            gtvp_dict['minDose'] = gtvp.doses[0]
            data['GTVp'] = gtvp_dict
        if self.gtvs[1].volume > 0:
            gtvn_dict = OrderedDict()
            gtvn = self.gtvs[1]
            gtvn_dict['x'] = gtvn.position[0]
            gtvn_dict['y'] = gtvn.position[1]
            gtvn_dict['z'] = gtvn.position[2]
            gtvn_dict['volume'] = gtvn.volume
            gtvn_dict['meanDose'] = gtvn.doses[1]
            gtvn_dict['maxDose'] = gtvn.doses[2]
            gtvn_dict['minDose'] = gtvn.doses[0]
            data['GTVn'] = gtvn_dict
        return(data)

    def check_missing_organs(self, doses):
        #check if any organs are missing using the dose file, and store them
        organs = set(Constants.organ_list[:])
        dose_organs = set(doses['ROI'].unique())
        diff = organs - dose_organs
        if len(diff) > 0:
            Constants.missing_organs[self.id] = {'organs': diff}
        return

    def get_doses_file_info(self, doses, dists):
        #rename the columns so they're consistent, now done in dataset so I can read which set it's from
#        doses.columns = Constants.centroid_file_names
        #move centroids so the center of the cloud is at zero?
        centroids = self.center_centroids(doses)
        centroids = centroids.set_index('ROI')
        #extract a secondary tumor (only gets the first one?)
        #several patients have no gtvp but a gtvn

        def getGTV(name):
            try:
                gtv = centroids.loc[name]
                if not isinstance(gtv.volume, float):
                    Constants.multiple_gtv.append(self.id)
                    gtv = gtv.iloc[0]
                volume = gtv.volume
                position = gtv[['x','y','z']].values
                doses = gtv[['min_dose','mean_dose','max_dose']]
                distances = self.get_tumor_distances(name, dists)
            except:
                volume = float(0)
                position = np.array([0,0,0])
                doses = np.array([0,0,0])
                distances = np.zeros((Constants.num_organs,))
            return( GTV(volume, position, doses, distances) )
        self.gtvs = [getGTV('GTVp'), getGTV('GTVn')]
        gtvn_count = 2
        while(True):
            gtv_name = 'GTVn' + str(gtvn_count)
            newGTV = getGTV(gtv_name)
            if newGTV.volume > 0:
                self.gtvs.append(newGTV)
                gtvn_count += 1
            else:
                break
        print(len(self.gtvs))
        #get the info the centers, volumes, nad doses for all the things
        centroid_matrix = np.zeros((Constants.num_organs,7)) #row = x,y,z,volume,dose
        for idx in range(0, Constants.num_organs):
            organ = Constants.organ_list[idx]
            try:
                organ_entry = centroids.loc[organ]
                centroid_matrix[idx, 0:3] = organ_entry[['x','y','z']].values
                centroid_matrix[idx, 3] = organ_entry.volume
                centroid_matrix[idx, 4] = organ_entry.mean_dose
                centroid_matrix[idx, 5] = organ_entry.min_dose
                centroid_matrix[idx, 6] = organ_entry.max_dose
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
        alphabetical_organ_list = sorted(Constants.organ_list)
        for row in range(0, Constants.num_organs):
            for col in range(row + 1, Constants.num_organs):
                organ1 = alphabetical_organ_list[row]
                organ2 = alphabetical_organ_list[col]
                try:
                    dist_matrix[row, col] = (dists.loc[organ1, organ2])['Eucledian Distance (mm)']
                except:
                    try:
                        dist_matrix[row, col] = (dists.loc[organ2, organ1])['Eucledian Distance (mm)']
                    except:
                        dist_matrix[row, col] = 0
                        print(self.id, ' ', organ1, ' ', organ2, ' missing')
        dist_matrix += np.transpose(dist_matrix)
        return(dist_matrix)

    def get_tumor_distances(self, name, dists):
        #gets the tumor-organ distances
        gtv_dists = np.zeros((Constants.num_organs,))
        dists = dists.set_index(['Reference ROI', 'Target ROI']).sort_index()
        for idx in range(0, Constants.num_organs):
            organ = Constants.organ_list[idx]
            try:
                tumor_row = dists.loc[name, organ]
                gtv_dists[idx] = tumor_row['Eucledian Distance (mm)']
            except:
                try:
                    tumor_row = dists.loc[organ, name]
                    gtv_dists[idx] = tumor_row['Eucledian Distance (mm)']
                except:
                    gtv_dists[idx] = float(0)
        return(gtv_dists)

    def get_main_tumor(self):
        #basically gives a proxy so we use only the most important tumor?
        tumor_volume = 0.0
        tumor_distances = np.zeros((Constants.num_organs,))
        tumor_position = np.zeros((3,))
        try:
            for gtv in self.gtvs:
                tumor_volume += gtv.volume
                tumor_distances += gtv.volume*gtv.dists
                tumor_position += gtv.volume*gtv.position
            tumor_distances /= tumor_volume
            tumor_position /= tumor_position
        except:
            print('error reading tumor volume for ', self.id, self.gtvs)
        return(tumor_volume, tumor_distances, tumor_position)