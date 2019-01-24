#"""Python 2.7.15"""
# changed 3.5.2 to 2.7.15
# to run...
# python2.7 CAMPRTdata.py patients/


from collections import OrderedDict
from scipy.stats.stats import pearsonr
from operator import itemgetter
import ssim

import os
import glob
import csv
import json

import numpy as np
import re


np.set_printoptions(threshold=np.nan)
np.seterr(divide='ignore', invalid='ignore')

#myssim = pyssim.initialize()

# replace this with reading in CSV file of partitions
# modify other code relating to checking partition/master organs as well
masterList = [
'Brainstem',
'Cricoid_cartilage',
'Cricopharyngeal_Muscle',
'Esophagus',
'Extended_Oral_Cavity',
'Genioglossus_M',
#'Glottic_Area',
'Hard_Palate',
'Hyoid_bone',
'IPC',
'Larynx',
'Lower_Lip',
'Lt_Ant_Digastric_M',
'Lt_Anterior_Seg_Eyeball',
'Lt_Brachial_Plexus',
'Lt_Lateral_Pterygoid_M',
'Lt_Masseter_M',
'Lt_Mastoid',
'Lt_Medial_Pterygoid_M',
'Lt_Parotid_Gland',
'Lt_Posterior_Seg_Eyeball',
'Lt_Sternocleidomastoid_M',
'Lt_Submandibular_Gland',
'Lt_thyroid_lobe',
'Mandible',
'MPC',
'Mylogeniohyoid_M',
'Rt_Ant_Digastric_M',
'Rt_Anterior_Seg_Eyeball',
'Rt_Brachial_Plexus',
'Rt_Lateral_Pterygoid_M',
'Rt_Masseter_M',
'Rt_Mastoid',
'Rt_Medial_Pterygoid_M',
'Rt_Parotid_Gland',
'Rt_Posterior_Seg_Eyeball',
'Rt_Sternocleidomastoid_M',
'Rt_Submandibular_Gland',
'Rt_thyroid_lobe',
'Soft_Palate',
'SPC',
'Spinal_Cord',
'Supraglottic_Larynx',
'Thyroid_cartilage',
'Tongue',
'Upper_Lip',
'GTVp',
'GTVn'
]

def RunTestCases(organ):

    # GTVn2 is needed, corresponds to second tumor within patient

    # 'GTV node', 'GTV-N', 'GTV_n', 'GTVn', 'GTVn1', 'GTVn2'
    if organ in ('GTV node', 'GTV-N', 'GTV_n', 'GTVn', 'GTVn1'):
        organ = 'GTVn'
    elif organ in ('GTV primary', 'GTV-P', 'GTV_p', 'GTVp'):
        organ = 'GTVp'

    return organ


def FillOrganData(f, organRef, pID):
    od = OrderedDict()

    reader = csv.reader(f)
    headers = next(reader)

    for row in reader:
        organ = row[0]

        organ = RunTestCases(organ)

        if organ in masterList:

            if organ not in organRef:  # list keeps reference to main
                organRef += [organ]

            od[organ] = OrderedDict()  # organ name

            od[organ]['x'] = float(row[1])         # x pos
            od[organ]['y'] = float(row[2])         # y pos
            od[organ]['z'] = float(row[3])         # z pos

            # when data is missing, print it out
            # REPORT EVERYTHING TO TEXANS

            if row[5] != "":
                od[organ]['volume'] = float(row[5])  # volume
            else:
                od[organ]['volume'] = 0.0
                print (str(pID) + " missing volume")

            if row[4] != "":
                od[organ]['meanDose'] = float(row[4])  # mean dose
            else:
                od[organ]['meanDose'] = -1.0  # originally assigned (None), but pearson correlation returns Nan which is not JSON friendly
                print (str(pID) + " missing mean dose")

            if row[6] != "":
                od[organ]['minDose'] = float(row[6])  # min dose
            else:
                od[organ]['minDose'] = -1.0
                print (str(pID) + " missing min dose")

            if row[7] != "":
                od[organ]['maxDose'] = float(row[7])  # max dose
            else:
                od[organ]['maxDose'] = -1.0
                print (str(pID) + " missing max dose")

    return od


def FillMatrix(f, organRef, pID):
    # od = OrderedDict()
    reader = csv.reader(f)
    headers = next(reader)

    rows = list(reader)

    hasGTVp = False
    hasGTVn = False

    for row in rows:

        organ1 = RunTestCases(row[0])
        organ2 = RunTestCases(row[1])

        row[0] = organ1
        row[1] = organ2


        if organ1 in masterList and organ2 in masterList:

            if organ1 == 'GTVp' or organ2 == 'GTVp':
                hasGTVp = True

            if organ1 == 'GTVn' or organ2 == 'GTVn':
                hasGTVn = True

            if organ1 not in organRef:  # list keeps reference to main
                organRef += [organ1]

            if organ2 not in organRef:  # list keeps reference to main
                organRef += [organ2]

    array2D = np.ones((len(organRef), len(organRef)))

    array2D_tDist = np.zeros((len(organRef), len(organRef)))

    for row in rows:
        organ1 = row[0]
        organ2 = row[1]

        if organ1 in masterList and organ2 in masterList:

            if row[2] == "":
                print (str(pID) + " missing distance")

            array2D[organRef.index(organ1), organRef.index(organ2)] = row[2]

            # flip over diagonal
            array2D[organRef.index(organ2), organRef.index(organ1)] = row[2]

            if hasGTVp:
                if organ1 == 'GTVp':
                    array2D_tDist[organRef.index(organ2), organRef.index(organ2)] = row[2]
                elif organ2 == 'GTVp':
                    array2D_tDist[organRef.index(organ1), organRef.index(organ1)] = row[2]

    return [array2D, array2D_tDist, hasGTVp, hasGTVn]


def GetPatientByInternalID(internalID, patients):
    for p in patients:
        if internalID == p['ID_internal']:
            return p

    return None

def GetOrganDose(organ, organList):
    for o in organList.items():
        if o[0] == organ:
            return o[1]['meanDose']

    return None

class Patient_Set():

    def __init__(self):
        self.patients = []
        self.count = 1

    def get_CSVs(self, file_location):
            # go through each file path from directories
        pIDs = []
        file_list = glob.glob(file_location + '**/*.csv')
        file_list.sort(key = lambda file: max([int(x) for x in re.findall("[0-9]+", file)]))
        for fpath in file_list:
            with open(fpath, 'r') as f:  # open current file
                fname = os.path.basename(f.name)
                pID = fname[0:fname.index('_')]

                if pID not in pIDs:  # new dictionary entry
                    pIDs.append(pID)
                    pEntry = OrderedDict()  # create new patient entry
                    pEntry['ID'] = pID  # converting to float breaks something?
                    pEntry['ID_int'] = int(pID)
                    pEntry['name'] = "Patient " + str(pID)
                    pEntry['tumorVolume'] = 0.0

                    if '_cent' in fname:  # parse centroids file
                        pEntry['organData'] = FillOrganData(f, self.organRef, pID)
                    else:
                        pEntry['organData'] = {}  # placeholder

                    pEntry['ID_internal'] = self.count

                    if '_dist' in fname:  # parse distances file
                        data = FillMatrix(f, self.organRef, pID)
                        pEntry['matrix'] = data[0]
                        pEntry['matrix_tumorDistances'] = data[1]
                        pEntry['hasGTVp'] = data[2]
                        pEntry['hasGTVn'] = data[3]
                    else:
                        pEntry['matrix'] = []  # placeholder
                        pEntry['matrix_tumorDistances'] = []
                        pEntry['hasGTVp'] = False
                        pEntry['hasGTVn'] = False

                    pEntry['matrix_ssim'] = []
                    pEntry['matrix_ssim_dist'] = []
                    pEntry['matrix_ssim_vol'] = []
                    pEntry['matrix_dose'] = []
                    pEntry['matrix_TumorVolume'] = []
                    pEntry['matrix_pos'] = []
                    pEntry['similarity'] = []
                    pEntry['scores'] = []
                    pEntry['similarity_ssim'] = []
                    pEntry['scores_ssim'] = []
                    pEntry['laterality'] = ""
                    pEntry['laterality_int'] = -1
                    pEntry['tumorSubsite'] = ""

                    self.patients.append(pEntry)

                    self.count += 1
                else:                # modify existing dictionary entry
                    pEntry = next(
                        (item for item in self.patients if item['ID'] == pID), None)  # find entry

                    if '_cent' in fname:  # parse centroids file
                        if pEntry != None:
                            pEntry['organData'] = FillOrganData(f, self.organRef, pID)

                    if '_dist' in fname:  # parse distances file
                        if pEntry != None:
                            data = FillMatrix(f, self.organRef, pID)
                            pEntry['matrix'] = data[0]
                            pEntry['matrix_tumorDistances'] = data[1]
                            pEntry['hasGTVp'] = data[2]
                            pEntry['hasGTVn'] = data[3]
        return

    def read_Laterality(self, file_string = "laterality.csv"):
        # alg compares tumor distances, then for laterality left and right are equal? is this right?
        # read in laterality
        with open(self.data_folder + file_string, 'r') as csvFile:
            reader = csv.reader(csvFile)
            header = next(reader)

            for row in reader:
                for p in self.patients:
                    if str(row[0]) == p['ID']:
                        p['laterality'] = str(row[1])
                        p['tumorSubsite'] = str(row[2])

                        if str(row[1]) in ["L", "R"]:
                            p['laterality_int'] = 0
                        elif str(row[1]) == "Bilateral":
                            p['laterality_int'] = 1
        return

    def read_total_dosage(self, file_string = "patient_info.csv'"):
        # read in total dose for each patient
        with open(self.data_folder + file_string, 'r') as csvFile:
            reader = csv.reader(csvFile)
            header = next(reader)

            for row in reader:
                for p in self.patients:
                    if str(row[0]) == p['ID']:
                        p['total_Dose'] = float(row[31])
        return

    def delete_tumors(self):
        # now deleting both tumors, use has_key or indexOf to make sure the right rows/columns are being deleted
        for currP in self.patients:
            if self.organRef[0] != "GTVn" and self.organRef[1] != "GTVp":
                print ("WARNING: GTV in wrong row/column. Incorrect format.")

            currP['matrix'] = np.delete(currP['matrix'], (0), axis=0) # delete first row
            currP['matrix'] = np.delete(currP['matrix'], (0), axis=1) # delete first column

            currP['matrix_tumorDistances'] = np.delete(currP['matrix_tumorDistances'], (0), axis=0) # delete first row
            currP['matrix_tumorDistances'] = np.delete(currP['matrix_tumorDistances'], (0), axis=1) # delete first column


            # DO IT AGAIN FOR GTVp

            currP['matrix'] = np.delete(currP['matrix'], (0), axis=0) # delete first row
            currP['matrix'] = np.delete(currP['matrix'], (0), axis=1) # delete first column

            currP['matrix_tumorDistances'] = np.delete(currP['matrix_tumorDistances'], (0), axis=0) # delete first row
            currP['matrix_tumorDistances'] = np.delete(currP['matrix_tumorDistances'], (0), axis=1) # delete first column
        return

    def populate_something(self):
        #TODO: Figure out what this block of code is doing
        for p in self.patients:
            #p['matrix'].resize((len(organRef), len(organRef)))
            # add padding to matrix so all matrices are the same size
            padSize = len(self.organRef) - p['matrix'].shape[0]
            p['matrix'] = np.lib.pad(p['matrix'], ((0, padSize), (0, padSize)), mode='constant')
            p['matrix_tumorDistances'] = np.lib.pad(p['matrix_tumorDistances'], ((0, padSize), (0, padSize)), mode='constant')

            organs = p['organData']

            matrixCopy = np.copy(p['matrix'])

            # initiliaze dose matrix
            p['matrix_dose'] = np.zeros((len(self.organRef), len(self.organRef)))

            # initiliaze tumor volume matrix
            p['matrix_TumorVolume'] = np.zeros((len(self.organRef), len(self.organRef)))

            if p['hasGTVp']:
                p['tumorVolume'] = p['organData']['GTVp']['volume']

            for organ in organs.items():  # populate diagonal of matrix with mean dose data

                if organ[0] != "GTVp" and organ[0] != "GTVn":
                    #populate dose matrix
                    # using dose matrix for total dose now
                    p['matrix_dose'][self.organRef.index(organ[0]), self.organRef.index(organ[0])] = p['total_Dose']
                    p['matrix_TumorVolume'][self.organRef.index(organ[0]), self.organRef.index(organ[0])] = p['tumorVolume']

                    # populate position matrix
                    ##posMatrix[0, organRef.index(organ[0])] = organ[1]['x']
                    ##posMatrix[1, organRef.index(organ[0])] = organ[1]['y']
                    ##posMatrix[2, organRef.index(organ[0])] = organ[1]['z']

            p['matrix_ssim'] = np.dot(matrixCopy, p['matrix_dose'])

            # SHOULD WE DELETE ROWS/COLS first before dot product?

            p['matrix_ssim_dist'] = np.dot(matrixCopy, p['matrix_tumorDistances'])
            p['matrix_ssim_vol'] = np.dot(matrixCopy, p['matrix_TumorVolume'])
        return

    def get_SSIM_score(self):
        # calculate ssim score
        for currP in self.patients:
            correlations = []
            ssimResults = []
            for nextP in self.patients:
                pCoeff_1 = pearsonr(currP['matrix'].flat, nextP['matrix'].flat)[0]
                ssimScor_totDose = ssim.compute_ssim(currP['matrix_ssim'], nextP['matrix_ssim'])
                ssimScor_dist = ssim.compute_ssim(currP['matrix_ssim_dist'], nextP['matrix_ssim_dist'])
                ssimScor_vol = ssim.compute_ssim(currP['matrix_ssim_vol'], nextP['matrix_ssim_vol'])

                if currP['laterality_int'] == nextP['laterality_int']:
                    ssimScor = (ssimScor_totDose + ssimScor_dist + ssimScor_vol + 1) / 4.0
                else:
                    ssimScor = (ssimScor_totDose + ssimScor_dist + ssimScor_vol + 0) / 4.0


                #ssimScor = 1.0
                if (ssimScor <= 0.0):
                    #ssimScor = 0
                    print(ssimScor)

                self.pSimMatrix[currP['ID_internal'], nextP['ID_internal']] = ssimScor
                self.pSimMatrix[0, currP['ID_internal']] = int(currP['ID_internal'])
                self.pSimMatrix[currP['ID_internal'], 0] = int(currP['ID_internal'])

                correlations.append((nextP['ID_internal'], pCoeff_1))
                ssimResults.append((nextP['ID_internal'], ssimScor))

            correlations = sorted(correlations, key=itemgetter(1), reverse=True)
            ssimResults = sorted(ssimResults, key=itemgetter(1), reverse=True)

            for score in correlations:
                currP["similarity"].append(score[0])
                currP["scores"].append(score[1])

            for score in ssimResults:
                currP["similarity_ssim"].append(score[0])
                currP["scores_ssim"].append(score[1])
        return

    def generate_differences_csv(self, file_name = 'differences.csv'):
        # Generate spreadsheet with predictions
        # missing organ data is replaced with -1
        spreadsheet = []
        header = ["ID", "N1_ID", "N2_ID", "N3_ID", "N4_ID", "N5_ID"]

        for organ in self.organRef:
            header.append(organ)

        header.append("Sum")
        header.append("Average")

        spreadsheet.append(header)

        for p in self.patients:
            row = []
            row.append(p["ID"])

            ranks = p["similarity_ssim"]
            scores = p["scores_ssim"]

            difference = []

            # neighbors
            for i in range(1, self.num_comparisons+1):
                neighbor = GetPatientByInternalID(ranks[i], self.patients)

                if neighbor != None:
                    row.append(neighbor["ID"])
                else:
                    row.append(-1)

            # organ differences
            for organ in self.organRef:

                organAverage = 0

                for i in range(1, self.num_comparisons + 1):

                    neighbor = GetPatientByInternalID(ranks[i], self.patients)

                    if neighbor != None:
                        organDose = GetOrganDose(organ, neighbor["organData"])

                        if organDose != None:
                            #organAverage += (organDose * scores[i])
                            organAverage += (organDose)

                actualDose = GetOrganDose(organ, p["organData"])

                if actualDose != None:
                    difference.append(round(abs(actualDose - (organAverage / 5.0)), 3))
                else:
                    difference.append(-1)

            for value in difference:
                row.append(value)

            row.append(np.sum(difference))
            row.append(round(np.average(difference), 3))
            spreadsheet.append(row)

        with open(self.write_folder + file_name, 'w+') as csvfile:
            csvWriter = csv.writer(csvfile, delimiter=',')
            csvWriter.writerows(spreadsheet)
        return

    def write_data(self):
        self.generate_differences_csv()

        #TODO figure out what these are
#        with open(self.write_folder + "matrix_p222_ssim_noDoses.csv", "w+") as my_csv:
#            csvWriter = csv.writer(my_csv, delimiter=',')
#            csvWriter.writerows(self.matrices[1]['matrix_ssim'])

        with open(self.write_folder + "pSimMatrix_noDoses.csv", "w+") as my_csv:
            csvWriter = csv.writer(my_csv, delimiter=',')
            csvWriter.writerows(self.pSimMatrix)

        with open(self.write_folder + 'patients_SSIM_wDoses_wDists.json', 'w+') as f:  # generate JSON
            json.dump( self.delete_matrices(self.patients), f, indent=4)

    def delete_matrices(self, patients):
        #saves the matrices in a seperate variable and then deletes then so they are saved separately
        #converts the numpy arrays to lists so that they are saveable
        for p in self.patients:  # json can't handle too many matrices, delete the matrices
            del p['matrix']
            del p['matrix_pos']
            del p['matrix_ssim']
            del p['matrix_ssim_dist']
            del p['matrix_ssim_vol']
            del p['matrix_dose']
            del p['matrix_tumorDistances']
            del p['matrix_TumorVolume']
        return patients

    def save_matrices(self):
        #saves data relevant for the ssim inputs to matrices.json
        ssim_inputs = {
            p['ID_internal']: {
                    'laterality': p['laterality'], #right, left, or center
                    'matrix_ssim': p['matrix_ssim'].tolist(), #dose for ssim?
                    'matrix_ssim_dist': p['matrix_ssim_dist'].tolist(), #distances for ssim
                    'matrix_ssim_vol': p['matrix_ssim_vol'].tolist() #volumes for ssim
                }
                for p in self.patients
            }

        #saves the original data before dot-products
        feature_arrays = {
            p['ID_internal']: {
                    'laterality': p['laterality'], #right, left, or center
                    'organ_distances': p['matrix'].tolist(), #matrix of organ-organ distances
                    'tumor_volumes': p['matrix_TumorVolume'][0,0], #GTVp Volume
                    'tumor_distances': np.diag(p['matrix_tumorDistances']).tolist(), #tumor to organ distance
                    'total_doses': p['matrix_dose'][0,0] #total dose
                }
                for p in self.patients
            }

        with open(self.write_folder + "matrices.json", 'w+') as f:
            json.dump(ssim_inputs, f, indent = 4)
        with open(self.write_folder + 'features.json', 'w+') as f:
            json.dump(feature_arrays, f, indent = 4)
        return


    def run(self, argv = "patients_v2\\", num_comparisons = 5, write = False, write_folder = "latest_results\\", data_folder = "data\\"):
        self.num_comparisons = num_comparisons
        self.write_folder = write_folder
        self.data_folder = data_folder
        self.organRef = []
        self.organRef += ['GTVn']
        self.organRef += ['GTVp']

        self.get_CSVs(argv)
        self.read_Laterality()
        self.read_total_dosage()
        self.delete_tumors()

        # make sure matrices are all same size [check]
        # fill matrix diagonal [check]
        # ask liz how to handle organ values missing for patients
        # what do do with negative distances
        # delete columns and rows corresponding to GTVn # ---------------------------------------------
        # ended up not using GTVn for similarity computation

        self.organRef.remove("GTVp")
        self.organRef.remove("GTVn")

        self.populate_something();

        self.pSimMatrix = np.zeros((len(self.patients) + 1, len(self.patients) + 1))
        self.get_SSIM_score()
        self.patients.sort(key = lambda x: x['ID_int'])
        if write: #writes the results to csv files, default false
            self.write_data()
        return

    def get_patients(self):
        return self.patients

    def get_matrices(self):
        return self.matrices

# __name__ == '__main__':
    # command-line argument specifies
    # patients parent directory
 #   main()#sys.argv[1])
data_set = Patient_Set()
data_set.run(write=False)
data_set.save_matrices()
#np.savetxt("latest_results\\all_ssim_scores.csv", data_set.pSimMatrix, delimiter=",")
