# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 17:22:31 2019

@author: Andrew
"""
from glob import glob
from re import findall
import numpy as np
import pandas as pd

#sorts by size of largest integer string, which is the id for our files
file_sort = lambda x: sorted(x, key =
                             lambda file: max([int(x) for x in findall("[0-9]+", file)])
                        )
distance_files = file_sort(glob('patients_v2\\' + '**/*distances.csv'))
dose_files = file_sort(glob('patients_v2\\' + '**/*meandoses.csv'))
#this would be faster if I zipped them maybe?
ids = [max([int(x) for x in findall('[0-9]+', file)]) for file in distance_files]

patient_distances = pd.read_csv(distance_files[1])
patient_dose = pd.read_csv(dose_files[1])