# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 10:12:11 2019

@author: Andrew
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from Constants import Constants
from patient_processing import Rankings, PatientSet, Patient

db = pickle.load(open('data\\patient_data_v2_only.p', 'rb'))
total_doses = db.total_doses