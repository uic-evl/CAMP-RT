# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 09:19:09 2019

@author: Andrew
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json

patients_file = "data\\patients_SSIM_wDoses_wDists.json"
ssim_scores_file = "latest_results\\all_ssim_scores.csv"

with open(patients_file) as file:
    data = json.load(file)

ssim_scores = pd.read_csv(ssim_scores_file)
s