# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 11:37:00 2019

@author: Andrew
"""
import PatientSet as ps
from Constants import Constants
import sys

def main(root):
    dataset = ps.PatientSet(root = root, outliers = Constants.v2_bad_entries)
    dataset.export()

if __name__ == '__main__':
    main(sys.argv[1])
    