#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 15:16:36 2019
@author: glavrent
This file test the python wrapping functions of the distance FORTRAN modules
"""

#load libraries 
import os
import shutil
import numpy as np

#change directories
os.chdir('/ssd/Research/Distance_metrics/testing/python')

#copy fortran compiled files
#fname_scr = '../../dist_metrics_progs/py_dist_metrics.so'
#fname_dst = 'py_dist_metrics.so'
#shutil.copyfile(fname_scr,fname_dst);
import py_dist_metrics as flt

#fault coordinates
flt_top_xyz  = np.array([[0.,  0.,  0.], 
                         [10., 10., 0],
                         [30., 10., 0.],
                         [30. ,40., 0.]])
flt_base_xyz = np.array([[-6.,  6.,-15.], 
                         [ 4., 16.,-15.],
                         [30., 16.,-15.],
                         [30. ,40.,-15.]])

#station coordinates
sta_xyz = np.array([0,10,0])
stas_xyz = np.array([[0,10,0],
                     [0,10,0],
                     [5,10,0]])

## Compute distance and projection points
# Compare with output of stand-alone function

#compute distances
dist_sta  = flt.py_wrapper_dist2sta(flt_top_xyz, flt_base_xyz, sta_xyz)
dist_stas = flt.py_wrapper_dist2stas(flt_top_xyz, flt_base_xyz, stas_xyz)

#compute distances and projection points
dist_prj_sta  = flt.py_wrapper_dist_prj2sta(flt_top_xyz, flt_base_xyz, sta_xyz) 
dist_prj_stas = flt.py_wrapper_dist_prj2stas(flt_top_xyz, flt_base_xyz, stas_xyz) 


