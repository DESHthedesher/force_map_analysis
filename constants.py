# -*- coding: utf-8 -*-
"""
Created on Thu May 22 15:16:37 2014

@author: Duncan
"""

#relating to path dependencies
PATH_TO_DATA = "/data/"
PATH_TO_OUTFILE = "/data/results"

#relating to hard data
SPRING_CONSTANT = 21.18*10**-3   #N/m
INVOLS = 138.63 #Defl InVols, nm/V

#relating to data analysis parameters
LEFT_TAIL_LIMIT = 0.75 #point the closest to the surface you trust the tip not to interact with the surface
RIGHT_TAIL_LIMIT = 0.92 #point the farthest from the surface you trust not to be funny for all force curves
SLIDING_WINDOW_SIZE = 15 #data points
DETECTION_LIMIT = 5 #number of STDEVS that constitite a "point which is probably off"
MIN_FIT_LENGTH = 3*SLIDING_WINDOW_SIZE #number of consecutive points for a fit region to qualify. SHould be at least 3 windows.
