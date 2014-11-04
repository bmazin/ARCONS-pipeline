#!/bin/python
'''
rebin.py

Created by Seth Meeker on 1-29-2013

Given arrays of wavelengths and fluxes (x and y) rebins to specified bin size by taking average value of input data within each bin

use: rebinnedData = rebin(x,y,binedges)

binedges typically can be imported from a FlatCal after being applied in an ObsFile

'''

import tables
import numpy as np
import matplotlib.pyplot as plt

def rebin(x,y,binedges):
    #must be passed binedges array since spectra will not be binned with evenly sized bins
    start = binedges[0]
    stop = x[-1]
    #calculate how many new bins we will have
    nbins = len(binedges)-1
    #create output arrays
    rebinned = np.zeros((nbins,2),dtype = float)
    n=0
    binsize=binedges[n+1]-binedges[n]
    while start+(binsize/2.0)<stop:
        #print start
        #print binsize
        #print stop
        #print n
        rebinned[n,0] = (start+(binsize/2.0))
        ind = np.where((x>start) & (x<start+binsize))
        rebinned[n,1] = np.mean(y[ind])
        start += binsize
        n+=1
        try:
            binsize=binedges[n+1]-binedges[n]
        except IndexError:
            break
    return rebinned
