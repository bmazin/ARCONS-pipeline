
# encoding: utf-8
"""
ARCONSclasses.py

Originally created by Ben Mazin on 2011-05-04.
Copyright (c) 2011 . All rights reserved.
"""

import numpy as np
from tables import *

class RawPulse(IsDescription):
    """The pytables derived class that hold raw pulse data on the disk.
    """
    starttime = Float64Col()    # start time of pulse data
    samprate = Float32Col()     # sample rate of the data in samples/sec
    npoints = Int32Col()        # number of data points in the pulse
    f0 =  Float32Col()          # resonant frequency data was taken at 
    atten1 =  Float32Col()      # attenuator 1 setting data was taken at
    atten2 =  Float32Col()      # attenuator 2 setting data was taken at
    Tstart = Float32Col()       # temp data was taken at
        
    I = Float32Col(2000)        # I pulse data, up to 5000 points.
    Q = Float32Col(2000)

class PulseAnalysis(IsDescription):     # contains final template info
    flag = Int16Col()                   # flag for quality of template.  If this could be a bad template set > 0
    count = Float32Col()                # number of pulses going into this template
    pstart = Int16Col()                 # index of peak of template
    phasetemplate = Float64Col(2000)
    phasenoise = Float64Col(800)
    phasenoiseidx = Float64Col(800)
    #optfilt = Complex128(800)   
    
    # fit quantities
    trise = Float32Col()                # fit value of rise time                    
    tfall = Float32Col()                # fit value of fall time
    
    # optimal filter parameters
    coeff = Float32Col(100)             # coefficients for the near-optimal filter
    nparam = Int16Col()                 # number of parameters in the filter

class BeamMap(IsDescription):
    roach = UInt16Col()                 # ROACH board number (0-15) for now!
    resnum = UInt16Col()                # resonator number on roach board (corresponds to res # in optimal pulse packets)
    f0 = Float32Col()                   # resonant frequency of center of sweep (can be used to get group name)
    pixel = UInt32Col()                 # actual pixel number - bottom left of array is 0, increasing up
    xpos = Float32Col()                 # physical X location in mm
    ypos = Float32Col()                 # physical Y location in mm
    scale = Float32Col(3)               # polynomial to convert from degrees to eV 

class ObsHeader(IsDescription):
    target = StringCol(80)
    datadir = StringCol(80)             # directory where observation data is stored
    calfile = StringCol(80)             # path and filename of calibration file
    beammapfile = StringCol(80)         # path and filename of beam map file
    version = StringCol(80)
    instrument = StringCol(80)
    telescope = StringCol(80)
    focus = StringCol(80)
    parallactic = Float64Col()
    ra = Float64Col()
    dec = Float64Col()
    alt = Float64Col()
    az = Float64Col()
    airmass = Float64Col()
    equinox = Float64Col()
    epoch = Float64Col()
    obslat = Float64Col()
    obslong = Float64Col()
    obsalt = Float64Col()
    timezone = Int32Col()
    localtime = StringCol(80)
    ut = Float64Col()
    lst = StringCol(80)
    jd = Float64Col()
    platescl = Float64Col()
    exptime = Int32Col()
