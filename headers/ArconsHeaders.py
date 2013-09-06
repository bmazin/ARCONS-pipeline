#!/usr/bin/env python
# encoding: utf-8
"""
ArconsHeaders.py

Created by Ben Mazin on 2012-07-26.
Copyright (c) 2012 . All rights reserved.

This header file contains the classes that hold the reduced photon list in HDF5 files.

**Under construction, subject to change!!** 5/3/2012 JvE

"""

import sys
import os
import numpy as np
from tables import *
from headers import TimeMask
from headers import pipelineFlags as flags 


class PhotonList(IsDescription):
    """The pytables derived class that stored reduced photons on disk.
       Each seperate observation file to be combined together should have its own separate PhotonList table in the h5 file.
       Note that this contains 3.3x the number of bytes as the raw photon packets, so expect a final file 3.3x larger than the original observation.
    """
    xPix = UInt8Col(dflt=255)                  # Xpixel the photon was recorded in, with the bottom left pixel being (0,0). Default value is an obviously bad value, since the detector is a lot smaller than 255 pixels.
    yPix = UInt8Col(dflt=255)                  # Ypixel the photon was recorded in, with the bottom left pixel being (0,0)
    xyPix = UInt16Col(dflt=65535)              # Packing of X and Y pixel into a single value for the purposes of query efficiency.
    ra = Float64Col(dflt=np.nan)               # Offset in arcseconds in RA from the observation center
    raErr = Float32Col(dflt=np.nan)            # Error in RA
    dec = Float64Col(dflt=np.nan)              # Offset in arcseconds in Dec from the observation center
    decErr = Float32Col(dflt=np.nan)           # Error in dec
    ha = Float32Col(dflt=np.nan)               # Hour angle of detector at time of photon arrival
    astrometryFlag = UInt8Col(dflt=255)        # Place holder for now
    arrivalTime = Float32Col(dflt=np.nan)      # Time in seconds since the beginning of the observation
    wavelength = Float32Col(dflt=np.nan)       # Wavelength of the photon in Angstroms        
    waveError = Float32Col(dflt=np.nan)        # Estimated 1-sigma Wavelength error in Angstroms 
    waveFlag = UInt8Col(dflt=flags.waveCal['undetermined'])  #Any wavelength calibration flags associated with the photon....
    totalWeight = Float32Col(dflt=np.nan)      # Combined weight of this photon.  Weight of 5.0 means the chance of this pixel getting from the top of the atmosphere to the pixel was 20%
    totalWeightErr = Float32Col(dflt=np.nan)   # Error in combined weight
    flatWeight = Float32Col(dflt=np.nan)       # Flat cal. weighting for this photon.
    flatWeightErr = Float32Col(dflt=np.nan)    # Error
    flatFlag = UInt8Col(dflt=flags.flatCal['undetermined'])     # Flat cal. flag.
    fluxWeight = Float32Col(dflt=np.nan)       # Flux cal. weighting for this photon.
    fluxWeightErr = Float32Col(dflt=np.nan)    # Error
    fluxFlag = UInt8Col(dflt=flags.fluxCal['undetermined'])    # Flux cal. flag. 255=undetermined.
    #QualityFlag = UInt8Col(dflt=255)           # Photon quality flag.  0=good, 1=pixel has incomplete wavelength coverage, 2=...; 255 indicates undetermined.
    timeMaskFlag = BoolCol(dflt=False)         # Set if photon is time-masked. True=bad, False=good.
    timeMaskReason = EnumCol(TimeMask.timeMaskReason, 'none', base='uint8')     # Time-mask reason ('none' if not time-masked)

