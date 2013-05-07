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
    Xpix = UInt8Col(dflt=255)                  # Xpixel the photon was recorded in, with the bottom left pixel being (0,0). Default value is an obviously bad value, since the detector is a lot smaller than 255 pixels.
    Ypix = UInt8Col(dflt=255)                  # Ypixel the photon was recorded in, with the bottom left pixel being (0,0)
    RA = Float32Col(dflt=np.nan)               # Offset in arcseconds in RA from the observation center
    RAErr = Float32Col(dflt=np.nan)            # Error in RA
    Dec = Float32Col(dflt=np.nan)              # Offset in arcseconds in Dec from the observation center
    DecErr = Float32Col(dflt=np.nan)           # Error in dec
    ArrivalTime = Float32Col(dflt=np.nan)      # Time in seconds since the beginning of the observation
    Wavelength = Float32Col(dflt=np.nan)       # Wavelength of the photon in Angstroms        
    WaveError = Float32Col(dflt=np.nan)        # Estimated 1-sigma Wavelength error in Angstroms 
    TotalWeight = Float32Col(dflt=np.nan)      # Combined weight of this photon.  Weight of 5.0 means the chance of this pixel getting from the top of the atmosphere to the pixel was 20%
    TotalWeightErr = Float32Col(dflt=np.nan)   # Error in combined weight
    FlatWeight = Float32Col(dflt=np.nan)       # Flat cal. weighting for this photon.
    FlatWeightErr = Float32Col(dflt=np.nan)    # Error
    FlatFlag = UInt8Col(dflt=flags.flatCal['undetermined'])     # Flat cal. flag.
    FluxWeight = Float32Col(dflt=np.nan)       # Flux cal. weighting for this photon.
    FluxWeightErr = Float32Col(dflt=np.nan)    # Error
    FluxFlag = UInt8Col(dflt=flags.fluxCal['undetermined'])    # Flux cal. flag. 255=undetermined.
    #QualityFlag = UInt8Col(dflt=255)           # Photon quality flag.  0=good, 1=pixel has incomplete wavelength coverage, 2=...; 255 indicates undetermined.
    TimeMaskFlag = BoolCol(dflt=False)         # Set if photon is time-masked. True=bad, False=good.
    TimeMaskReason = EnumCol(TimeMask.timeMaskReason, 'none', base='uint8')     # Time-mask reason ('none' if not time-masked)

