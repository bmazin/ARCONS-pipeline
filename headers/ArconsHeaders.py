#!/usr/bin/env python
# encoding: utf-8
"""
ArconsHeaders.py

Created by Ben Mazin on 2012-07-26.
Copyright (c) 2012 . All rights reserved.

This header file contains the classes that hold the reduced photon list in HDF5 files.

Updated, 4/19/2013, JvE.

"""

import sys
import os
from tables import *
from cosmic import TimeMask
 
class PhotonList(IsDescription):
    """The pytables derived class that stored reduced photons on disk.
       Each seperate observation file to be combined together should have its own separate PhotonList table in the h5 file.
       Note that this contains 3.3x the number of bytes as the raw photon packets, so expect a final file 3.3x larger than the original observation.
    """
    Xpix = UInt8Col()               # Xpixel the photon was recorded in, with the bottom left pixel being (0,0)
    Ypix = UInt8Col()               # Ypixel the photon was recorded in, with the bottom left pixel being (0,0)
    RA = Float32Col()               # Offset in arcseconds in RA from the observation center
    RAErr = Float32Col()            # Error in RA
    Dec = Float32Col()              # Offset in arcseconds in Dec from the observation center
    DecErr = Float32Col()           # Error in dec
    ArrivalTime = Float32Col()      # Time in seconds since the beginning of the observation
    Wavelength = Float32Col()       # Wavelength of the photon in Angstroms        
    WaveError = Float32Col()        # Estimated 1-sigma Wavelength error in Angstroms 
    TotalWeight = Float32Col()      # Combined weight of this photon.  Weight of 5.0 means the chance of this pixel getting from the top of the atmosphere to the pixel was 20%
    TotalWeightErr = Float32Col()   # Error in combined weight
    FlatWeight = Float32Col()       # Flat cal. weighting for this photon.
    FlatWeightErr = Float32Col()    # Error
    FluxWeight = Float32Col()       # Flux cal. weighting for this photon.
    FluxWeightErr = Float32Col()    # Error
    Flag = UInt8Col()               # Photon quality flag.  0=good, 1=pixel has incomplete wavelength coverage, 2=...
    TimeMaskFlag = BoolCol()        # Set if photon is time-masked. True=bad, False=good.
    TimeMaskReason = EnumCol(TimeMask.timeMaskReason, 'none', base='uint8')     # Time-mask reason ('none' if not time-masked)

    
