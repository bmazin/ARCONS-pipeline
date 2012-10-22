from numpy import *
import numpy as np
def cpsAtLick(mag, filterWidth):
    diameter = 120.0 * 2.54 # in cm
    area = pi * diameter*diameter/4
    eAtmosphere = 0.98
    ePrimary = 0.80
    eSecondary = 0.80
    eFilter = 0.80
    eVacuumWindow = 0.9
    eDetector = .10
    eTotal = eAtmosphere*ePrimary*eSecondary*eFilter*eVacuumWindow*eDetector
    # Vega of Magnitude 0 is 1000 cps per cm**2 per angstrom at top of atmosphere
    fluxTopOfAtmosphere = 1000 * 10**(-0.4*mag) * area * filterWidth
    flux = fluxTopOfAtmosphere * eTotal
    return flux,eTotal,fluxTopOfAtmosphere


