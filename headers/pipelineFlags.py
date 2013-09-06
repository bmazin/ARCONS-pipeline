'''
Author: Julian van Eyken    Date: Apr 30 2013
Definitions of all data flags used by the pipeline.
Currently dictionaries to map flag descriptions to integer values.
May update to use pytables Enums at some point down the road....
'''

#Flat cal. flags:
flatCal = {
           'good':0,                #No flagging.
           'infWeight':1,           #Spurious infinite weight was calculated - weight set to 1.0
           'zeroWeight':2,          #Spurious zero weight was calculated - weight set to 1.0
           'belowWaveCalRange':10,  #Derived wavelength is below formal validity range of calibration
           'aboveWaveCalRange':11,  #Derived wavelength is above formal validity range of calibration
           'undefined':20,          #Flagged, but reason is undefined.
           'undetermined':99,       #Flag status is undetermined.
           }

#Flux cal. flags
fluxCal = {
           'good':0,                #No flagging.
           'infWeight':1,           #Spurious infinite weight was calculated - weight set to 1.0
           'LEzeroWeight':2,        #Spurious less-than-or-equal-to-zero weight was calculated - weight set to 1.0
           'nanWeight':3,           #NaN weight was calculated.
           'belowWaveCalRange':10,  #Derived wavelength is below formal validity range of calibration
           'aboveWaveCalRange':11,  #Derived wavelength is above formal validity range of calibration
           'undefined':20,          #Flagged, but reason is undefined.
           'undetermined':99        #Flag status is undetermined.
           }
