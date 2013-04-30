#!/bin/python

'''
Author: Paul Szypryt		Date: March 6, 2013

Takes a psf fitted npz array and creates a Lomb-Scargle Periodogram. Outputs the eclipse period and
frequency and compares it to the expected period for the SDSS-J0926 object.  Plots the light curve 
and periodogram.  Optimal for low frequencies (<1 Hz).
'''

import numpy as np
import scipy as sp
from scipy.signal import spectral
import matplotlib.pyplot as plt
from util import utils

# Function that calculates the Poisson noise level for a given frequency.  Will eventually want to include cosmic rays.
def poissonLevel(frequency,countRate,deadtime,timeBinSize,frequencyNumber):
    poissonNoise = 2*(1-2*countRate*deadtime*(1-deadtime/(2*timeBinSize))) - 2*((frequencyNumber-1)/frequencyNumber)*countRate*deadtime*(deadtime/timeBinSize)*np.cos(2*np.pi*timeBinSize*frequency)
    # Cosmic ray term: 2*countRateVLE*countRate*(deadtimeVLE**2)*(np.sin(np.pi*deadtimeVLE*frequency)/(np.pi*deadtimeVLE*frequency))**2
    return poissonNoise

# Load relevant data from psf fitted npz array
FileName = '/home/pszypryt/oldLightCurves/sdss_data/20121208/Blue-Fit.npz'
#FileName = '/home/pszypryt/sdss_data/20121211/Blue-ImageStackFit.npz'
t = np.load(FileName)
params = t['params']
jd = t['jd']*86400
amps = params[:,1]
widths = params[:,4]

# Calculate light curve, cut out bad data sequences, and normalize
curve = amps*widths**2

totalCounts= curve.sum()
totalTime= (jd[len(jd)-1]-jd[0])*86400
countRate= totalCounts/totalTime

#scaled_curve = curve
scaled_curve = (curve-curve.mean())/curve.std()
time = jd

# Create array of frequencies to check for fourier components, function requires angular frequencies
freqs=np.linspace(1,4000,10000)/86400
#freqs = np.logspace(1,3.63,num=10000)
angular_freqs=2*np.pi*freqs

# Create periodogram using frequencies, times, and normalized curve
periodogram = spectral.lombscargle(time, scaled_curve, angular_freqs)

# Calculate eclipse period and frequency, and compare it to expected value
eclipse_period = 2*np.pi/(angular_freqs[np.argmax(periodogram)])
eclipse_frequency = 1/eclipse_period
expected_period = 0.01966127 # in days
print 'Eclipse period =',eclipse_period,'days.'
print 'Eclipse frequency =',eclipse_frequency, 'cycles/day.'
print 'Percent error = ' + str(100*(eclipse_period-expected_period)/expected_period) + '%'

# Create a figure with light curve in top plot and periodogram in bottom plot
fig = plt.figure()
# Plot light curve
ax = fig.add_subplot(211)
ax.plot(time,scaled_curve,'b.')
ax.set_title('Light Curve')
plt.xlabel('Julian Date')
plt.ylabel('Scaled Counts')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(16)
# Plot fourier transform
ax = fig.add_subplot(212)
ax.plot(freqs,periodogram)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title('Periodogram')
plt.xlabel('Frequency (Cycles/Day)')
plt.ylabel('Power')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(16)
plt.show()
