#!/bin/python

'''
Author: Paul Szypryt		Date: March 6, 2013

Takes a psf fitted npz array and creates a Lomb-Scargle Periodogram. Outputs the eclipse period and
frequency and compares it to the expected period for the SDSS-J0926 object.  Plots the light curve 
and periodogram.
'''

import numpy as np
import scipy as sp
from scipy.signal import spectral
import matplotlib.pyplot as plt
from util import utils


# Load relevant data from psf fitted npz array
FileName = '/home/pszypryt/sdss_data/20121208/Blue-Fit.npz'
t = np.load(FileName)
params = t['params']
jd = t['jd']
amps = params[:,1]
widths = params[:,4]

# Calculate light curve, cut out bad data sequences, and normalize
curve = amps*widths**2
scaled_curve = (curve-curve.mean())/curve.std()
#curve = curve[3100:4800]
#jd = jd[3100:4800]
time = jd

# Create array of frequencies to check for fourier components, function requires angular frequencies
freqs=np.linspace(0.1,250,10000)
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
plt.xlabel('Time (Days)')
plt.ylabel('Scaled Counts')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(16)
# Plot fourier transform
ax = fig.add_subplot(212)
ax.plot(freqs,periodogram)
ax.set_title('Periodogram')
plt.xlabel('Frequency (Cycles/Day)')
plt.ylabel('Transform Component')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(16)
plt.show()
