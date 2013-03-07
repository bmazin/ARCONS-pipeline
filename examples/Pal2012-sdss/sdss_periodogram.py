import numpy as np
import scipy as sp
from scipy.signal import spectral
import matplotlib.pyplot as plt
from util import utils

FileName = '/home/pszypryt/sdss_data/20121208/Red-Fit.npz'

FoldPeriod = 0.01966127 #This is in fractions of a day
t = np.load(FileName)
params = t['params']
jd = t['jd']
amps = params[:,1]
widths = params[:,4]
xpos = params[:,2]
ypos = params[:,3]

curve = amps*widths**2
#curve = curve[0:390]
#jd = jd[0:390]

#curve /= np.median(curve)


time = jd

# Check to see if this scaling is done correctly
scaled_curve = (curve-curve.mean())/curve.std()

freqs=np.linspace(0.1,250,10000)
angular_freqs=2*np.pi*freqs
periodogram = spectral.lombscargle(time, scaled_curve, angular_freqs)
eclipse_period = 2*np.pi/(angular_freqs[np.argmax(periodogram)])
eclipse_frequency = 1/eclipse_period
print 'Eclipse period =',eclipse_period,'days.'
print 'Eclipse frequency =',eclipse_frequency, 'cycles/day.'

'''
# generates 100 evenly spaced points between 1 and 1000
time = np.linspace(0, 10, 100)

# computes the sine value of each of those points
mags = np.sin(time)

# scales the sine values so that the mean is 0 and the variance is 1 (the documentation specifies that this must be done)
scaled_mags = (mags-mags.mean())/mags.std()

# generates 1000 frequencies between 0.01 and 1
freqs = np.linspace(0.01, 1, 1000)

# computes the Lomb Scargle Periodogram of the time and scaled magnitudes using each frequency as a guess
periodogram = spectral.lombscargle(time, scaled_mags, freqs)

# returns the inverse of the frequence (i.e. the period) of the largest periodogram value
print 1/freqs[np.argmax(periodogram)]
'''

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
ax.set_title('Fourier Transform')
plt.xlabel('Frequency (Cycles/Day)')
plt.ylabel('Transform Component')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(16)
plt.show()
