import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from timing.photonTiming import processTimestamps
import datetime
from astropy.time import Time
from matplotlib import dates

obsTimeStr = '2016-07-04 03:00:00'
obsTime = Time(obsTimeStr,scale='utc')
obsMjd = obsTime.mjd

#Binary orbit parameters
pb = 0.38196663639948138217 #orbit period [days]
pbdot = 2.0945492216604939659e-11 #time derivative of orbit period [days/day]
tasc = 51260.200390532708553 #reference time, time of ascending node [mjd]
#ascending node is phase 0, minimum brightness is phase 0.25, maximum is phase 0.75
minPhaseOffset = 0.25

phaseAtObs = (obsMjd-tasc)/pb
nearestAscendingNodeTime = Time(np.round(phaseAtObs)*pb+tasc,format='mjd',scale='utc')

timeNowHours = (obsMjd - nearestAscendingNodeTime.mjd)*24.

sampledPhase,sampledMagR = np.loadtxt('lc2.txt',unpack=True)
sampledPhase = np.append(sampledPhase,1-sampledPhase[::-1])+minPhaseOffset
sampledMagR = np.append(sampledMagR,sampledMagR[::-1])

lcInterpFunc = interp1d(x=sampledPhase,y=sampledMagR,kind='quadratic')

phases = np.linspace(0.0,1,100)+minPhaseOffset
lightCurve = lcInterpFunc(phases)
phases = np.concatenate([phases-1,phases])
lightCurve = np.concatenate([lightCurve,lightCurve])
hours = phases*pb*24.#+timeNowHours
dhours = [nearestAscendingNodeTime.datetime+datetime.timedelta(hours=hour) for hour in hours]
fhours = dates.date2num(dhours)

hourFormat = dates.DateFormatter('%H')
fig,ax = plt.subplots(1,1)
#ax.plot(sampledPhase,sampledMagR,'b')
ax.plot(fhours,lightCurve,'r')
ax.axvline(dates.date2num(obsTime.datetime))
ax.xaxis.set_major_locator(dates.HourLocator(interval=2))
ax.xaxis.set_major_formatter(hourFormat)
#ax.axvline(0)
ax.set_title('B1957+20 light curve')
ax.set_ylabel('R')
ax.set_ylim([25,19.4])
ax.set_xlabel(obsTime.iso+' (UTC)')

plt.show()
