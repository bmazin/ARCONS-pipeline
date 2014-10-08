import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

filename = '/Scratch/dataProcessing/flatTests/20121211.txt'
#cols: tstamp  counts  rawFwhm illumFwhm   flatFwhm    shotFwhm    flatShotRatio
firstTwilightTstamp = 134024
twilightExpTime = 60. #sec
obsExpTime = 300. #sec

table = np.loadtxt(filename,skiprows=1)

timestamps = table[:,0]
counts = table[:,1]
raw = table[:,2]
illum = table[:,3]
flat = table[:,4]
ratio = table[:,6]
counts[timestamps >= firstTwilightTstamp] /= twilightExpTime
counts[timestamps < firstTwilightTstamp] /= obsExpTime

timestamps = ['{:06d}'.format(int(x)) for x in timestamps]
dates = [datetime.strptime(s,'%H%M%S') for s in timestamps]

fig,ax = plt.subplots(1,1)
ax2 = ax.twinx()

ax.plot(dates,raw,'k-',label='raw',marker='.')
ax.plot(dates,illum,'g-',label='illum',marker='.')
ax.plot(dates,flat,'b-',label='flat',marker='.')
ax2.plot(dates,counts,'r--',label='count rate',marker='.')

astroTwilightDate = datetime.strptime('131255','%H%M%S')
nauticalTwilightDate = datetime.strptime('134309','%H%M%S')
firstTwilightDate = datetime.strptime(str(firstTwilightTstamp),'%H%M%S')
ax.axvline(astroTwilightDate,linestyle='--',color='k')
ax.axvline(nauticalTwilightDate,linestyle='--',color='k')
#ax.plot(dates,ratio,'r--',label='ratio')
fig.autofmt_xdate()
ax.legend(loc='upper right')
ax2.legend(loc='upper left')
ax.set_title('Array Uniformity on 20121211')
ax.set_ylabel('FWHM of counts (%)')
ax2.set_ylabel('count rate (cps)')
ax.set_xlabel('time')

plt.show()
