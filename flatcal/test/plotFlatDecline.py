import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

date = 20140924
#flatLabel = '20121212-074700'
flatLabel = '20140924'
wvlStart = 7000
wvlStop = 9000
#filename = '/Scratch/dataProcessing/flatTests/{}.txt'.format(flatLabel)
filename = '/Scratch/dataProcessing/flatTests/{}_{}-{}.txt'.format(flatLabel,wvlStart,wvlStop)
#cols: tstamp  counts  rawFwhm illumFwhm   flatFwhm    shotFwhm    flatShotRatio
if date == 20121210:
    firstTwilightTstamp = 135016 #20121210
    astroTwilight = 131215
    nauticalTwilight = 134228
elif date == 20121211:
    firstTwilightTstamp = 134024 #20121211
    astroTwilight = 131255
    nauticalTwilight = 134309
elif date == 20140924:
    firstTwilightTstamp = 124254 #20121211
    astroTwilight = 121412
    nauticalTwilight = 124321
    
table = np.loadtxt(filename,skiprows=1)

timestamps = table[:,0]
counts = table[:,1]
raw = table[:,2]
illum = table[:,3]
flat = table[:,4]
shot = table[:,5]
ratio = table[:,6]
exptime = table[:,7]
counts = 1.*counts / exptime

timestamps = ['{:06d}'.format(int(x)) for x in timestamps]
dates = [datetime.strptime(s,'%H%M%S') for s in timestamps]

fig,ax = plt.subplots(1,1)
ax2 = ax.twinx()

ax.plot(dates,raw,'k-',label='raw',marker='.')
ax.plot(dates,illum,'g-',label='illum',marker='.')
ax.plot(dates,flat,'b-',label='flat',marker='.')
ax.plot(dates,shot,'m-',label='shot',marker='.')
ax2.plot(dates,counts,'r--',label='count rate',marker='.')

astroTwilightDate = datetime.strptime(str(astroTwilight),'%H%M%S')
nauticalTwilightDate = datetime.strptime(str(nauticalTwilight),'%H%M%S')
firstTwilightDate = datetime.strptime(str(firstTwilightTstamp),'%H%M%S')
ax.axvline(astroTwilightDate,linestyle='--',color='k')
ax.axvline(nauticalTwilightDate,linestyle='--',color='k')
#ax.plot(dates,ratio,'r--',label='ratio')
fig.autofmt_xdate()
ax.legend(loc='upper right')
ax2.legend(loc='upper left')
ax.set_title('Array Uniformity on {}'.format(date))
ax.set_ylabel('FWHM of counts (%)')
ax2.set_ylabel('count rate (cps)')
ax.set_xlabel('time')

plt.show()
