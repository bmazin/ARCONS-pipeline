import numpy as np
import matplotlib.pyplot as plt
from util import utils

t1 = np.load('/home/pszypryt/sdss_data/20121208/RedTest-Fit.npz')
t2 = np.load('/home/pszypryt/sdss_data/20121208/BlueTest-Fit.npz')
params1 = t1['params']
params2 = t2['params']
jd = t1['jd']
amps1 = params1[:,1]
amps2 = params2[:,1]
widths1 = params1[:,4]
xpos1 = params1[:,2]
ypos1 = params1[:,3]
widths2 = params2[:,4]
xpos2 = params2[:,2]
ypos2 = params2[:,3]
jd2 = (jd/0.01966127)%1.

fig = plt.figure()
ax = fig.add_subplot(211)
curve1 = amps1*widths1**2
curve2 = amps2*widths2**2

curve=curve2/curve1

curve1/=np.median(curve1)
curve2/=np.median(curve2)

#curve=curve2/curve1

ax.plot(jd,curve1,'r.',jd,curve2,'b.')
ax2 = fig.add_subplot(212)
ax2.plot(jd,curve,'k.')

plt.show()
