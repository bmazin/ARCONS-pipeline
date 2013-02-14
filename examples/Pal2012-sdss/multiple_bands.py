import numpy as np
import matplotlib.pyplot as plt
from util import utils

t1 = np.load('/home/pszypryt/sdss_data/20121208/Blue-Fit.npz')
t2 = np.load('/home/pszypryt/sdss_data/20121208/Green-Fit.npz')
t3 = np.load('/home/pszypryt/sdss_data/20121208/Red-Fit.npz')
t4 = np.load('/home/pszypryt/sdss_data/20121208/IR-Fit.npz')

jd = t1['jd']
jd2 = (jd/0.01966127)%1.

params1 = t1['params']
params2 = t2['params']
params3 = t3['params']
params4 = t4['params']

amps1 = params1[:,1]
amps2 = params2[:,1]
amps3 = params3[:,1]
amps4 = params4[:,1]

widths1 = params1[:,4]
widths2 = params2[:,4]
widths3 = params3[:,4]
widths4 = params4[:,4]

xpos1 = params1[:,2]
xpos2 = params2[:,2]
xpos3 = params3[:,2]
xpos4 = params4[:,2]

ypos1 = params1[:,3]
ypos2 = params2[:,3]
ypos3 = params3[:,3]
ypos4 = params4[:,3]

curve1 = amps1*widths1**2
curve2 = amps2*widths2**2
curve3 = amps3*widths3**2
curve4 = amps4*widths4**2

#curve1/=2000
#curve2/=1000
#curve3/=1000
#curve4/=2000

fig = plt.figure()
ax = fig.add_subplot(111)

curve=curve2/curve1

scale_factor = np.median(curve)
#curve2/=scale_factor

#curve1/=np.median(curve1)
#curve2/=np.median(curve2)

ax.plot(jd,curve1,'b.',label='3000-5000')
ax.plot(jd,curve2,'g.',label='5000-6000')
ax.plot(jd,curve3,'r.',label='6000-7000')
#ax.plot(jd,curve4,'m.',label='7000-9000')
ax.legend(numpoints =1)
plt.xlabel('JD')
plt.ylabel('Photon Counts')

plt.show()
