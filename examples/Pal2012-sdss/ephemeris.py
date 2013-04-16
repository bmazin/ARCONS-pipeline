import numpy as np
import matplotlib.pyplot as plt
from util import utils
import mpfit
import scipy.optimize as optimize

FileName = '/Scratch/dataProcessing/SDSS_J0926/AllData/Dec8SIfitpsfRed.npz'
NumFrames = 1700
IntTime = 3
TimeCutOff = 300

t = np.load(FileName)

params = t['params']
jd = t['jd']

amps = params[:,1]
widths = params[:,4]
#xpos = params[:,2]
#ypos = params[:,3]

jd = jd[:TimeCutOff]
amps = amps[:TimeCutOff]
widths = widths[:TimeCutOff]

fig = plt.figure()
ax = fig.add_subplot(111)
curve = amps*widths**2

ax.plot(jd,curve,'g')
ax.set_title(FileName)
#plt.show()

x=jd
data=curve

def gaussian(x, pars):
    center, width, height, back = pars
    width = float(width)
    return back - height*np.exp(-(((center-x)/width)**2)/2)
 
# Define an error function between data and a Gaussian
def errorfunction(params, data, x):
    errorfunction = data - gaussian(x, params)
    return errorfunction
 
# Find an optimal Guassian fit for the data, return parameters of that Gaussian
def fitgaussian(data,x):
    params=(x.mean(),2.*(x[10]-x[0]),data.max()-data.min(), 390.)
    p, cov_x, infodict, mesg, success = optimize.leastsq(errorfunction, params, args=(data, x), full_output=1)
    print cov_x, infodict
    print mesg
    print success
    print params
    return p

xpeakguess=np.where(data ==data.min())[0]+1
#print xpeakguess
#print jd[xpeakguess]
xfitstart=max([xpeakguess-50,0])
xfitend=min([xpeakguess+50,len(x)])

p=fitgaussian(data[xfitstart:xfitend],x[xfitstart:xfitend])
print p
#print p[0]

ax.plot(jd,gaussian(jd, p))
plt.show()

#Gaussian function
#def gauss(x, a, x0, sigma, b):
#    return a*np.exp(-(x-x0)**2/(2*sigma**2))+b
#p0=[-100,2456271.010,.001,450]

#ax.plot(jd,gauss(jd,p0[0],p0[1],p0[2],p0[3]))

#popt, pcov = curve_fit(gauss, jd, curve, p0)
#print popt
#print pcov

#ax.plot(jd,gauss(jd,popt[0],popt[1],popt[2],popt[3]),'r.')
#plt.show()

