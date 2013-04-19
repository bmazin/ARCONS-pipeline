import numpy as np
import matplotlib.pyplot as plt
from util import utils
import mpfit
import scipy.optimize as optimize
from scipy import pi,sqrt,exp
from scipy.special import erf

cutoff = 82
tempName = '/Scratch/dataProcessing/SDSS_J0926/AllData/Dec8BlueTruncTemplateUpdated560.npz'
FittedtempName = '/Scratch/dataProcessing/SDSS_J0926/AllData/Dec8BlueTruncTemplateFitted560.npz'
FittedtempNameM = '/Scratch/dataProcessing/SDSS_J0926/AllData/Dec8BlueTruncTemplateMirrorFitted560.npz'
template = np.load(tempName)
data = tempValues = template['template']
x = tempjd = template['jd']
#x = (x-x[0])*1000
data1 = np.hstack((data[:cutoff],data[:cutoff][::-1]))
data1 = data1[:len(x)]
print len(data)
#x = x[:cutoff]


def gauss(x, x0, sigma, a, b):
    return b-a*np.exp(-(x-x0)**2/(2*sigma**2))

def poly5(x, a0, a1, a2, a3, a4, a5):
    return a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4 + a5*x**5

def pdf(x):
    return np.exp(-x**2/2)

def cdf(x):
    return (1 + erf(x/sqrt(2)))

def skew(x,x0,sigma,a,b,c): #e=0,w=1,a=0):
    t = (x-x0)/sigma
    return b-a*pdf(t)*cdf(c*t)

 
# Find an optimal Guassian fit for the data, return parameters of that Gaussian
def fitgaussian(x,data):
    params=(x.mean(),2.*(x[10]-x[0]),data.max()-data.min(),1)
    popt, pcov = optimize.curve_fit(gauss, x, data, p0=params)
    print params
    return popt

def fitskewgaussian(x,data):
    Sparams=(x.mean(),2.*(x[10]-x[0]),data.max()-data.min(),1, 10)
    Spopt, Spcov = optimize.curve_fit(skew, x, data, p0=Sparams)
    print Sparams
    return Spopt

def fitPoly(x,data):
    coeff5=(-.1,0,1,.1,.1,.1)
    coeffOpt5, pcov5 = optimize.curve_fit(poly5, x, data, p0=coeff5)
    print coeff5
    print coeffOpt5
    return coeffOpt5
    
#popt=fitgaussian(x,data)
#print popt

MirrorPopt=fitgaussian(x,data1)
print MirrorPopt

#coeffOpt5 = fitPoly(x,data)
#print coeffOpt

Sparams = fitskewgaussian(x,data)
print Sparams

jdbin = (x[len(x)-1]-x[0])/float(len(x))
jd = np.arange(x[0],x[len(x)-1],jdbin)

FittedTemp = []
FittedTempM = []
for point in jd:
    Value = skew(point, Sparams[0], Sparams[1], Sparams[2], Sparams[3], Sparams[4])
    FittedTemp.append(Value)
    ValueM = gauss(point, MirrorPopt[0],MirrorPopt[1],MirrorPopt[2],MirrorPopt[3])
    FittedTempM.append(ValueM)

fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(x,gaussian(x, p))
ax.set_title(tempName)
#ax.plot(x,gauss(x, popt[0],popt[1],popt[2],popt[3]), label="gauss")
ax.plot(x,gauss(x, MirrorPopt[0],MirrorPopt[1],MirrorPopt[2],MirrorPopt[3]), label="Mirrored gauss")
ax.plot(x,skew(x, Sparams[0], Sparams[1], Sparams[2], Sparams[3], Sparams[4]), label="skewed gaussian")
#ax.plot(x,poly5(x,coeffOpt5[0],coeffOpt5[1],coeffOpt5[2],coeffOpt5[3],coeffOpt5[4],coeffOpt5[5]),label="ploy5")
ax.plot(x,data,'k')
ax.plot(jd,FittedTemp,'r.')
ax.plot(jd,FittedTempM,'k.')
#ax.plot(x,data1)
plt.legend()
plt.show()

np.savez(FittedtempName,template=FittedTemp,jd=jd) 
np.savez(FittedtempNameM,template=FittedTempM,jd=jd) 

print FittedtempName

#def gaussian(x, pars):
#    center, width, height, back = pars
#    width = float(width)
#    return back - height*np.exp(-(((center-x)/width)**2)/2)
 
# Define an error function between data and a Gaussian
#def errorfunction(params, data, x):
#    errorfunction = data - gaussian(x, params)
#    return errorfunction

#    p, cov_x, infodict, mesg, success = optimize.leastsq(errorfunction, params, args=(x,data), full_output=1)
#    print cov_x, infodict
#    print mesg
#    print success

#def poly8(x, a0, a1, a2, a3, a4, a5, a6, a7, a8):
#    return a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4 + a5*x**5 + a6*x**6 + a7*x**7 + a8*x**8
#def poly7(x, a0, a1, a2, a3, a4, a5, a6, a7):
#    return a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4 + a5*x**5 + a6*x**6 + a7*x**7
#def poly6(x, a0, a1, a2, a3, a4, a5, a6):
#    return a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4 + a5*x**5 + a6*x**6

#    coeff8=(-.1,0,1,.1,.1,.1,1,1,1)
#    coeffOpt8, pcov8 = optimize.curve_fit(poly8, x, data, p0=coeff8)
#    coeff7=(-.1,0,1,.1,.1,.1,1,1)
#    coeffOpt7, pcov7 = optimize.curve_fit(poly7, x, data, p0=coeff7)
#    coeff6=(-.1,0,1,.1,.1,.1,1)
#    coeffOpt6, pcov6 = optimize.curve_fit(poly6, x, data, p0=coeff6)
#    print coeff8
#    print coeffOpt8
#    print coeff7
#    print coeffOpt7
#    print coeff6
#    print coeffOpt6

#ax.plot(x,poly8(x,coeffOpt8[0],coeffOpt8[1],coeffOpt8[2],coeffOpt8[3],coeffOpt8[4],coeffOpt8[5],coeffOpt8[6],coeffOpt8[7],coeffOpt8[8]),label="ploy8")
#ax.plot(x,poly7(x,coeffOpt7[0],coeffOpt7[1],coeffOpt7[2],coeffOpt7[3],coeffOpt7[4],coeffOpt7[5],coeffOpt7[6],coeffOpt7[7]),label="ploy7")
#ax.plot(x,poly6(x,coeffOpt6[0],coeffOpt6[1],coeffOpt6[2],coeffOpt6[3],coeffOpt6[4],coeffOpt6[5],coeffOpt6[6]),label="ploy6")
