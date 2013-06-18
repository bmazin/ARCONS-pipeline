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

#This program takes a template from "dataTemplate.py" which folds, bins, and averages the data and fits this template. A skewed gaussian is the best model for sdss j0926+3624 eclipses. The program can fit using a gaussian or polynomial with minor changes. Then use ConvolutionWithTemp.py to use the template.


template = np.load(tempName)
data = tempValues = template['template']
x = tempjd = template['jd']
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

#coeffOpt5 = fitPoly(x,data)
#print coeffOpt

Sparams = fitskewgaussian(x,data)
print Sparams

jdbin = (x[len(x)-1]-x[0])/float(len(x))
jd = np.arange(x[0],x[len(x)-1],jdbin)

FittedTemp = []

for point in jd:
    Value = skew(point, Sparams[0], Sparams[1], Sparams[2], Sparams[3], Sparams[4])
    FittedTemp.append(Value)

fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(x,gaussian(x, p))
ax.set_title(tempName)
#ax.plot(x,gauss(x, popt[0],popt[1],popt[2],popt[3]), label="gauss")
ax.plot(x,skew(x, Sparams[0], Sparams[1], Sparams[2], Sparams[3], Sparams[4]), label="skewed gaussian")
#ax.plot(x,poly5(x,coeffOpt5[0],coeffOpt5[1],coeffOpt5[2],coeffOpt5[3],coeffOpt5[4],coeffOpt5[5]),label="ploy5")
ax.plot(x,data,'k')
ax.plot(jd,FittedTemp,'r.')
plt.legend()
plt.show()

np.savez(FittedtempName,template=FittedTemp,jd=jd) 

print FittedtempName

