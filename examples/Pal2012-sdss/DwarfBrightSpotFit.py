import numpy as np
import matplotlib.pyplot as plt
from util import utils
import mpfit
import scipy.optimize as optimize
from scipy import pi,sqrt,exp
from scipy.special import erf
import scipy.integrate as integrate

#This program tries to fit the white dwarf and bright spot component of the light curve. The bright spot is not modeled correctly. DwarfFit.py is a more sophisticated program that only fits the white dwarf component.

npzDataFile = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlue.npz'

#Dec8 eclipse 1: [0:300], TFp0=(1,.7,10.5,-2.5,1,2,1,-5,-1.5,1,1,2),    sigma= [10]*110+[5]*130+[10]*60
#Dec8 eclipse 2: [600:900], TFp0=(1,2,10.5,-2.5,1,1,.8,-3,-1.5,1,1,2),     sigma= [10]*60+[1]*130+[10]*110
cutoffmin = 600
cutoffmax = 900
DataFile = np.load(npzDataFile)
params = DataFile['params']
jd = DataFile['jd']
#params[:,0] are the height offset of the gaussian (so the number of background pulses)
amps = params[:,1]
widths = params[:,4]
xpos = params[:,2]
ypos = params[:,3]
curve = 2*pi*amps*widths**2

jd = jd[cutoffmin:cutoffmax]
x = (jd-jd[0])*1000
data = curve[cutoffmin:cutoffmax]/np.average(curve[cutoffmin:cutoffmax])*10


def zero_if_negative(x):
    if isinstance(x, (int, long, float, complex))==True:
        if x < 0:
            x=0
    else: 
        for item in range(len(x)):
             if x[item] < 0 or item < 50 or item>250:
                 x[item]=0
    return x

#theta =< pi/2 which implies mu = cos(theta) => 0
#epsilon is a limb-darkening coefficient
#trapezoidally weighted averaging?
#try singular-value decomposition (SVD)= all components contribute linearly to the light curve
def WhiteDwarf(x,a,b,c,amp,epsilon):
    theta=a*x+b
    mu=np.cos(theta)
    mu=zero_if_negative(mu)
    return amp*mu*(1-epsilon+epsilon*mu)+c # Intensity is proportional to mu*(1-epsilon+epsilon*mu)

# d is the distance along the line defining the bright spot
#beta and gamma are power-law exponents which allow some flexibility in how the brightness varies
#l is a scalelength
def BrightSpot(x,BSa,BSb,l,beta,gamma):
    d=BSa*x+BSb
    d=zero_if_negative(d)
    S = (d/l)**beta*exp(-(d/l)**gamma) #surface brightness of the elements of the bright spot
#    Rspot = l*(beta/gamma)**(1/gamma) #max surface brightness
    return S 

def BSIntegrate(x,width,BSa,BSb,c,BSamp,l,beta,gamma):
    BSIntensity = []
#    a0List = np.arange(x[0],x[len(x)-1],.1)
    a0List = x
    for a0 in a0List:
        a0min = a0-width
        I= integrate.quad(BrightSpot,a0min,a0,args=(BSa,BSb,l,beta,gamma))
        I=BSamp*I[0]+c
        BSIntensity.append(I)
    return BSIntensity

def TotalFunction(x,a,b,c,amp,epsilon,width,BSa,BSb,BSamp,l,beta,gamma): #12 parameters
    return WhiteDwarf(x,a,b,c,amp,epsilon)+BSIntegrate(x,width,BSa,BSb,0,BSamp,l,beta,gamma)

def fitTotalFunction(x,data):
    sigma= [10]*60+[1]*130+[10]*110
    TFpopt, TFpcov = optimize.curve_fit(TotalFunction, x, data, p0=TFp0,sigma=sigma)
    print TFp0
    return TFpopt

#      a, b, c,  amp,epsilon,width,BSa,BSb,BSamp,l,beta,gamma
TFp0=(1,2,10.5,-2.5,1,1,.8,-3,-1.5,1,1,2)

fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(x,BrightSpot(x, BSp0[0],BSp0[1],BSp0[2],BSp0[3],BSp0[4],BSp0[5],BSp0[6]),'r')
#ax.plot(a0List,BSIntensity/np.average(BSIntensity)*10,'k')


TFparams = fitTotalFunction(x,data)
print TFparams

#TFparams=TFp0
TF0 = TotalFunction(x, TFp0[0],TFp0[1],TFp0[2],TFp0[3],TFp0[4],TFp0[5],TFp0[6],TFp0[7],TFp0[8],TFp0[9],TFp0[10],TFp0[11])
TF = TotalFunction(x, TFparams[0],TFparams[1],TFparams[2],TFparams[3],TFparams[4],TFparams[5],TFparams[6],TFparams[7],TFparams[8],TFparams[9],TFparams[10],TFparams[11])
WD = WhiteDwarf(x, TFparams[0],TFparams[1],TFparams[2],TFparams[3],TFparams[4])

index_min=TF.argmin()
x_min=x[index_min]
eclipseTime = x_min/1000.+jd[0]
print eclipseTime

WDindex_min=WD.argmin()
WDx_min=x[WDindex_min]
WDeclipseTime = WDx_min/1000.+jd[0]
print WDeclipseTime

index_min0=TF0.argmin()
x_min0=x[index_min0]
eclipseTime0 = x_min0/1000.+jd[0]
print eclipseTime0

ax.plot([jd[60],jd[190]],[np.average(data)]*2,'o')
ax.set_title(npzDataFile + ' eclipse 1, Fitted Min Of Total Light Curve = %f\n "error" = White Dwarf Fit Min-Total Light Fit Min = %f'%(eclipseTime0,np.abs(WDeclipseTime-eclipseTime)))
ax.plot(jd,TF,'r', label = 'Total Light Curve Fitted')
ax.plot(jd,WD,'b', label = 'White Dwarf Light Curve Fitted')
ax.plot(jd,BSIntegrate(x, TFparams[5],TFparams[6],TFparams[7],TFparams[2],TFparams[8],TFparams[9],TFparams[10],TFparams[11]),'g', label = 'Bright Spot Light Curve Fitted')

ax.plot(jd,WhiteDwarf(x, TFp0[0],TFp0[1],TFp0[2],TFp0[3],TFp0[4]),'b--',label = 'dashed lines use initial parameters')
ax.plot(jd,BSIntegrate(x, TFp0[5],TFp0[6],TFp0[7],TFp0[2],TFp0[8],TFp0[9],TFp0[10],TFp0[11]),'g--')
ax.plot(jd,TF0,'r--')
ax.legend()
ax.plot(jd,data,'k')
plt.show()

