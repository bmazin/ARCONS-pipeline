import numpy as np
import matplotlib.pyplot as plt
from util import utils
import mpfit
import scipy.optimize as optimize
from scipy import pi,sqrt,exp
from scipy.special import erf
import scipy.integrate as integrate

npzDataFile = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueUpdated.npz'

cutoffmin = 300
cutoffmax = 1100
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

check = []
def zero_if_negative(x):
    if isinstance(x, (int, long, float, complex))==True:
        check.append(x)
        if x < 0:
            x=0
    else: 
        for item in range(len(x)):
             if x[item] < 0:
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

print len(x)
print len(check)

def TotalFunction(x,a,b,c,amp,epsilon,width,BSa,BSb,BSamp,l,beta,gamma): #12 parameters
    return WhiteDwarf(x,a,b,c,amp,epsilon)+BSIntegrate(x,width,BSa,BSb,0,BSamp,l,beta,gamma)

#p=(a,b,c,amp,epsilon,width,BSa,BSb,BSamp,l,beta,gamma)

def mpTotalFunction(p, fjac=None, x=None, y=None, err=None): #12 parameters
    model= WhiteDwarf(x,p[0],p[1],p[2],p[3],p[4])+BSIntegrate(x,p[5],p[6],p[7],0,p[8],p[9],p[10],p[11])
    status = 0
    return([status, (y-model)/err])

def mpfitTotalFunction():
    parinfo = [{'value':0., 'limited':[0,0], 'limits':[0.,0.]}]*12
    parinfo[10]['limited'][0] = 1
    parinfo[10]['limits'][0]  = 0.1
    parinfo[11]['limited'][0] = 1
    parinfo[11]['limits'][0]  = 0.1
    values = [1,.7,10.5,-2.5,1,2,1,-5,-1.5,1,1,2]
    for i in range(5): parinfo[i]['value']=values[i]

    errs = np.sqrt(data)                         # Poisson counts 
    errs[np.where(errs == 0.)] = 1.
    quiet=True

    fa = {'x':x,'y':data,'err':errs}

    m = mpfit.mpfit(mpTotalFunction, functkw=fa, parinfo=parinfo, maxiter=1000, quiet=quiet)
            
    mpp = m.params                                #The fit params
    mpperr = m.perror
    print mpperr
    chi2gauss = m.fnorm
    redchi2gauss2 = chi2gauss/len(x)
    return mpp


def fitTotalFunction(x,data):
    TFpopt, TFpcov = optimize.curve_fit(TotalFunction, x, data, p0=TFp0)
    print TFp0
    return TFpopt

#a,b,c,amp,epsilon,width,BSa,BSb,BSamp,l,beta,gamma
TFp0=(1,.7,10.5,-2.5,1,2,1,-5,-1.5,1,1,2)

fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(x,BrightSpot(x, BSp0[0],BSp0[1],BSp0[2],BSp0[3],BSp0[4],BSp0[5],BSp0[6]),'r')
#ax.plot(a0List,BSIntensity/np.average(BSIntensity)*10,'k')

TFparams = fitTotalFunction(x,data)
print TFparams

TF0 = TotalFunction(x, TFp0[0],TFp0[1],TFp0[2],TFp0[3],TFp0[4],TFp0[5],TFp0[6],TFp0[7],TFp0[8],TFp0[9],TFp0[10],TFp0[11])
TF = TotalFunction(x, TFparams[0],TFparams[1],TFparams[2],TFparams[3],TFparams[4],TFparams[5],TFparams[6],TFparams[7],TFparams[8],TFparams[9],TFparams[10],TFparams[11])

index_min=TF.argmin()
x_min=x[index_min]
eclipseTime = x_min/1000.+jd[0]
print eclipseTime

index_min0=TF0.argmin()
x_min0=x[index_min0]
eclipseTime0 = x_min0/1000.+jd[0]
print eclipseTime0

ax.set_title(npzDataFile + ' eclipse 1, Fitted Min Of Total Light Curve = %f'%(eclipseTime))
ax.plot(jd,TotalFunction(x, TFparams[0],TFparams[1],TFparams[2],TFparams[3],TFparams[4],TFparams[5],TFparams[6],TFparams[7],TFparams[8],TFparams[9],TFparams[10],TFparams[11]),'r', label = 'Total Light Curve Fitted')
ax.plot(jd,WhiteDwarf(x, TFparams[0],TFparams[1],TFparams[2],TFparams[3],TFparams[4]),'b', label = 'White Dwarf Light Curve Fitted')
ax.plot(jd,BSIntegrate(x, TFparams[5],TFparams[6],TFparams[7],TFparams[2],TFparams[8],TFparams[9],TFparams[10],TFparams[11]),'g', label = 'Bright Spot Light Curve Fitted')

ax.plot(jd,WhiteDwarf(x, TFp0[0],TFp0[1],TFp0[2],TFp0[3],TFp0[4]),'b--',label = 'dashed lines use initial parameters')
ax.plot(jd,BSIntegrate(x, TFp0[5],TFp0[6],TFp0[7],TFp0[2],TFp0[8],TFp0[9],TFp0[10],TFp0[11]),'g--')
ax.plot(jd,TotalFunction(x, TFp0[0],TFp0[1],TFp0[2],TFp0[3],TFp0[4],TFp0[5],TFp0[6],TFp0[7],TFp0[8],TFp0[9],TFp0[10],TFp0[11]),'r--')
ax.legend()
ax.plot(jd,data,'k')
#plt.show()








#mpparms = mpfitTotalFunction()
#print mpparms

#ax.plot(x,WhiteDwarf(x, mpparms[0],mpparms[1],mpparms[2],mpparms[3],mpparms[4]),'b:')
#ax.plot(x,BSIntegrate(x, mpparms[5],mpparms[6],mpparms[7],mpparms[2],mpparms[8],mpparms[9],mpparms[10],mpparms[11]),'g:')

#ax.plot(x,TotalFunction(x, mpparms[0],mpparms[1],mpparms[2],mpparms[3],mpparms[4],mpparms[5],mpparms[6],mpparms[7],mpparms[8],mpparms[9],mpparms[10],mpparms[11]),'r:')

#ax.plot(x,mpTotalFunction(mpparms),'k')
plt.show()


def fitWhiteDwarf(x,data):
    WDp0=(1,1,10,-1,1)
    popt, pcov = optimize.curve_fit(WhiteDwarf, x, data, p0=WDp0)
    print WDp0
    return popt

def fitBrightSpot(x,data):
    print BSp0
    BSpopt, BSpcov = optimize.curve_fit(BrightSpot, x, data, p0=BSp0)
    print BSpopt
    return BSpopt

def fitBSIntensity(x,data):
    print BSInt0
    BSIntopt, BSIntcov = optimize.curve_fit(BSIntegrate, x, data, p0=BSInt0)
    print BSpopt
    return BSpopt

#WDparams = fitWhiteDwarf(x,data)
#print WDparams

#BSparams = fitBSIntensity(x,data)

#BSparams = fitBrightSpot(x,data)
#print BSparams

#ax.plot(x,fitBSIntensity(x, BSparams[0],BSparams[1],BSparams[2],BSparams[3],BSparams[4],BSparams[5],BSparams[6]))
#BSp0=(1,-5.5,10.5,-1,1,1,2)
#BSInt0 = (2,)+BSp0
#BSIntensity=BSIntegrate(x,BSInt0[0],BSInt0[1],BSInt0[2],BSInt0[3],BSInt0[4],BSInt0[5],BSInt0[6],BSInt0[7])
