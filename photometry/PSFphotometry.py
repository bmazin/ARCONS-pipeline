#This code will perform PSF fitting photometry on a target and one guide star
from Photometry import Photometry
from plot3DImage import *

from util.ObsFile import ObsFile
from util.FileName import FileName
from util.fitFunctions import model_list




def fitData2D(data,errs,parameter_guess,parameter_lowerlimit,parameter_upperlimit,model,verbose=False):
    """
        Runs mpfit.py. Gets all the parameters in the right format, calls mpfit, parses output
        
        Inputs:
            data - 2d array of data (0 for dead pixel)
            errs - 2d array of errors (np.inf for dead pixel, 1 for pixel with no counts)
            parameter_guess - Best guess for parameters
            parameter_lowerlimit - Strict lower limit in parameter space
            parameter_upperlimit - Strict upper limit in parameter space
                - None means no limit
                - If lower and upper limits are the same then the parameter is fixed
            model - [String] model used to fit data. Must be contained in model_list. See /util/fitFunctions.py
            verbose - Show runtime comments
            
        Outputs:
            parameter_fit - array of parameters for best fit
            redchi2gauss2 - the reduced chi^2 of the fit
            mpperr - array of errors on fit parameters
    """
    parinfo = []
    for k in range(len(parameter_guess)):
        lowLimit=True
        highLimit=True
        if parameter_lowerlimit[k]==None: lowLimit=False
        if parameter_upperlimit[k]==None: highLimit=False

        fix_guess = False
        if parameter_lowerlimit[k]==parameter_upperlimit[k] and parameter_lowerlimit[k]!=None: fix_guess=True

        par = {'n':k,'value':parameter_guess[k],'limits':[parameter_lowerlimit[k], parameter_upperlimit[k]],'limited':[lowLimit,highLimit],'fixed':fix_guess,'error':0}
        parinfo.append(par)

    #errs = np.sqrt(z_Arr)                         # Poisson counts 
    #errs[np.where(np.invert(errs > 0.))] = 1.
    fa = {'data':data,'err':errs}
    quiet=True

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        #m = mpfit.mpfit(model_list[model], functkw=fa, ftol=1.e-15, xtol=1.e-20, parinfo=parinfo, maxiter=2000, quiet=quiet)
        m = mpfit.mpfit(model_list[model], functkw=fa, parinfo=parinfo, maxiter=1000, quiet=quiet)

    mpp = m.params                                #The fit params
    mpperr = m.perror
    chi2gauss = m.fnorm
    redchi2gauss2 = 1.0*chi2gauss/np.sum(np.invert(np.isnan(data)))

    if verbose:
        print "Status: "+str(m.status)+" after "+str(m.niter)+" iterations"
        print "mpperr: "+str(mpperr)
        print "reduced Chi^2: "+str(redchi2gauss2)
        #print "fnorm: "+str(m.fnorm)
        if mpperr==None:
            print m.errmsg

    parameter_fit=np.copy(np.asarray(parameter_guess))
    for k,p in enumerate(mpp):
        parinfo[k]['value'] = p
        parameter_fit[k]=p

    if verbose:
        print_guesses(parameter_guess, parameter_lowerlimit, parameter_upperlimit, parameter_fit)

    return parameter_fit, redchi2gauss2, mpperr

def print_guesses(parameter_guess, parameter_lowerlimit, parameter_upperlimit, parameter_fit):
    """
        Prints a nicely formatted list of parameters:
            lower_limit < guess < upper_limit --> best_fit
    """
    print 'parameters -'
    for k in range(len(parameter_guess)):
        print '\t'+str(k)+': '+str(parameter_lowerlimit[k])+' < '+str(parameter_guess[k])+' < '+ str(parameter_upperlimit[k]) +' --> '+str(parameter_fit[k])


class PSFphotometry(Photometry):

    def __init__(self,path='/Scratch/DisplayStack/RUN_TEMPLATE/TARGET_TEMPLATE',verbose=False,showPlot=False):
        '''
        Inputs:
            path - path to the display stack target info
            verbose - show error messages
            showPlot - show and pause after each PSF fit
        '''
        self.verbose=verbose
        self.showPlot=showPlot
        
        super(PSFphotometry,self).__init__(path)


    def getImage(self):
        '''
        This code needs to be implemented
        '''
        return self.testImage()
        
    def testImage(self):
        n = 10
        deadTime = deadTime=100.e-6
        firstSec = 0
        integrationTime = 10
        
        print self.obsFNs[10].obs()
        obs = ObsFile(self.obsFNs[10].obs())
        obs.loadBestWvlCalFile()
        obs.setWvlCutoffs(wvlLowerLimit=3500, wvlUpperLimit=8000)
        obs.loadFlatCalFile(FileName(obsFile=obs).flatSoln())
        obs.loadHotPixCalFile(self.obsFNs[10].timeMask())
        
        im_dict = obs.getPixelCountImage(firstSec=firstSec, integrationTime=integrationTime,weighted=True, fluxWeighted=False, getRawCount=False)
        im = im_dict['image']
        #Correct for dead time
        w_deadTime = 1.0-im_dict['rawCounts']*deadTime/im_dict['effIntTimes']
        im = im/w_deadTime
        
        return im,im_dict['effIntTimes']

    def PSFfit(self,image,expTime,target_centroid,*arg):
        '''
        Inputs:
            image - 2D array of data
            expTime - 2d array of effective exposure times (same size as image)
        '''
        
        print len(arg)









if __name__ == '__main__':
    path = '/Scratch/DisplayStack/PAL2014/HAT_P1'
    verbose=True
    showPlot=True
    
    phot = PSFphotometry(path=path,verbose=verbose,showPlot=showPlot)
    image,expTime = phot.getImage()
    pop(plotFunc=lambda fig,axes: plot3DImage(fig,axes,image))
    phot.PSFfit(image,expTime,(0.0,0.0),(0.0,1.0))












