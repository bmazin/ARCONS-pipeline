#This code will perform PSF fitting photometry on a target and one guide star
import warnings

from Photometry import Photometry
from plot3DImage import *

from util.ObsFile import ObsFile
from util.FileName import FileName
from util.fitFunctions import model_list
from util.mpfit import mpfit




def fitData2D(data,errs,parameter_guess,parameter_lowerlimit,parameter_upperlimit,model,verbose=False):
    """
        Runs mpfit.py. Gets all the parameters in the right format, calls mpfit, parses output.
        
        This version scales your guess to 1.0, then fits for a percent difference (unitless). 
        Then it scales back (puts the units back) and returns the fitted parameters.  
        
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

        #p_ll = None if parameter_lowerlimit[k]==None else 1.0*parameter_lowerlimit[k]/parameter_guess[k]
        #p_ul = None if parameter_upperlimit[k]==None else 1.0*parameter_upperlimit[k]/parameter_guess[k]
        #par = {'n':k,'value':1.,'limits':[p_ll, p_ul],'limited':[lowLimit,highLimit],'fixed':fix_guess,'error':0}
        par = {'n':k,'value':parameter_guess[k],'limits':[parameter_lowerlimit[k],parameter_upperlimit[k]],'limited':[lowLimit,highLimit],'fixed':fix_guess,'error':0}
        parinfo.append(par)

    #print parinfo

    #errs = np.sqrt(z_Arr)                         # Poisson counts 
    #errs[np.where(np.invert(errs > 0.))] = 1.
    fa = {'data':data,'err':errs}
    quiet=True

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        #m = mpfit.mpfit(model_list[model], functkw=fa, ftol=1.e-15, xtol=1.e-20, parinfo=parinfo, maxiter=2000, quiet=quiet)
        #m = mpfit(model_list[model](parameter_guess), functkw=fa, parinfo=parinfo, maxiter=1000, quiet=quiet)
        m = mpfit(model_list[model](np.ones(len(parameter_guess))), functkw=fa, parinfo=parinfo, maxiter=1000, quiet=quiet)

    mpp = m.params                                #The fit params
    #mpperr = (np.asarray(m.perror)*np.asarray(parameter_guess)).tolist()
    mpperr = m.perror
    chi2gauss = m.fnorm
    dof = np.sum(np.isfinite(errs))     #degrees of freedom
    redchi2gauss2 = 1.0*chi2gauss/dof

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
        #parameter_fit[k]=p*parameter_guess[k]
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

    def __init__(self,path='/Scratch/DisplayStack/RUN_TEMPLATE/TARGET_TEMPLATE',model = 'multiple_2d_circ_gauss_func', verbose=False,showPlot=False):
        '''
        Inputs:
            path - path to the display stack target info
            verbose - show error messages
            showPlot - show and pause after each PSF fit
        '''
        self.verbose=verbose
        self.showPlot=showPlot
        self.model=model
        
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
        integrationTime = 30
        
        print self.obsFNs[10].obs()
        obs = ObsFile(self.obsFNs[10].obs())
        obs.loadBestWvlCalFile()
        obs.setWvlCutoffs(wvlLowerLimit=3500, wvlUpperLimit=8000)
        obs.loadFlatCalFile(FileName(obsFile=obs).flatSoln())
        #obs.loadHotPixCalFile(self.obsFNs[10].timeMask(),reasons=['hot pixel','cold pixel','dead pixel'])
        
        im_dict = obs.getPixelCountImage(firstSec=firstSec, integrationTime=integrationTime,weighted=True, fluxWeighted=False, getRawCount=False)
        im = im_dict['image']
        #Correct for dead time
        w_deadTime = 1.0-im_dict['rawCounts']*deadTime/im_dict['effIntTimes']
        im = im/w_deadTime
        #Correct for exposure time
        im = im*integrationTime/im_dict['effIntTimes']
        #Remove any funny values
        im[np.invert(np.isfinite(im))]=0.

        
        return im,im_dict['effIntTimes']
        
    def getCentroid(self):
        '''
        Should return a list of (col,row) tuples indicating the location of the stars in the field.
        The first location should be the target star. The rest are reference stars.
        
        Needs to be implemented to grab info from centroid file in /Scratch/DisplayStack/
        '''

        return self.testCentroid()
        
    def testCentroid(self):
        target = (6.,30.)
        ref = (31.,26.)
        
        #return [target]
        return [target,ref]
        
    def guess_parameters(self, image, centroid, aper_radius):
        '''
        Inputs:
            image - 2D image of data (0 for dead pixel, shouldn't be any nan's or infs)
            centroid - list of (col,row) tuples. The first tuple is the target location. The next are reference stars in the field
            radius - double or list of doubles of the same length as centroid. Number of pixels around the star to be used in estimating parameters. -1 for the whole array.
        '''
        if self.model=='multiple_2d_circ_gauss_func':
            #p[0] = background
            #p[1] = amplitude
            #p[2] = x_offset    column
            #p[3] = y_offset    row
            #p[4] = sigma
            #And so on for the 2nd and 3rd gaussians etc...
            #A+Be^-((x-xo)^2+(y-y0)^2)/2s^2 + Ce^-((x-x1)^2+(y-y1)^2)/2d^2 + ...
            
            bkgdPercentile = 30.0
            #overallBkgd = np.percentile(image[np.where(np.isfinite(image))],bkgdPercentile)
            overallBkgd=1.
            parameter_guess = [overallBkgd]
            parameter_lowerlimit = [0.0]
            parameter_upperlimit = [np.mean(image[np.where(np.isfinite(image))])]
            
            for star_i in range(len(centroid)):
                #p_guess = [1.,centroid[star_i][0],centroid[star_i][1],0.01]
                
                try:
                    radius=aper_radius[star_i]
                except TypeError:
                    radius=aper_radius
                
                x_guess = centroid[star_i][0]
                x_ll = 0.
                x_ul = len(image[0])    #number of columns
                y_guess = centroid[star_i][1]
                y_ll = 0.
                y_ul = len(image)       #number of rows
                if radius > 0.:
                    x_ll = max(centroid[star_i][0] - radius, x_ll)
                    x_ul = min(centroid[star_i][0] + radius, x_ul)
                    y_ll = max(centroid[star_i][1] - radius, y_ll)
                    y_ul = min(centroid[star_i][1] + radius, y_ul)
                    
                pixLoc = np.where(np.isfinite(image))
                if radius > 0.:
                    x_arr=np.tile(range(len(image[0])),(len(image),1))
                    y_arr=np.tile(range(len(image)),(len(image[0]),1)).transpose()
                    d_arr = np.sqrt((x_arr - x_guess)**2 + (y_arr - y_guess)**2)
                    
                    pixLoc = np.where(np.isfinite(image) * d_arr<=radius)
                    
                amp_guess = np.amax(image[pixLoc]) - overallBkgd
                #amp_guess = 10000.
                #amp_ul = amp_guess + 5.*np.sqrt(amp_guess)
                #amp_ul = 2.*amp_guess
                amp_ul = None
                amp_ll = max(amp_guess - 3.*np.sqrt(amp_guess),0.)
                #amp_ll = amp_guess/2.
                amp_ll = 0.
                
                sig_guess = 1.5
                sig_ll = 0.5
                sig_ul = 3.*sig_guess
                if radius > 0. and sig_ul > radius: sig_ul = radius
                
                p_guess = [amp_guess, x_guess, y_guess, sig_guess]
                p_ll = [amp_ll, x_ll, y_ll, sig_ll]
                p_ul = [amp_ul, x_ul, y_ul, sig_ul]
                parameter_guess+=p_guess
                parameter_lowerlimit+=(p_ll)
                parameter_upperlimit+=(p_ul)
                
        #print_guesses(parameter_guess, parameter_lowerlimit, parameter_upperlimit, parameter_guess)
        return parameter_guess,parameter_lowerlimit,parameter_upperlimit

    def PSFfit(self,image,expTime,centroid,aper_radius):
        '''
        Inputs:
            image - 2D array of data (0 for dead pixel, shouldn't be any nan's or infs)
                  - Should be fully calibrated, dead time corrected, and scaled up to the effective integration time
            expTime - 2d array of effective exposure times (same size as image)
        '''
        image[np.invert(np.isfinite(image))]=0.
        errs = np.sqrt(image)
        #errs[np.where(image==0.)]=1.
        #errs[np.where(expTime==0.)]=np.inf
        errs[np.where(image==0.)]=np.inf
        #errs[0:25,:] = np.inf
        #errs[:,11:] = np.inf
        
        
        parameter_guess,parameter_lowerlimit,parameter_upperlimit = self.guess_parameters(image,centroid,aper_radius)
        
        models = model_list[self.model](parameter_guess)(p=np.ones(len(parameter_guess)),data=image,return_models=True)
        guess = models[0]
        for m in models[1:]:
            guess+=m
        pop(plotFunc=lambda fig,axes: plot3DImage(fig,axes,image,errs=errs,fit=guess),title="Guess")

        #p_guess = np.ones(len(parameter_guess))
        #p_ll = np.asarray(parameter_lowerlimit,dtype=np.float)/np.asarray(parameter_guess,dtype=np.float)
        #p_ul = np.asarray(parameter_upperlimit,dtype=np.float)/np.asarray(parameter_guess,dtype=np.float)
        parameter_fit, redchi2gauss2, mpperr = fitData2D(np.copy(image),np.copy(errs),parameter_guess,parameter_lowerlimit,parameter_upperlimit,self.model,verbose=self.verbose)

        models2 = model_list[self.model](parameter_fit)(p=np.ones(len(parameter_fit)),data=image,return_models=True)
        guess2 = models2[0]
        for m in models2[1:]:
            guess2+=m
        pop(plotFunc=lambda fig,axes: plot3DImage(fig,axes,image,errs=errs,fit=guess2),title="Fitted")






if __name__ == '__main__':
    path = '/Scratch/DisplayStack/PAL2014/HAT_P1'
    verbose=True
    showPlot=True
    
    phot = PSFphotometry(path=path,verbose=verbose,showPlot=showPlot)
    image,expTime = phot.getImage()
    centroid = phot.getCentroid()
    
    
    phot.PSFfit(image,expTime,centroid,5.)












