'''
Author: Alex Walter
Date: Dec 3, 2014

This code performs PSF fitting photometry on a target and any number of guide stars in the field
'''

import warnings
import tables
from tables import *
import os

from photometry.Photometry import Photometry
from photometry.plot3DImage import *
from util.ObsFile import ObsFile
from util.FileName import FileName
from util.fitFunctions import model_list
from util.mpfit import mpfit
from util.popup import *




def readPsfPhotometryFile(filename):
    photometryFile = tables.openFile(filename, mode='r')
    photometryTable = photometryFile.getNode(photometryFile.root,'PSFphotometry')._f_getChild('photometry')
    startTimes = photometryTable.col('time')
    intTimes = photometryTable.col('integration_time')
    flux = photometryTable.col('flux')
    parameters = photometryTable.col('parameters')
    perrors = photometryTable.col('perrors')
    redChi2 = photometryTable.col('redChi2')
    flag = photometryTable.col('flag')
    photometryFile.close()
    
    return {'startTimes': startTimes, 'intTimes': intTimes, 'flux': flux, 'parameters': parameters, 'perrors': perrors, 'redChi2': redChi2, 'flag': flag}
    

def writePsfPhotometryFile(fluxDict_list, im_dict, filename,verbose = False):
    if os.path.isfile(filename):
        warnings.warn("Overwriting photometry file: "+str(filename),UserWarning)
        
    photometryFile = tables.openFile(filename, mode='w')
    photometryGroup = photometryFile.createGroup(photometryFile.root, 'PSFphotometry', 'Table of fit parameters for PSF')
    
    num_params = int(len(fluxDict_list[0]['parameters']))
    num_stars = int(len(fluxDict_list[0]['flux']))
    PSFphotometry_Description = {
        "time"              : Float64Col(),             # Start time of images
        "integration_time"  : Float64Col(),             # integration time of individual images
        "flux"              : Float64Col(num_stars),    # array of flux values. [target_flux, ref0_flux, ...]
        "parameters"        : Float64Col(num_params),   # parameters used to fit data
        "perrors"           : Float64Col(num_params),   # the errors on the fits
        "redChi2"           : Float64Col(),             # reduced chi^2 of the fits
        "flag"              : UInt16Col()}              # flag to indicate if fit is good (0)  
        
    photometryTable = photometryFile.createTable(photometryGroup, 'photometry', PSFphotometry_Description, title='PSF fitting solution')
    
    for i in range(len(im_dict['startTimes'])):
        row = photometryTable.row
        row['time'] = im_dict['startTimes'][i]
        row['integration_time'] = im_dict['intTimes'][i]
        row['flux'] = fluxDict_list[i]['flux']
        row['parameters'] = fluxDict_list[i]['parameters']
        row['perrors'] = fluxDict_list[i]['mpperr']
        row['redChi2'] = fluxDict_list[i]['redChi2']
        row['flag'] = fluxDict_list[i]['flag']
        row.append()
        
    # flush the table's I/O buffer to write the data to disk
    photometryTable.flush()
    if verbose:
        print "Wrote to: "+filename
    # close the file, flush all remaining buffers
    photometryFile.close()


def fitData2D(data,errs,parameter_guess,parameter_lowerlimit,parameter_upperlimit,model,parameter_ties=None,verbose=False):
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
            parameter_ties - array of length parameter_guess
                           - [0,0,1,0,1,2,0,2,0] means parameters 3 and 5 should be tied and parameters 6 and 8 should be tied.
                           - don't tie parameters with flags <=0
            verbose - Show runtime comments
            **kwargs - keyword parameters for model
            
        Outputs:
            parameter_fit - array of parameters for best fit
            redchi2gauss2 - the reduced chi^2 of the fit
            mpperr - array of errors on fit parameters
    """
    #Create list of parameters in format that mpfit likes
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
    
    #Tie some parameters together if specified
    if parameter_ties!=None and len(parameter_ties)==len(parameter_guess):
        for p_flag in np.unique(np.asarray(parameter_ties)[np.where(np.asarray(parameter_ties)>0)]):
            p_to_tie = np.where(np.asarray(parameter_ties)==p_flag)[0]
            if len(p_to_tie)>1:
                for p in p_to_tie[1:]:
                    parinfo[p]['tied'] = 'p['+str(p_to_tie[0])+']'

    #put parameters for fitting function in dictionary
    fa = {'data':data,'err':errs}
    quiet=True

    #Run mpfit, catch annoying warnings
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

    #Parse out fitted parameters
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

    def __init__(self,image,centroid,expTime=None,model = 'multiple_2d_circ_gauss_func', verbose=False,showPlot=False):
        '''
        Inputs:
            image - 2D array of data (0 for dead pixel, shouldn't be any nan's or infs)
                  - Should be fully calibrated, dead time corrected, and scaled up to the effective integration time
            centroid - list of (col,row) tuples. The first tuple is the target location. The next are reference stars in the field 
            expTime - 2d array of effective exposure times (same size as image)
            model - model used for fit. Options in util.fitFunctions
            verbose - show error messages
            showPlot - show and pause after each PSF fit
        '''
        self.verbose=verbose
        self.showPlot=showPlot
        self.model=model
        
        super(PSFphotometry,self).__init__(image=image,centroid=centroid,expTime=expTime)


        
    def guess_parameters(self, aper_radius=9,tie_sigmas=True):
        '''
        Inputs:
            aper_radius - double or list of doubles of the same length as self.centroid. Number of pixels around the star to be used in estimating parameters. -1 for the whole array.
            tie_sigmas - By default, tells mpfit to tie sigmas for multiple stars to the same value
            
        Outputs:
            parameter_guess - List of best guesses for fit parameters
            parameter_lowerlimit - Lower limit for each parameter to constrain fit
            parameter_upperlimit - Upper limit
            parameter_ties - Parameters to be tied should have matching numbers in their indices in this array. 
                           - eg. [0,0,1,0,1] would tie parameters 3 and 5 together
                           - anything <=0 is ignored
                           - Should be same length as parameter_guess
                           - Or can return parameter_ties=None
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
            overallBkgd = np.percentile(self.image[np.where(np.isfinite(self.image) & (self.image > 0.))],bkgdPercentile)   #This doesn't seem to work for some reason...
            #print 'bkgd ',overallBkgd
            #overallBkgd=1.
            parameter_guess = [overallBkgd]
            parameter_lowerlimit = [0.0]
            parameter_upperlimit = [np.mean(self.image[np.where(np.isfinite(self.image) & (self.image > 0.))])]
            parameter_ties=[0.]
            
            for star_i in range(len(self.centroid)):
                #p_guess = [1.,self.centroid[star_i][0],self.centroid[star_i][1],0.01]
                
                try:
                    radius=aper_radius[star_i]
                except TypeError:
                    radius=aper_radius
                
                x_guess = self.centroid[star_i][0]
                x_ll = 0.
                x_ul = len(self.image[0])    #number of columns
                y_guess = self.centroid[star_i][1]
                y_ll = 0.
                y_ul = len(self.image)       #number of rows
                if radius > 0.:
                    x_ll = max(self.centroid[star_i][0] - radius, x_ll)
                    x_ul = min(self.centroid[star_i][0] + radius, x_ul)
                    y_ll = max(self.centroid[star_i][1] - radius, y_ll)
                    y_ul = min(self.centroid[star_i][1] + radius, y_ul)
                    
                pixLoc = np.where(np.isfinite(self.image))
                if radius > 0.:
                    x_arr=np.tile(range(len(self.image[0])),(len(self.image),1))
                    y_arr=np.tile(range(len(self.image)),(len(self.image[0]),1)).transpose()
                    d_arr = np.sqrt((x_arr - x_guess)**2 + (y_arr - y_guess)**2)
                    
                    pixLoc = np.where(np.isfinite(self.image) * d_arr<=radius)
                    
                amp_guess = np.amax(self.image[pixLoc]) - overallBkgd
                #amp_guess = 10000.
                #amp_ul = amp_guess + 5.*np.sqrt(amp_guess)
                #amp_ul = 2.*amp_guess
                amp_ul = None
                amp_ll = max(amp_guess - 3.*np.sqrt(amp_guess),0.)
                #amp_ll = amp_guess/2.
                amp_ll = 0.
                
                sig_guess = 1.8
                sig_ll = 0.8
                sig_ul = 2.5
                if radius > 0. and sig_ul > radius: sig_ul = radius
                
                p_guess = [amp_guess, x_guess, y_guess, sig_guess]
                p_ll = [amp_ll, x_ll, y_ll, sig_ll]
                p_ul = [amp_ul, x_ul, y_ul, sig_ul]
                parameter_guess+=p_guess
                parameter_lowerlimit+=(p_ll)
                parameter_upperlimit+=(p_ul)
                if tie_sigmas==True:
                    parameter_ties+=[0.,0.,0.,1]
                else:
                    parameter_ties+=[0.,0.,0.,0.]
            
                
        #print_guesses(parameter_guess, parameter_lowerlimit, parameter_upperlimit, parameter_guess)
        return parameter_guess,parameter_lowerlimit,parameter_upperlimit,parameter_ties

    def PSFfit(self,aper_radius=-1, tie_sigmas=True):
        '''
        Inputs:
            aper_radius - double or list of doubles of the same length as self.centroid. Number of pixels around the star to be used in estimating parameters. -1 for the whole array.
            tie_sigmas - By default, tells mpfit to tie sigmas for multiple stars to the same value
            
        Returns: Dictionary with keywords
            flux - array of flux values. [target_flux, ref0_flux, ...]
            parameters - parameters used for fit
            mpperr - error on parameters from mpfit
            redChi2 - reduced chi^2 of fit
            flag - flag indicating if fit failed. 0 means success
        '''
        self.image[np.invert(np.isfinite(self.image))]=0.
        errs = np.sqrt(self.image)
        #errs[np.where(self.image==0.)]=1.
        #errs[np.where(expTime==0.)]=np.inf
        errs[np.where(self.image==0.)]=np.inf
        #errs[0:25,:] = np.inf
        #errs[:,11:] = np.inf
        
        
        parameter_guess,parameter_lowerlimit,parameter_upperlimit,parameter_ties = self.guess_parameters(aper_radius=aper_radius,tie_sigmas=tie_sigmas)
        
        models = model_list[self.model](parameter_guess)(p=np.ones(len(parameter_guess)),data=self.image,return_models=True)
        guess = models[0]
        for m in models[1:]:
            guess+=m
        if self.showPlot:
            pop(plotFunc=lambda popupOb: plot3DImage(popupOb.fig,popupOb.axes,self.image,errs=errs,fit=guess),title="Guess")

        #p_guess = np.ones(len(parameter_guess))
        #p_ll = np.asarray(parameter_lowerlimit,dtype=np.float)/np.asarray(parameter_guess,dtype=np.float)
        #p_ul = np.asarray(parameter_upperlimit,dtype=np.float)/np.asarray(parameter_guess,dtype=np.float)
        parameter_fit, redchi2gauss2, mpperr = fitData2D(np.copy(self.image),np.copy(errs),parameter_guess,parameter_lowerlimit,parameter_upperlimit,self.model,parameter_ties=parameter_ties,verbose=self.verbose)

        models2 = model_list[self.model](parameter_fit)(p=np.ones(len(parameter_fit)),data=self.image,return_models=True)
        fitModelImg = models2[0]
        for m in models2[1:]:
            fitModelImg+=m
        if self.showPlot:
            pop(plotFunc=lambda popupOb: plot3DImage(popupOb.fig,popupOb.axes,self.image,errs=errs,fit=fitModelImg),title="Fitted")
        
        
        flag = 2.0*self.fitHitLimit(parameter_fit,parameter_guess,parameter_lowerlimit,parameter_upperlimit)
        #return self.getFlux(parameter_fit)
        return {'flux': self.getFlux(parameter_fit), 'parameters': parameter_fit, 'mpperr':mpperr,'redChi2':redchi2gauss2,'flag':flag,'fitModelImg':fitModelImg}

    def fitHitLimit(self,parameter_fit,parameter_guess,parameter_lowerlimit,parameter_upperlimit):
        """
            This function just checks if the fit is railed against one of its parameter limits
            
            Returns:
                True if fit has a parameter that hit its upper or lower limit --> Bad fit
                False otherwise
        """
        fixed_guess = (parameter_lowerlimit==parameter_upperlimit) * (parameter_lowerlimit!=np.asarray([None]*len(parameter_lowerlimit)))
        
        s1=np.sum((parameter_fit==parameter_lowerlimit)*(np.logical_not(fixed_guess)))
        s2=np.sum((parameter_fit==parameter_upperlimit)*(np.logical_not(fixed_guess)))
        s3=np.sum((parameter_fit==parameter_guess)*(np.logical_not(fixed_guess)))
        if s1>0 or s2>0 or s3>0:
            return True
        else:
            return False


    def getFlux(self,parameter_fit):
        '''
        Inputs:
            parameter_fit - parameters of fit to calculate flux from
        Returns:
            flux - array of flux values. [target_flux, ref0_flux, ...]
        '''
        flux=[]
        if self.model=='multiple_2d_circ_gauss_func':
            for i in range(len(parameter_fit[1:])/4):
                star_flux = 2.*np.pi*parameter_fit[1+4*i]*(parameter_fit[4+4*i])**2
                flux.append(star_flux)
            
        return np.asarray(flux)








