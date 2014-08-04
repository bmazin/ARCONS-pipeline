'''
Author: Alex Walter                         Date: March 18, 2014
Based on Danica Marsden's earlier code

Reads in a wavelength calibration file (blue, red and IR laser light
) that has a beam map solution, and makes a histogram of photon phase
amplitudes for each pixel. Bad pixels are flagged.  Gaussians are fit
to the histograms.

The peaks are fit with a polynomial and that gives a phase amplitude
<-> wavelength correspondance.

Returns the polynomial fit parameters for each pixel, the error in
Angstroms (the width of the Gaussian), and the start and stop wave-
lengths, as well as the array of flagged pixels.

To run inside another script:

paramFile = sys.argv[1]
wavelengthCal(paramFile)

where paramFile has the same parameters as e.g. waveCal.dict

Depends on ObsFile.py, FileName.py, hotPixels.py, waveCal_diagnostics.py, readDict.py, mpfit.py, smooth.py, fitFunctions.py

see README.txt for more info
'''

#!/bin/env python

import time
import sys, os
import tables
from tables import *
import numpy as np
from scipy.special import erfc
import matplotlib as mpl
import matplotlib.pylab as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.offsetbox as offsetbox

from util.ObsFile import ObsFile
from util.readDict import readDict
from util.FileName import FileName
import hotpix.hotPixels as hp
import smooth
from fitFunctions import *
import mpfit
from waveCal_diagnostics import *

def getCalFileNames(paramFile,startswith='cal_', endswith='.h5',getAll=False,**kwargs):
    """
    Reads the parameter file with filename paramFile and then finds the correct cal_*.h5 observations in the mkidDataDir directory
    """
    try:
        params = readDict(paramFile)
        params.readFromFile(paramFile)
    except:
        params = paramFile

    mkidDataDir=params['mkidDataDir']
    if mkidDataDir==None:
        mkidDataDir=os.getenv('MKID_DATA_DIR', default="/ScienceData")
    intermDir=params['intermDir']
    if intermDir is None:
        intermDir = os.getenv('INTERM_DIR', default="/Scratch")
    outdir=params['outdir']
    if outdir is None or outdir is '':
        outdir = '/waveCalSolnFiles'
    if params['filename']!=None:
        return [FileName(obsFile = params['filename'], mkidDataDir=mkidDataDir,intermDir=intermDir,wavecalSolnDir=outdir,**kwargs)], params

    #Load in cal files from param file
    run = params['run']
    if run==None:
        raise IOError
    sunsetDate = params['sunsetDate']
    tStampList = params['calTimeStamps']
    if getAll:
        sunsetDate=None
    walkPath = None
    if sunsetDate==None or tStampList==None:
        if startswith=='cal_' or startswith=='obs_' or startswith=='flat_':
            walkPath=mkidDataDir+os.sep+run
        else:
            walkPath=intermDir+outdir+os.sep+run
        if not startswith.startswith('master') and sunsetDate != None:
            walkPath+=os.sep+sunsetDate
    if walkPath==None:
        return [FileName(run=run, date=sunsetDate,tstamp=tStamp,mkidDataDir=mkidDataDir,intermDir=intermDir,wavecalSolnDir=outdir) for tStamp in tStampList], params
    else:
        calFNs=[]
        #print 'walking: ',walkPath
        for root,dirs,files in os.walk(walkPath):
            for f in files:
                if f.startswith(startswith) and f.endswith(endswith) and not f.endswith('_drift'+endswith):
                    calFNs.append(root+os.sep+f)
        return [FileName(obsFile=fn,mkidDataDir=mkidDataDir,intermDir=intermDir,wavecalSolnDir=outdir) for fn in calFNs], params


def fitData(xArr,yArr,parameter_guess,parameter_lowerlimit,parameter_upperlimit,model,cut_off_phase=None,make_plot=False,verbose=False):
    #model=str(model)
    ##Model
    model_list = {
        'parabola': parabola,
        'gaussian': gaussian,
        'fourgaussian': fourgaussian,
        'threegaussian_exp': threegaussian_exp,
        'threegaussian_exppow': threegaussian_exppow,
        'threegaussian_moyal': threegaussian_moyal,
        'threegaussian_power': threegaussian_power,
        'threegaussian_lorentzian': threegaussian_lorentzian}

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

    if cut_off_phase!=None:
        cut_off_index = np.argmin(np.abs(np.asarray(xArr)-cut_off_phase))
        xArr=xArr[:cut_off_index]
        yArr=yArr[:cut_off_index]

    errs = np.sqrt(yArr)                         # Poisson counts 
    errs[np.where(errs == 0.)] = 1.
    fa = {'x':xArr,'y':yArr,'err':errs}
    quiet=True

    #m = mpfit.mpfit(model_list[model], functkw=fa, ftol=1.e-15, xtol=1.e-20, parinfo=parinfo, maxiter=2000, quiet=quiet)
    m = mpfit.mpfit(model_list[model], functkw=fa, parinfo=parinfo, maxiter=1000, quiet=quiet)

    mpp = m.params                                #The fit params
    mpperr = m.perror
    chi2gauss = m.fnorm
    redchi2gauss2 = 1.0*chi2gauss/len(xArr)

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

    if make_plot:
        plt.figure()
        if len(xArr)<20 and len(xArr)>1:
            plt.plot(xArr,yArr,'.-')
            xArr=np.arange(xArr[0],xArr[-1],(xArr[2]-xArr[1])/10.)
        else:
            plt.plot(xArr,yArr,'-')
        modelArr=model_list[model](p=parameter_fit,x=xArr,return_models=True)
        fullModel=np.zeros(len(xArr))
        for model_part in modelArr:
            plt.plot(xArr,model_part,'g')
            fullModel = fullModel+model_part
        plt.plot(xArr,fullModel,'r')
        plt.xlim(-600,0)
        plt.ylim(0,500)
        if not debug:
            plt.show()

    return parameter_fit, redchi2gauss2, mpperr

def print_guesses(parameter_guess, parameter_lowerlimit, parameter_upperlimit, parameter_fit):
    print 'parameters -'
    for k in range(len(parameter_guess)):
        print '\t'+str(k)+': '+str(parameter_lowerlimit[k])+' < '+str(parameter_guess[k])+' < '+ str(parameter_upperlimit[k]) +' --> '+str(parameter_fit[k])

class waveCal:
    def __init__(self,calFN,params,save_pdf=True,verbose=False,debug=False):
        """
        opens cal file, prepares hot pixel mask

        Inputs:
            calFN - FileName object to analyze
            params - dictionary of parameters
            save_pdf - if True, save plots in PDF
            verbose - if False, suppresses comments
            debug - if False, catch errors
        """
        self.params=params
        self.save_pdf=save_pdf
        self.verbose=verbose
        self.debug = debug
        self.calFN=calFN
        self.run=calFN.run
        self.laserCalFile=ObsFile(calFN.cal())
        timeAdjustFileName=FileName(run=self.run).timeAdjustments()
        try:
            self.laserCalFile.loadTimeAdjustmentFile(timeAdjustFileName,verbose=verbose)
        except:
            pass
        if not os.path.exists(self.calFN.timeMask()):
            print 'Running hotpix for ',self.calFN.cal()
            hp.findHotPixels(self.calFN.cal(),self.calFN.timeMask(),fwhm=np.inf,useLocalStdDev=True,nSigmaHot=3.0,maxIter=10)
            print "Laser cal file pixel mask saved to %s"%(self.calFN.timeMask())
        self.laserCalFile.loadHotPixCalFile(self.calFN.timeMask())

        if not debug:
            self.setupOutDirs(calFN)

        self.n_rows = self.laserCalFile.nRow
        self.n_cols = self.laserCalFile.nCol

        #data for calsoln file
        self.calsoln_tableNames=['pixelrow','pixelcol','roach','pixelnum','polyfit','sigma','solnrange','wave_flag']
        self.pixRowArr=[]
        self.pixColArr=[]
        self.pixRoachArr=[]
        self.pixNumArr=[]
        self.polyFitArr=[]
        self.sigmaArr=[]
        self.solnRangeArr=[]
        self.waveFlagArr=[]

        #data for drift study
        self.drift_tableNames=['pixelrow','pixelcol','gaussparams','perrors']
        self.drift_pixRowArr=[]
        self.drift_pixColArr=[]
        self.drift_gaussParams=[]
        self.drift_perrors=[]

        if (not debug) and (self.save_pdf):
            self.pp = PdfPages(self.outpath+params['figdir']+os.sep+'calsol_'+calFN.tstamp+'_fits.pdf')
            mpl.rcParams['font.size'] = 4
            self.n_plots_x = 3
            self.n_plots_y = 4
            self.n_plots_per_page = self.n_plots_x*self.n_plots_y
            self.plotcounter = 0

        self.rescounter = 0

        self.non_alloc_pix=self.laserCalFile.getNonAllocPixels()
        #self.dead_pix=self.laserCalFile.getDeadPixels(weighted=False,getRawCount=True)             #This takes too long...


        self.sensitivity=[]
        self.noisecutoff=[]
        self.deltaNoise=[]
        self.blueSig=[]
        self.noiseSig=[]


    def __del__(self):
        pass

    def setupOutDirs(self,calFN):
        intermDir=self.params['intermDir']
        outdir=self.params['outdir']
        if intermDir is None or intermDir is '':
            intermDir = os.getenv('MKID_PROC_PATH', default="/Scratch")
        if outdir is None or outdir is '':
            outdir = '/waveCalSolnFiles'
        self.outpath=intermDir+outdir+os.sep+calFN.run+os.sep+calFN.date


        try:
            os.makedirs(self.outpath)
        except:
            pass
        try:
            os.mkdir(self.outpath+self.params['figdir'])
        except:
            pass
        try:
            os.mkdir(self.outpath+self.params['driftdir'])
        except:
            pass


    def finish_pixel(self, pixelrow,pixelcol,polyfit=[-1]*3,sigma=-1,solnrange=[-1]*2,fitparams=None,mpperr=None,failFlag=None):
        if failFlag==None or failFlag>13:
            raise RuntimeError("Invalid failFlag!")
        if self.verbose and failFlag!=0:
            print "Failed: "+str(failFlag)

        self.pixRowArr.append(pixelrow)
        self.pixColArr.append(pixelcol)
        pixel_name=self.laserCalFile.beamImage[pixelrow][pixelcol]
        self.pixRoachArr.append(int(pixel_name.split('/')[1][1:]))
        self.pixNumArr.append(int(pixel_name.split('/')[2][1:]))
        self.polyFitArr.append(polyfit)
        self.sigmaArr.append(sigma)
        self.solnRangeArr.append(solnrange)
        self.waveFlagArr.append(failFlag)

        if failFlag==0:
            self.drift_pixRowArr.append(pixelrow)
            self.drift_pixColArr.append(pixelcol)
            self.drift_gaussParams.append(fitparams)
            self.drift_perrors.append(mpperr)

    def write_waveCal(self):

        WaveCalSoln_Description = {
            "roach"     : UInt16Col(),      # ROACH board number
            "pixelnum"  : UInt16Col(),      # pixel number on the roach
            "pixelrow"  : UInt16Col(),      # physical x location - from beam map
            "pixelcol"  : UInt16Col(),      # physical y location 
            "polyfit"   : Float64Col(3),    # polynomial to convert from phase amplitude to wavelength float 64 precision
            "sigma"     : Float64Col(),     # 1 sigma (Gaussian width) in eV, for blue peak
            "solnrange" : Float32Col(2),    # start and stop wavelengths for the fit in Angstroms
            "wave_flag" : UInt16Col()}      # flag to indicate if pixel is good (0), unallocated (1), dead (2), or failed during wave cal fitting (2+)   

        cal_file_name = self.outpath+"/calsol_" + self.calFN.tstamp + '.h5'
        try:
            waveCalFile = tables.openFile(cal_file_name,mode='w')
        except:
            pass
        # Create a group that branches from the root node, save cal data table in this group
        calgroup = waveCalFile.createGroup(waveCalFile.root, 'wavecal', 'Table of calibration parameters for each pixel')

        # Create a table instance under group calgroup with node name "calsoln"
        # the WaveCalSoln class declared before is the description parameter to define the columns of the table
        caltable = waveCalFile.createTable(calgroup, 'calsoln', WaveCalSoln_Description, title='Wavelength Cal Table')

        for pix in range(len(self.pixNumArr)):
            row = caltable.row
            row['pixelrow'] = self.pixRowArr[pix]
            row['pixelcol'] = self.pixColArr[pix]
            row['roach'] = self.pixRoachArr[pix]
            row['pixelnum'] = self.pixNumArr[pix]
            row['wave_flag'] = self.waveFlagArr[pix]
            row['polyfit'] = self.polyFitArr[pix]
            row['sigma'] = self.sigmaArr[pix]
            row['solnrange'] = self.solnRangeArr[pix]
            #row['fitmodel'] = self.params['model_type']
            row.append()

        caltable.attrs.model_type = self.params['model_type']

        if self.verbose:
            print "Wrote to: "+cal_file_name

        # flush the table's I/O buffer to write the data to disk
        caltable.flush()
        # close the file, flush all remaining buffers
        waveCalFile.close()

    def write_waveCal_drift(self):
        num_params = len(self.drift_gaussParams[0])
        DriftObj_Description = {
            "pixelrow"      : UInt16Col(),                  # physical x location - from beam map
            "pixelcol"      : UInt16Col(),                  # physical y location 
            "gaussparams"   : Float64Col(num_params),        # parameters used to fit data
            "perrors"       : Float64Col(num_params)}        # the errors on the fits

        drift_file_name = self.outpath+self.params['driftdir']+os.sep+"calsol_" + self.calFN.tstamp + '_drift.h5'
        try:
            driftCalFile = tables.openFile(drift_file_name,mode='w')
        except:
            pass
        
        driftgroup = driftCalFile.createGroup(driftCalFile.root, 'params_drift', 'Table of parameters for drift study')
        drifttable = driftCalFile.createTable(driftgroup, 'driftparams', DriftObj_Description, title='Drift Params Table')

        #for pix in range(len(self.drift_pixRowArr)):
        for pix in range(self.rescounter):
            row = drifttable.row
            row['pixelrow']=self.drift_pixRowArr[pix]
            row['pixelcol']=self.drift_pixColArr[pix]
            row['gaussparams']=self.drift_gaussParams[pix]
            row['perrors']=self.drift_perrors[pix]
            row.append()

        if self.verbose:
            print "Wrote to: "+drift_file_name

        drifttable.attrs.model_type = self.params['model_type']

        # flush the table's I/O buffer to write the data to disk
        drifttable.flush()
        # close the file, flush all remaining buffers
        driftCalFile.close()

    def guessBluePeak(self,n_inbin_total,phase_bins_total):
        last_ind_temp = np.where(n_inbin_total>self.params['min_amp'])[0][-1]
        sigma4_guess=np.abs(phase_bins_total[last_ind_temp]/(self.params['threshold_sigma']))
        #print sigma4_guess
        window = 'hanning'
        windowsize = int(self.params['bin_smooth'])   
 
        parab_smooth = smooth.smooth(n_inbin_total, windowsize, window)
        smoothed_data = np.array(parab_smooth, dtype=float)
        first_ind=np.where(smoothed_data>self.params['min_amp'])[0][0]
        #last_ind = np.where(n_inbin_total>self.params['min_amp'])[0][-1]
        last_ind = np.argmax(n_inbin_total)
        if first_ind+50>=last_ind:
            raise RuntimeError("No Blue Peak, try lowering 'min_amp'")
        smoothed_data=smoothed_data[first_ind:last_ind]
        phase_bins=phase_bins_total[first_ind:last_ind]
        n_inbin=n_inbin_total[first_ind:last_ind]
        
        gradients=np.diff(smoothed_data)>0
        min_Ind=[]
        max_Ind=[]
        for k in range(len(gradients)-1):
            if gradients[k]==True and gradients[k+1]==False:
                #if smoothed_data[k]
                max_Ind.append(k)
            if gradients[k]==False and gradients[k+1]==True and len(max_Ind)>0:
                min_Ind.append(k)
        
        try:
            #print phase_bins[max_Ind]
            #print phase_bins[min_Ind]
            blueMaxPhase_guess=phase_bins[max_Ind[0]]
            smallest_min_after_first_max = min_Ind[np.argmin(smoothed_data[min_Ind])]
            blueMinPhase_guess=phase_bins[smallest_min_after_first_max]
            blueMaxAmp_guess=n_inbin[max_Ind[0]]
        except:
            if self.verbose:
                print "Can't find Blue Peak local maxima"
                print "\tmaximums: "+str(phase_bins[max_Ind])
                print "\tminimums: "+str(phase_bins[min_Ind])
            raise
        
        #blueMaxAmp_guess=smoothed_data[max_Ind[0]]
        #blueSigma_guess=np.abs(blueMaxPhase_guess/params['blueR_guess'])/2.0
        blueSigma_guess = sigma4_guess

        bluePeak_guess=[blueSigma_guess,blueMaxPhase_guess,blueMaxAmp_guess]
        bluePeak_lowerlimit=[blueSigma_guess/10.,blueMaxPhase_guess-np.abs(blueMaxPhase_guess-blueMinPhase_guess),2.0]
        bluePeak_upperlimit=[blueSigma_guess*10.,blueMaxPhase_guess+np.abs(blueMaxPhase_guess-blueMinPhase_guess),2.*blueMaxAmp_guess]

        return bluePeak_guess, bluePeak_lowerlimit, bluePeak_upperlimit, phase_bins[smallest_min_after_first_max]

    def check_IR_in_Noise(self,parameter_fit):
        IR_offset = parameter_fit[7]
        IR_amp = parameter_fit[8]
        model_list = {
            'parabola': parabola,
            'gaussian': gaussian,
            'fourgaussian': fourgaussian,
            'threegaussian_exp': threegaussian_exp,
            'threegaussian_exppow': threegaussian_exppow,
            'threegaussian_moyal': threegaussian_moyal,
            'threegaussian_power': threegaussian_power,
            'threegaussian_lorentzian': threegaussian_lorentzian}
        model=self.params['model_type']
        modelArr=model_list[model](p=parameter_fit,x=[IR_offset],return_models=True)
        noise_amp_at_IR = modelArr[-1][0]

        if noise_amp_at_IR >= IR_amp:
            return True
        else:
            return False

    def guessBlueRedIrPeaks(self,n_inbin,phase_bins,effIntTime,bluePeak_Fit):
        ##Model
        model_list = {
            'parabola': parabola,
            'gaussian': gaussian,
            'fourgaussian': fourgaussian,
            'threegaussian_exp': threegaussian_exp,
            'threegaussian_exppow': threegaussian_exppow,
            'threegaussian_moyal': threegaussian_moyal,
            'threegaussian_power': threegaussian_power,
            'threegaussian_lorentzian': threegaussian_lorentzian}
        model=self.params['model_type']

        ##Already found blue peak
        sigma1_guess=bluePeak_Fit[0]
        x_offset1_guess=bluePeak_Fit[1]
        amplitude1_guess=bluePeak_Fit[2]
        blue_sig_ind = np.argmin(np.abs(phase_bins - (x_offset1_guess + 1.5*sigma1_guess)))
        #print phase_bins[blue_sig_ind]
        

        ##Red Peak
        ##Guess sigmas are about the same and offset proportional to energy
        sigma2_guess=sigma1_guess
        x_offset2_guess=(x_offset1_guess)*self.params['bluelambda']/(1.0*self.params['redlambda'])
        Red_ind=np.argmin(np.abs(phase_bins-x_offset2_guess))
        red_amp_1=np.mean(n_inbin[Red_ind:Red_ind+10])
        amplitude2_guess = red_amp_1
        #correct the location if red peak is actually to left
        red_amp_2_ind=np.argmax(n_inbin[blue_sig_ind:Red_ind])+blue_sig_ind
        #print phase_bins[blue_sig_ind:Red_ind]
        #print phase_bins[red_amp_2_ind]
        red_amp_2=n_inbin[red_amp_2_ind]-2.0*np.sqrt(n_inbin[red_amp_2_ind])
        if red_amp_2 > red_amp_1:
            Red_ind = red_amp_2_ind
            x_offset2_guess = phase_bins[Red_ind]
            amplitude2_guess = red_amp_2
            if self.verbose:
                print "Red peak moved left! "+str((x_offset1_guess)*self.params['bluelambda']/(1.0*self.params['redlambda']))+' --> '+str(x_offset2_guess)
                print "Compensating IR peak: "+str((x_offset1_guess)*self.params['bluelambda']/(1.0*self.params['irlambda']))+' --> '+str((x_offset2_guess)*self.params['redlambda']/(1.0*self.params['irlambda']))
        else:
            if amplitude1_guess > amplitude2_guess:
                amplitude2_guess=amplitude1_guess
        red_sig_ind = np.argmin(np.abs(phase_bins - (x_offset2_guess + 1.0*sigma2_guess)))

        ##IR Peak
        sigma3_guess=sigma1_guess
        #x_offset3_guess=(x_offset1_guess)*self.params['bluelambda']/(1.0*self.params['irlambda'])
        x_offset3_guess=(x_offset2_guess)*self.params['redlambda']/(1.0*self.params['irlambda'])
        #print str((x_offset1_guess)*self.params['bluelambda']/(1.0*self.params['irlambda'])) + ' --> '+str(x_offset3_guess)
        IR_ind=np.argmin(np.abs(phase_bins-x_offset3_guess))
        if IR_ind+20>len(n_inbin)-10:
            IR_ind=max(IR_ind,len(n_inbin)-25)
        IR_amp_1=np.mean(n_inbin[IR_ind:IR_ind+20])
        amplitude3_guess=IR_amp_1
        #correct the amplitude if IR is actually to left
        if red_sig_ind < IR_ind:
            IR_amp_2=max(n_inbin[red_sig_ind:IR_ind])
            IR_amp_2=IR_amp_2-1.5*np.sqrt(IR_amp_2)
            if IR_amp_2 > IR_amp_1:
                amplitude3_guess = IR_amp_2
                if self.verbose:
                    print "IR peak moved left!"
        else:
            if amplitude1_guess > amplitude3_guess:
                amplitude3_guess=amplitude1_guess

        parameter_guess = [sigma1_guess,x_offset1_guess,amplitude1_guess,sigma2_guess,x_offset2_guess,amplitude2_guess,sigma3_guess,x_offset3_guess,amplitude3_guess]
        parameter_lowerlimit = [sigma1_guess*0.8,x_offset1_guess-0.5*sigma1_guess,amplitude1_guess*0.80,
                                sigma2_guess*0.3, x_offset2_guess-1.2*sigma2_guess,amplitude2_guess*0.5,
                                sigma3_guess*0.3, x_offset3_guess-1.0*sigma3_guess,amplitude3_guess*0.5]
        parameter_upperlimit = [sigma1_guess*1.05,x_offset1_guess+0.5*sigma1_guess,amplitude1_guess*1.5,
                                sigma2_guess*1.3, x_offset2_guess+1.2*sigma2_guess,amplitude2_guess*1.7,
                                sigma3_guess*1.3, x_offset3_guess+1.2*sigma3_guess,amplitude3_guess*1.6]

        if model == None:
            return parameter_guess, parameter_lowerlimit, parameter_upperlimit

        if self.params['noise_guess_lower']!=None and self.params['noise_guess_upper']!=None:
            parameter_lowerlimit=np.concatenate((parameter_lowerlimit,self.params['noise_guess_lower']))
            parameter_upperlimit=np.concatenate((parameter_upperlimit,self.params['noise_guess_upper']))
        else:
            ##Noise tail
            if model=='fourgaussian':    #gaussian
                scale_factor4_guess=np.abs(phase_bins[-1]/(self.params['threshold_sigma']))
                x_offset4_guess=0.0
                num_noise=effIntTime*self.params['sample_rate']-np.sum(n_inbin)*self.params['num_points_per_photon']
                amplitude4_guess=num_noise*np.sqrt(2.0/np.pi)/scale_factor4_guess*1.0/erfc((-self.params['threshold_sigma']*scale_factor4_guess-x_offset4_guess)/(np.sqrt(2.0)*scale_factor4_guess))

                parameter_guess=np.concatenate((parameter_guess,[scale_factor4_guess,x_offset4_guess,amplitude4_guess]))
                gauss_lower=[scale_factor4_guess*0.8,phase_bins[-1]+scale_factor4_guess,max(n_inbin)]
                gauss_upper=[scale_factor4_guess*1.2,None,None]
                parameter_lowerlimit=np.concatenate((parameter_lowerlimit,gauss_lower))
                parameter_upperlimit=np.concatenate((parameter_upperlimit,gauss_upper))
            elif model=='threegaussian_exp':  #exponential
                x_offset4_guess = phase_bins[-1]
                amplitude4_guess = max(n_inbin[-6:])
                scale_factor4_guess = self.params['noise_scale_factor']

                parameter_guess=np.concatenate((parameter_guess,[scale_factor4_guess,x_offset4_guess,amplitude4_guess]))
                exp_lower=[None,x_offset3_guess,amplitude4_guess*0.6]
                exp_upper=[None,phase_bins[-1],amplitude4_guess*1.2]
                parameter_lowerlimit=np.concatenate((parameter_lowerlimit,exp_lower))
                parameter_upperlimit=np.concatenate((parameter_upperlimit,exp_upper))
            elif model=='threegaussian_exppow':  #exponential power
                x_offset4_guess = phase_bins[-1]
                amplitude4_guess = max(n_inbin[-6:])
                scale_factor4_guess = self.params['noise_scale_factor']
                pow4_guess = 1.2

                parameter_guess=np.concatenate((parameter_guess,[scale_factor4_guess,x_offset4_guess,amplitude4_guess,pow4_guess]))
                exppow_lower=[None,x_offset3_guess,amplitude4_guess*0.6,1.0]
                exppow_upper=[None,phase_bins[-1],amplitude4_guess*1.2,1.99]
                parameter_lowerlimit=np.concatenate((parameter_lowerlimit,exppow_lower))
                parameter_upperlimit=np.concatenate((parameter_upperlimit,exppow_upper))
            elif model=='threegaussian_moyal':    #moyal
                scale_factor4_guess=np.abs(phase_bins[-1]/(self.params['threshold_sigma']))
                x_offset4_guess=0.0
                num_noise=effIntTime*self.params['sample_rate']-np.sum(n_inbin)*self.params['num_points_per_photon']
                amplitude4_guess=max(3000.,max(n_inbin)+100.)

                parameter_guess=np.concatenate((parameter_guess,[scale_factor4_guess,x_offset4_guess,amplitude4_guess]))
                moyal_lower=[None,None,max(n_inbin)]
                moyal_upper=[None,None,None]
                parameter_lowerlimit=np.concatenate((parameter_lowerlimit,moyal_lower))
                parameter_upperlimit=np.concatenate((parameter_upperlimit,moyal_upper))
            elif model=='threegaussian_power':  #power
                x_offset4_guess=min((x_offset2_guess+phase_bins[-1])/2.0,phase_bins[-6])
                scale_factor4_guess = 2.2
                amplitude4_guess = 1.0

                parameter_guess=np.concatenate((parameter_guess,[scale_factor4_guess,x_offset4_guess,amplitude4_guess]))
                #pow_lower=[0.0,None,0.0]
                #pow_upper=[None,phase_bins[-5],None]
                pow_lower=[scale_factor4_guess,None,0.0]
                pow_upper=[scale_factor4_guess,phase_bins[-5],None]

                parameter_lowerlimit=np.concatenate((parameter_lowerlimit,pow_lower))
                parameter_upperlimit=np.concatenate((parameter_upperlimit,pow_upper))
            elif model=='threegaussian_lorentzian':  #lorentzian
                scale_factor4_guess=np.abs(phase_bins[-1]/(self.params['threshold_sigma']))
                x_offset4_guess=0.0
                num_noise=effIntTime*self.params['sample_rate']-np.sum(n_inbin)*self.params['num_points_per_photon']
                amplitude4_guess=max(3000.,max(n_inbin)+100.)

                parameter_guess=np.concatenate((parameter_guess,[scale_factor4_guess,x_offset4_guess,amplitude4_guess]))
                lorentz_lower=[scale_factor4_guess*0.8,phase_bins[-1]+scale_factor4_guess,max(n_inbin)]
                lorentz_upper=[scale_factor4_guess*1.2,None,None]
                parameter_lowerlimit=np.concatenate((parameter_lowerlimit,lorentz_lower))
                parameter_upperlimit=np.concatenate((parameter_upperlimit,lorentz_upper))
            else:
                raise RuntimeError("No model: "+model)
        try:
            model_list[model](parameter_guess,x=[-500],return_models=True)
            if len(parameter_guess)!=len(parameter_lowerlimit) or len(parameter_guess)!=len(parameter_upperlimit):
                raise RuntimeError("Invalid parameter guesses!")
        except:
            if self.verbose:
                print "Invalid parameter guesses!"
                print_guesses(parameter_guess, parameter_lowerlimit, parameter_upperlimit, [None]*len(parameter_guess))
            raise

        #cut_off_ind = len(n_inbin)-10+np.argmax(n_inbin[-10:])
        #print 'At ground: '+str(phase_bins[-1])
        #print '2.5 sig threshold: '+str(-self.params['threshold_sigma']*sigma1_guess)
        #print 'cut off: '+str(phase_bins[cut_off_ind])
        cut_off_ind = len(n_inbin)-self.params['noise_fall']+np.argmax(n_inbin[-self.params['noise_fall']:])
        return parameter_guess,parameter_lowerlimit,parameter_upperlimit, phase_bins[cut_off_ind]


    def fitparabola(self,n_inbin,phase_bins,fit_params,make_plot=False):
        model_list = {
            'parabola': parabola,
            'gaussian': gaussian,
            'fourgaussian': fourgaussian,
            'threegaussian_exp': threegaussian_exp,
            'threegaussian_exppow': threegaussian_exppow,
            'threegaussian_moyal': threegaussian_moyal,
            'threegaussian_power': threegaussian_power,
            'threegaussian_lorentzian': threegaussian_lorentzian}
        model=self.params['model_type']

        wavelengths = [self.params['bluelambda'],self.params['redlambda'], self.params['irlambda']]     #Angstroms
        energies = [self.params['h'] * self.params['c'] / (x * self.params['ang2m']) for x in wavelengths]             #eV
        laser_amps=np.asarray(fit_params[[1,4,7]])
        #energies.append(0.)
        #laser_amps = np.concatenate((fit_params[[1,4,7]],[0.0]))

        # parameter_guess [constant, linear term (slope), quadratic term (perturbation to straight line)]
        parameter_guess = [0.0,(energies[0]-energies[1])*1.0/(laser_amps[0]-laser_amps[1]), -10.**-5.]
        parameter_lowerlimit=[None]*3
        parameter_upperlimit=[None]*3
        if fit_params[8]==0.0:      #No IR solution
            energies=np.delete(energies,2)
            laser_amps = np.delete(laser_amps,2)
            # fix quadratic term to zero so we only fit a straight line
            parameter_guess[-1] = 0.
            parameter_lowerlimit[-1] = 0.0
            parameter_upperlimit[-1] = 0.0        
        else:
            guess_blue_E_from_red_IR = (energies[1]-energies[2])*1.0/(laser_amps[1]-laser_amps[2])*(laser_amps[0]-laser_amps[1])+energies[1]
            if guess_blue_E_from_red_IR < energies[0]:      #check if curving up instead of down
                parameter_guess[-1]*=-1.0
        
        parameter_fit, redchi2gauss2, mpperr = fitData(laser_amps,energies,parameter_guess,parameter_lowerlimit,parameter_upperlimit,model='parabola',make_plot=make_plot,verbose=self.verbose)
        #mpperr gives 0 if parameter is at limit or fixed or not varied
        num_param_fails = np.sum(parameter_fit==parameter_guess)+np.sum(parameter_fit==parameter_lowerlimit) + np.sum(parameter_fit==parameter_upperlimit)
        if parameter_fit[2]==0.0:       # no quadratic term, just fitting line
            num_param_fails = 0
        if mpperr==None or num_param_fails>0.0:
            if self.verbose:
                print "Unable to fit parabola"
                print_guesses(parameter_guess,[None]*3,[None]*3,parameter_fit)
                print "energies: "+str(energies)
                print "laser_amps: "+str(laser_amps)
            raise RuntimeError("Parabola Fit Failure. mpperr: "+str(mpperr))

        blue_ind = np.argmin(np.abs(phase_bins-fit_params[1]))
        left_ind = np.where(n_inbin[:blue_ind]<self.params['min_amp'])[0][-1]

        try:
            if fit_params[8]!=0.0:
                start_ind=np.argmin(np.abs(phase_bins-fit_params[7]))
            else:
                start_ind=np.argmin(np.abs(phase_bins-fit_params[4]))
            modelArr=model_list[model](p=fit_params,x=phase_bins[start_ind:],return_models=True)
            noise_model=modelArr[-1]
            last_laser_model=modelArr[-2]
            if fit_params[8]==0.0:
                last_laser_model=modelArr[-3]
            right_ind = start_ind+np.argmin(np.abs(noise_model-last_laser_model))
            
        except:
            if self.verbose:
                print "Unable to get model"
            right_ind = len(phase_bins)-6
            raise

        #print left_ind
        #print phase_bins[left_ind]
        #print start_ind
        #print phase_bins[start_ind]
        #print right_ind
        #print phase_bins[right_ind]        

        e_fromphase=(parabola(parameter_fit,x=phase_bins[left_ind:right_ind],return_models=True))[0]
        lambda_fromphase = (params['h'] * params['c'] / params['ang2m']) / e_fromphase
    
        #solution range in Angstroms
        soln_range = np.asarray([lambda_fromphase[-1],lambda_fromphase[0]])
        #blue sigma in eV
        blue_peak = (parabola(parameter_fit,x=np.asarray([fit_params[1]]),return_models=True))[0][0]
        blue_sigma = np.abs(blue_peak - (parabola(parameter_fit,x=np.asarray([fit_params[1]+fit_params[0]]),return_models=True))[0][0])
        #print blue_sigma
        return parameter_fit, soln_range, blue_sigma

    def fitHitLimit(self,parameter_fit,parameter_guess,parameter_lowerlimit,parameter_upperlimit):
        fixed_guess = (parameter_lowerlimit==parameter_upperlimit) * (parameter_lowerlimit!=np.asarray([None]*len(parameter_lowerlimit)))
        s1=np.sum((parameter_fit==parameter_lowerlimit)*(np.logical_not(fixed_guess)))
        s2=np.sum((parameter_fit==parameter_upperlimit)*(np.logical_not(fixed_guess)))
        s3=np.sum((parameter_fit==parameter_guess)*(np.logical_not(fixed_guess)))
        if s1>0 or s2>0 or s3>0:
            return True
        else:
            return False

    def findWaveLengthSoln(self):

        print 'Starting file ', self.calFN.cal(), ' at time: ', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
        for i in range(self.n_rows):
            for j in range(self.n_cols):
#        for i in [12,13]:
#            for j in range(self.n_cols):

#        for i in [12]:
#            for j in [15,16]:
#            for j in range(3,19):
                self.findWaveLengthSoln_pix(i,j)

        self.plot_pix_pdf(None,None,None,None,None,None,None)       #Save any remaining PDF figures and close PDF file
        plt.close('all')                                            #Close any remaining matplotlib figures

        #self.makeDiagnositcPlots()
        
        print '\nFound ', self.rescounter, ' pixels with wavelength calibration solutions. Time: ', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())

    def findWaveLengthSoln_pix(self,i,j,fit_guesses=None):
        '''
            Finds wavelength solution for pixel i,j. Uses fit_guesses as first guess of parameters. If None, it automatically finds these
            
            fit_guesses: [sigma_blue,x_offset_blue,amplitude_blue,sigma_red,x_offset_red,amplitude_red,sigma_IR,x_offset_IR,amplitude_IR,scale_factor_noise,x_offset_noise,amplitude_noise,phase_cutoff]
        '''

        if self.verbose:
            print '('+str(i)+', '+str(j)+') --> '+self.laserCalFile.beamImage[i][j]
        if self.non_alloc_pix[i,j]==0:
            self.finish_pixel(i,j,failFlag=1)
            return False
        #if self.dead_pix[i,j]==0:
        #    self.failure(2)
        #    continue
        dataDict=self.laserCalFile.getTimedPacketList(i,j,timeSpacingCut=self.params['danicas_cut'])
        peakHeights=np.asarray(dataDict['peakHeights'])*1.0
        ## less than 'min_amp' per second average count rate
        if dataDict['effIntTime']==0.0 or len(peakHeights)<=(dataDict['effIntTime']*self.params['min_count_rate']):
            # See if it's hot the whole time, this might be a bad hotpixel mask
            self.laserCalFile.switchOffHotPixTimeMask()
            dataDict=self.laserCalFile.getTimedPacketList(i,j,timeSpacingCut=self.params['danicas_cut'])
            peakHeights=np.asarray(dataDict['peakHeights'])*1.0
            self.laserCalFile.switchOnHotPixTimeMask()
            if dataDict['effIntTime']==0.0 or len(peakHeights)<=(dataDict['effIntTime']*self.params['min_count_rate']):
                self.finish_pixel(i,j,failFlag=2)
                if self.verbose:
                    print "Dead Pixel"
                return False
        baselines=np.asarray(dataDict['baselines'])*1.0
        peakHeights-=baselines
        biggest_photon = int(min(peakHeights))
        n_inbin,phase_bins=np.histogram(peakHeights,bins=np.abs(biggest_photon),range=(biggest_photon,0))
        phase_bins=(phase_bins+(phase_bins[1]-phase_bins[0])/2.0)[:-1]
        try:
            last_ind = np.where(n_inbin>self.params['min_amp'])[0][-1]
        except IndexError:
            last_ind=len(n_inbin)-1
        ## Cut out all the zeros on the right
        n_inbin = n_inbin[:last_ind]
        phase_bins = phase_bins[:last_ind]
        
        ## Try to fit Blue peak first. Catch fit error
        if fit_guesses==None:
            try:
                parameter_guess,parameter_lowerlimit,parameter_upperlimit, cut_off_phase = self.guessBluePeak(n_inbin,phase_bins)
                bluePeak_Fit, redchi2blue, mpperr = fitData(phase_bins,n_inbin,parameter_guess,parameter_lowerlimit,parameter_upperlimit,model='gaussian',cut_off_phase=cut_off_phase,make_plot=self.debug,verbose=self.verbose)
            except:
                self.finish_pixel(i,j,failFlag=3)
                self.plot_pix_pdf(i,j,n_inbin,phase_bins,None,None,3)
                if self.verbose:
                    print "System Info: "+str(sys.exc_info()[0])
                if self.debug:
                    raise
                else:
                    return False

            ## If the fit didn't work then cut solution
            if (mpperr == None):
                if self.verbose:
                    print "Error fitting blue peak"
                self.finish_pixel(i,j,failFlag=4)
                self.plot_pix_pdf(i,j,n_inbin,phase_bins,bluePeak_Fit,redchi2blue,4)
                return False

            ## Cut if the reduced chi^2 indicates a bad fit
            if redchi2blue==None or redchi2blue>self.params['max_chi2_blue']:
                if self.verbose:
                    print "Chi^2 too high: "+str(redchi2blue)
                self.finish_pixel(i,j,failFlag=5)
                self.plot_pix_pdf(i,j,n_inbin,phase_bins,bluePeak_Fit,redchi2blue,5)
                return False
        
        ## Try to fit blue, red, IR, noise peaks concurrently. Catch any errors during fit
        try:
            if fit_guesses==None:
                parameter_guess,parameter_lowerlimit,parameter_upperlimit, cut_off_phase = self.guessBlueRedIrPeaks(n_inbin,phase_bins,dataDict['effIntTime'],bluePeak_Fit)
            else:
                parameter_guess = fit_guesses[:-1]
                parameter_lowerlimit = [None]*(len(fit_guesses)-1)
                parameter_upperlimit = [None]*(len(fit_guesses)-1)
                cut_off_phase = fit_guesses[-1]
            parameter_fit, redchi2, mpperr = fitData(phase_bins,n_inbin,parameter_guess,parameter_lowerlimit,parameter_upperlimit,model=self.params['model_type'],cut_off_phase=cut_off_phase,make_plot=self.debug,verbose=self.verbose)
        except:
            self.finish_pixel(i,j,failFlag=6)
            self.plot_pix_pdf(i,j,n_inbin,phase_bins,bluePeak_Fit,None,6)
            if self.verbose:
                print "System Info: "+str(sys.exc_info()[0])
            if self.debug:
                raise
            else:
                return False

        ## Cut if not a good fit and try again without IR peak
        try_2_peak_fit = False
        if mpperr == None:
            if self.verbose:
                print "mpperr == None. Now fitting without IR peak"
            try_2_peak_fit = True
        elif self.fitHitLimit(parameter_fit,parameter_guess,parameter_lowerlimit,parameter_upperlimit):
            if self.verbose:
                print "Blue, Red, IR, or Noise peak hit fit limits. Now fitting without IR peak"
            try_2_peak_fit = True
        #elif np.sum((parameter_fit==parameter_lowerlimit)) + np.sum((parameter_fit==parameter_upperlimit)) != 0.0:
        ##elif min(mpperr[6:])<=0.0:
        #    if self.verbose:
        #        print "Blue, Red, IR, or Noise peak hit fit limits. Now fitting without IR peak"
        #    try_2_peak_fit = True
        #elif np.sum((parameter_fit==parameter_guess)) != 0.0:
        #    if self.verbose:
        #        print "Some parameters weren't varied. Now fitting without IR peak"
        #    try_2_peak_fit = True
        elif redchi2==None:
            if self.verbose:
                print "reduced chi^2 == None. Now fitting without IR peak"
            try_2_peak_fit = True
        elif redchi2>self.params['max_chi2_all']:
            if self.verbose:
                print "reduced chi^2 is too large: "+str(redchi2)+". Now fitting without IR peak"
            try_2_peak_fit = True
        elif parameter_fit[7]>phase_bins[-5]:
            if self.verbose:
                print "IR peak location outside of phase range. Now fitting without IR peak"
            try_2_peak_fit = True
        elif self.check_IR_in_Noise(parameter_fit)==True:
            if self.verbose:
                print "IR peak is smaller than noise tail. Now fitting without IR peak"
            try_2_peak_fit = True
        if try_2_peak_fit:
            try:
                parameter_guess[8]=0        #set amplitude3=0
                parameter_lowerlimit[6:9]=parameter_guess[6:9]  #fix parameters
                parameter_upperlimit[6:9]=parameter_guess[6:9]
                parameter_fit, redchi2, mpperr = fitData(phase_bins,n_inbin,parameter_guess,parameter_lowerlimit,parameter_upperlimit,model=self.params['model_type'],cut_off_phase=cut_off_phase,make_plot=self.debug,verbose=self.verbose)
            except:
                self.finish_pixel(i,j,failFlag=7)
                self.plot_pix_pdf(i,j,n_inbin,phase_bins,parameter_fit,redchi2gauss2,7)
                if self.verbose:
                    print "System Info: "+str(sys.exc_info()[0])
                if self.debug:
                    raise
                else:
                    return False

            #Check if failed again
            if mpperr == None:
                if self.verbose:
                    print "mpperr == None. Fit Failed"
                self.finish_pixel(i,j,failFlag=8)
                self.plot_pix_pdf(i,j,n_inbin,phase_bins,parameter_fit,redchi2,8)
                return False
            elif self.fitHitLimit(parameter_fit,parameter_guess,parameter_lowerlimit,parameter_upperlimit):
                if self.verbose:
                    print "Blue, Red, or Noise peak hit fit limits. Fit Failed"
                self.finish_pixel(i,j,failFlag=9)
                self.plot_pix_pdf(i,j,n_inbin,phase_bins,parameter_fit,redchi2,9)
                return False
            #elif np.sum(parameter_fit==parameter_lowerlimit) + np.sum(parameter_fit==parameter_upperlimit)-6 != 0.0:
            ##elif min(mpperr[9:])<=0.0 or min(mpperr[:6])<=0.0:
            #    if self.verbose:
            #        print "Blue, Red, or Noise peak hit fit limits. Fit Failed"
            #    self.finish_pixel(i,j,failFlag=9)
            #    self.plot_pix_pdf(i,j,n_inbin,phase_bins,parameter_fit,redchi2,9)
            #    continue
            #elif np.sum((parameter_fit==parameter_guess))-3 != 0.0:
            #    if self.verbose:
            #        print "Some parameters weren't varied. Fit Failed"
            #    self.finish_pixel(i,j,failFlag=10)
            #    self.plot_pix_pdf(i,j,n_inbin,phase_bins,parameter_fit,redchi2,10)
            elif redchi2==None:
                if self.verbose:
                    print "reduced chi^2 == None. Fit Failed"
                self.finish_pixel(i,j,failFlag=11)
                self.plot_pix_pdf(i,j,n_inbin,phase_bins,parameter_fit,redchi2,11)
                return False
            elif redchi2>self.params['max_chi2_all']:
                if self.verbose:
                    print "reduced chi^2 is too large: "+str(redchi2)+". Fit Failed"
                self.finish_pixel(i,j,failFlag=12)
                self.plot_pix_pdf(i,j,n_inbin,phase_bins,parameter_fit,redchi2,12)
                return False

        if not try_2_peak_fit and cut_off_phase > parameter_fit[7]+parameter_fit[6]:
            self.sensitivity.append(bluePeak_Fit[1])
            self.noisecutoff.append(cut_off_phase)
            self.deltaNoise.append(phase_bins[-1]-cut_off_phase)

        ## Now fit parabola to get wavelength <-> phase amplitude correspondance
        try:
            polyfit, solnrange, sigma = self.fitparabola(n_inbin,phase_bins,parameter_fit,make_plot=self.debug)
        except:
            if self.verbose:
                print "Failed while fitting wavelength cal to peaks!"
            self.finish_pixel(i,j,failFlag=13)
            self.plot_pix_pdf(i,j,n_inbin,phase_bins,parameter_fit,redchi2,13)
            if self.debug:
                raise
            else:
                return False

        ## Fit Worked!
        self.rescounter+=1
        self.finish_pixel(i,j,polyfit, sigma, solnrange, parameter_fit, mpperr, failFlag=0)
        self.plot_pix_pdf(i,j,n_inbin,phase_bins,parameter_fit,redchi2,0)     #Add pixel to PDF figure

        return True






    def makeDiagnositcPlots(self):

        mpl.rcParams['font.size'] = 10
        #blueSigArr=np.asarray(self.blueSig)
        #noiseSigArr=np.asarray(self.noiseSig)
        sensArr=np.asarray(self.sensitivity)
        cutoffArr=np.asarray(self.noisecutoff)
        deltanoiseArr=np.asarray(self.deltaNoise)
        
        plt.figure()
        plt.plot(sensArr,cutoffArr,'.b')
        plt.xlabel("sensitivity (location of blue peak)")
        plt.ylabel("phase cutoff")

        plt.figure()
        plt.plot(sensArr,deltanoiseArr,'.b')
        plt.xlabel("sensitivity (location of blue peak)")
        plt.ylabel("extended fall off")
        '''


        plt.figure()
        plt.plot(sensArr,blueSigArr,'.b',label="blue sigma")
        plt.plot(sensArr,noiseSigArr,'.r',label="noise sigma")
        plt.xlabel("sensitivity (location of blue peak)")
        plt.ylabel("Sigma")
        plt.legend()
        
        plt.figure()
        plt.plot(sensArr,blueSigArr/noiseSigArr,'.',label="blue/noise sigma ratio")
        plt.xlabel("sensitivity (location of blue peak)")
        plt.ylabel("Sigma ratio")
        plt.legend()

        plt.figure()
        plt.plot(blueSigArr,noiseSigArr,'.r',label="noise Sig")
        plt.plot(blueSigArr,blueSigArr,'-k',label="1")
        plt.xlabel("blue sigma")
        plt.ylabel("noise sigma")
        plt.legend()

        plt.figure()
        plt.plot(sensArr,sensArr/blueSigArr,'.b')
        plt.xlabel("sensitivity (location of blue peak)")
        plt.ylabel("R (location of blue peak / sigma blue peak)")
        '''

        plt.show()

    def plot_pix_pdf(self,iRow,iCol,n_inbin,phase_bins,fit_params,redchi2=None,failFlag=0):

        if self.debug:
            plt.show()
            return

        if not self.save_pdf:
            return

        if iRow==None or iCol==None:
            if self.plotcounter % self.n_plots_per_page!=0:
                try:
                    self.pp.savefig(self.fig_pdf)
                except:
                    pass
                #self.fig_pdf.clf()
            self.pp.close()
            return

        if (self.plotcounter % self.n_plots_per_page == 0):
            self.fig_pdf = plt.figure(figsize=(8.25, 10), dpi=100)

        plt.subplot(self.n_plots_y,self.n_plots_x,self.plotcounter%self.n_plots_per_page+1)
        plt.plot(phase_bins,n_inbin, label=str(failFlag))                    
        titlestring = '('+str(iRow)+', '+str(iCol)+') --> '+self.laserCalFile.beamImage[iRow][iCol]
        plt.title(titlestring)
        if fit_params!=None:
            plt.xlim(fit_params[1]-2.5*fit_params[0],phase_bins[-1])
        else:
            plt.xlim(-600., 0.0)
        #plt.ylim(0.,np.max(max_vals)+60.)


        model=self.params['model_type']
        if failFlag!=0 and failFlag <=6:
            model='gaussian'
        ##Model
        model_list = {
            'parabola': parabola,
            'gaussian': gaussian,
            'fourgaussian': fourgaussian,
            'threegaussian_exp': threegaussian_exp,
            'threegaussian_exppow': threegaussian_exppow,
            'threegaussian_moyal': threegaussian_moyal,
            'threegaussian_power': threegaussian_power,
            'threegaussian_lorentzian': threegaussian_lorentzian}
        
        if fit_params != None:
            modelArr=model_list[model](p=fit_params,x=phase_bins,return_models=True)
            fullModel=np.zeros(len(n_inbin))
            ylimit_list=[]
            for model_part in modelArr:
                fullModel = fullModel+model_part
                ylimit_list.append(max(model_part))

        if failFlag==0:
            plt.plot(phase_bins,fullModel,'r',label="$\chi^{2}_{red}$: "+str(redchi2))
            ylimit=max(ylimit_list[:-1])*2.0
            plt.ylim(0,ylimit)

        elif fit_params != None:
            plt.plot(phase_bins,fullModel,'g',label="$\chi^{2}_{red}$: "+str(redchi2))
            ylimit=max(ylimit_list)*2.0
            if len(fit_params)>3:
                for model_part in modelArr:
                    plt.plot(phase_bins,model_part,'k')
                ylimit=max(ylimit_list[:-1])*2.0
            plt.ylim(0,ylimit)

        a=plt.legend(loc=2)
        if failFlag==0:
            num_peaks_fit=3
            if fit_params[8]==0:
                num_peaks_fit=2
            txt=offsetbox.TextArea("Fit "+str(num_peaks_fit)+" peaks")
            box=a._legend_box
            box.get_children().append(txt)
            box.set_figure(box.figure)
            


        if ((self.plotcounter +1) % self.n_plots_per_page == 0):
            self.pp.savefig(self.fig_pdf)
            #self.fig_pdf.clf()

        self.plotcounter+=1
            

                


if __name__ == '__main__':
    
    try:
        paramFile = sys.argv[1]
    except IndexError:
        paramFile=os.getenv('ARCONS_PARAMS_PATH', default=os.path.expanduser('~')+'/ARCONS-pipeline/params/')+'waveCal.dict'
        #paramFile = '/home/abwalter/ARCONS-pipeline/params/waveCal.dict'
        print "Loading parameters from: "+paramFile
    calFNs, params = getCalFileNames(paramFile)
    #calFNs, params = getCalFileNames(paramFile,wavecal='obs_',timeMaskDir='/ChargeDrift')
    for calFN in calFNs:
        print calFN.cal()
    debug=False
    verbose=False
    for calFN in calFNs:
        waveCalObject=waveCal(calFN,params,verbose=verbose,debug=debug)
        try:
            waveCalObject.findWaveLengthSoln()
        except:
            if hasattr(waveCalObject,'pp'):
                waveCalObject.pp.close()
            if len(waveCalObject.pixRowArr) != 0:
                print "last completed resonator: ("+str(waveCalObject.pixRowArr[-1])+', '+str(waveCalObject.pixColArr[-1])+')'
            raise
        if not debug:
            if waveCalObject.rescounter > 0:
                waveCalObject.write_waveCal()
                waveCalObject.write_waveCal_drift()
                try:
                    diag_obj=waveCal_diagnostic(calFN,params,save=True)
                    diag_obj.make_R_array()
                    diag_obj.plot_R_array()
                    diag_obj.plot_nlaser_array()
                    diag_obj.plot_R_hist()
                    #diag_obj.plot_parabolafit_hist()
                    #diag_obj.plot_parameterfit_hist()
                    #diag_obj.plot_sigmas()
                    #diag_obj.plot_amps()
                    #plt.show()
                    plt.close('all')
                    del diag_obj
                except:
                    pass
            else:
                print "No solutions found!"
        del waveCalObject



