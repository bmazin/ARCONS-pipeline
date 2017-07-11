'''
Author: Alex Walter
Date: Sept 14, 2014

This code is designed to visualize QE data combined with wavecal data.
Run QECalibration.py first.
'''

import sys,os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from multiprocessing import Process
from popup import *
from util.ObsFile import ObsFile
from util.readDict import readDict
from util.FileName import FileName
from fitFunctions import *
from waveCal_diagnostics import *
from waveCal import fitData

def f(self,phase_bins=[0],n_inbin=[0],Peak_Fit=[1,0,0]):
    self.axes.plot(phase_bins,n_inbin,'.')
    
    xArr=np.arange(phase_bins[0],phase_bins[-1],(phase_bins[1]-phase_bins[0])/10.)
    modelArr=(gaussian(p=Peak_Fit,x=xArr,return_models=True))[0]
    self.axes.plot(xArr,modelArr,'r-')
    
def plot_QE_phaseHeights(self,list_of_fits,list_of_phase,list_of_ninbin):
    cmap=matplotlib.cm.jet
    wavelengths = self.parent.QE_data[:,0]
    for k in range(len(wavelengths)):
        color=cmap((k+1.)/len(wavelengths))
        phase_bins = list_of_phase[k]
        n_inbin = list_of_ninbin[k]
        self.axes.plot(phase_bins,n_inbin,'.',color=color)
        
    for k in range(len(wavelengths)):
        color=cmap((k+1.)/len(wavelengths))
        Peak_Fit = list_of_fits[k]
        if not (Peak_Fit[0]==1 and Peak_Fit[1]==0 and Peak_Fit[2]==0):
            phase_bins = list_of_phase[k]
            xArr=np.arange(phase_bins[0],phase_bins[-1],(phase_bins[1]-phase_bins[0])/10.)
            modelArr=(gaussian(p=Peak_Fit,x=xArr,return_models=True))[0]
            self.axes.plot(xArr,modelArr,'-',color=color)

def guessFirstPeak(self,n_inbin_total,phase_bins_total):
    """
        Guess location of first peak in phase histogram. First finds where channelizer trigger suddenly catches photons and 
        estimates the energy resolution from the size of the phase noise. Next, smoothes data and finds location of peaks and valleys. 
        The first peak surrounded by two valleys is evidently the blue laser peak. 
        
        Inputs:
            n_inbin_total - total number of photons counted in each phase bin
            phase_bins_total - A list of phase values making up the phase bins. ADC/DAC units
            
        Outputs:
            bluePeak_guess - parameters for blue peak gaussian
            bluePeak_lowerlimit - lower limit of parameters
            bluePeak_upperlimit - upper limit of parameters
            phase_cut_off - phase at edge of blue peak. This allows us to ignore the red/IR/noise peaks while fitting the blue peak
            
        Uses parameters from params dictionary:
            min_amp - Cuts data with fewer than min_amp counts at the beginning and end of the phase height histogram
            threshold_sigma - The number of sigmas used for photon triggering in channelizer. 
                              If the energy resolution is determined by the phase noise we can guess the sigma of the blue peak
            bin_smooth - number of bins to use in the smoothing window
    """
    params = self.params
    params['min_amp'] = 2
    params['bin_smooth'] = 60
    last_ind_temp = np.where(n_inbin_total>params['min_amp'])[0][-1]
    sigma4_guess=np.abs(phase_bins_total[last_ind_temp]/(params['threshold_sigma']))
    #print sigma4_guess
    window = 'hanning'
    windowsize = int(params['bin_smooth'])   

    parab_smooth = smooth.smooth(n_inbin_total, windowsize, window)
    smoothed_data = np.array(parab_smooth, dtype=float)
    first_ind=np.where(smoothed_data>params['min_amp'])[0][0]
    #last_ind = np.where(n_inbin_total>self.params['min_amp'])[0][-1]
    last_ind = np.argmax(n_inbin_total)
    if first_ind+50>=last_ind:
        raise RuntimeError("No Blue Peak, try lowering 'min_amp'")

    

    #pop = pixel_PopUp(parent=self,plot_func=f,phase_bins=phase_bins_total,n_inbin=n_inbin_total,smoothed_data=smoothed_data)
    #pop.show()
    
    #smoothed_data=smoothed_data[first_ind:last_ind]
    #phase_bins=phase_bins_total[first_ind:last_ind]
    #n_inbin=n_inbin_total[first_ind:last_ind]
    
    phase_bins=phase_bins_total
    n_inbin=n_inbin_total
    

    
    #return 0,0,0,0
    
    gradients=np.diff(smoothed_data)>0
    min_Ind=[]
    max_Ind=[]
    for k in range(len(gradients)-1):
        if gradients[k]==True and gradients[k+1]==False:
            if smoothed_data[k] >= 1.4:
                max_Ind.append(k)
        if gradients[k]==False and gradients[k+1]==True and len(max_Ind)>0:
            min_Ind.append(k)

    #while max_Ind[0] <= min_Ind[0]:
    #    max_Ind = max_Ind[1:]
    #print phase_bins[min_Ind],phase_bins[max_Ind]
    try:
        #print phase_bins[max_Ind]
        #print phase_bins[min_Ind]
        blueMaxPhase_guess=phase_bins[max_Ind[0]]
        smallest_min_after_first_max = min_Ind[np.argmin(smoothed_data[min_Ind])]
        blueMinPhase_guess=phase_bins[smallest_min_after_first_max]
        blueMaxAmp_guess=smoothed_data[max_Ind[0]]
    except:
        #print "Can't find Blue Peak local maxima"
        #print "\tmaximums: "+str(phase_bins[max_Ind])
        #print "\tminimums: "+str(phase_bins[min_Ind])
        raise
    
    #print blueMinPhase_guess,blueMaxPhase_guess,blueMaxAmp_guess
    
    #blueMaxAmp_guess=smoothed_data[max_Ind[0]]
    #blueSigma_guess=np.abs(blueMaxPhase_guess/params['blueR_guess'])/2.0
    blueSigma_guess = sigma4_guess

    bluePeak_guess=[blueSigma_guess,blueMaxPhase_guess,blueMaxAmp_guess]
    bluePeak_lowerlimit=[blueSigma_guess/10.,blueMaxPhase_guess-np.abs(blueMaxPhase_guess-blueMinPhase_guess),2.0]
    bluePeak_upperlimit=[blueSigma_guess*10.,blueMaxPhase_guess+np.abs(blueMaxPhase_guess-blueMinPhase_guess),2.*blueMaxAmp_guess]

    return bluePeak_guess, bluePeak_lowerlimit, bluePeak_upperlimit, phase_bins[smallest_min_after_first_max]

def fit_QE_gaussians(self,row,col):

    list_of_fits = []
    list_of_phase = []
    list_of_ninbin = []

    for wave_ind in range(len(QE_param)):
        
        startTime = QE_param[wave_ind,2]
        endTime = QE_param[wave_ind,3]
        #print QE_data[wave_ind,0],startTime,endTime
        #continue
        duration = endTime-startTime
        
        dataDict=self.QE_obs.getTimedPacketList(row,col,firstSec=startTime,integrationTime=duration,timeSpacingCut=self.params['danicas_cut'])
        peakHeights=np.asarray(dataDict['peakHeights'])*1.0
        ## less than 'min_amp' per second average count rate            
        if dataDict['effIntTime']==0.0 or len(peakHeights)<=(dataDict['effIntTime']*self.params['min_count_rate']):
            #print 'Not enough data for phase histogram'
            list_of_fits.append([1,0,0])
            list_of_phase.append([0,1])
            list_of_ninbin.append([0,0])
            continue
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
        list_of_phase.append(phase_bins)
        list_of_ninbin.append(n_inbin)
        try:
            parameter_guess,parameter_lowerlimit,parameter_upperlimit, cut_off_phase = guessFirstPeak(self,n_inbin,phase_bins)
            #print cut_off_phase
        except:
            list_of_fits.append([1,0,0])
            #print 0
            continue
        #print parameter_guess
        #continue
        Peak_Fit, redchi2, mpperr = fitData(phase_bins,n_inbin,parameter_guess,parameter_lowerlimit,parameter_upperlimit,model='gaussian',cut_off_phase=cut_off_phase)
        #print Peak_Fit
        
        #Check that we didn't just fit the noise tail
        if Peak_Fit[1]+50 > phase_bins[-1]:
            #print '('+str(col)+', '+str(row)+') @ '+str(QE_data[wave_ind,0]) +' fit noise: '+str(Peak_Fit)
            #print 'cut_off_phase: '+str(cut_off_phase)
            win_title='('+str(col)+', '+str(row)+') @ '+str(QE_data[wave_ind,0]) +' fitted noise: '+str(Peak_Fit)
            pop = pixel_PopUp(parent=self,plot_func=f,win_title=win_title,phase_bins=phase_bins,n_inbin=n_inbin,Peak_Fit=Peak_Fit)
            pop.show()
            list_of_fits.append([1,0,0])
            continue
        
        list_of_fits.append(np.copy(Peak_Fit))
        
    for wave_ind in range(len(QE_param)):
        pop = pixel_PopUp(parent=self,plot_func=f,phase_bins=list_of_phase[wave_ind],n_inbin=list_of_ninbin[wave_ind],Peak_Fit=list_of_fits[wave_ind])
        pop.show()
    
    
    return np.asarray(list_of_fits), list_of_phase, list_of_ninbin
    
def plot_pixel_soln(self,row,col,list_of_fits):
    #self.fig, axarr = plt.subplots(2, sharex=True)
    #self.axes=axarr[0]
    #axes2=axarr[1]
    self.fig.delaxes(self.axes)
    gs = gridspec.GridSpec(2,1,height_ratios=[4,1])
    self.fig.subplots_adjust(hspace=0.,wspace=0.)
    self.axes=self.fig.add_subplot(gs[0])

    QE_wavelengths = self.parent.QE_data[:,0]
    QE_energies = [self.parent.params['h'] * self.parent.params['c'] / (x*10. * self.parent.params['ang2m']) for x in QE_wavelengths]
    QE_x_offsets = np.asarray(list_of_fits)[:,1]
    
    waveFN = FileName(obsFile=self.parent.wave_obs)
    calFile=tables.openFile(waveFN.calSoln(),mode='r')
    pixRow=calFile.root.wavecal.calsoln.cols.pixelrow[:]
    pixCol=calFile.root.wavecal.calsoln.cols.pixelcol[:]
    cal_ind = np.where((pixRow==row) * (pixCol==col))[0][0]
    polyfit=calFile.root.wavecal.calsoln.cols.polyfit[cal_ind]
    
    if np.all(polyfit==-1) and np.all(QE_x_offsets==0):
        raise ValueError
    
    calFile.close()
    driftFile=tables.openFile(waveFN.calDriftInfo(),mode='r')
    drift_row = driftFile.root.params_drift.driftparams.cols.pixelrow[:]    
    drift_col = driftFile.root.params_drift.driftparams.cols.pixelcol[:]  
    try:  
        drift_ind=(np.where(np.multiply(drift_row==row, drift_col==col)))[0][0]
        fit_params = driftFile.root.params_drift.driftparams.cols.gaussparams[drift_ind]
        cal_x_offsets = fit_params[[1,4,7]]
        cal_x_offsets[1]+=cal_x_offsets[0]
        cal_x_offsets[2]+=cal_x_offsets[1]
        if fit_params[8]==0.:
            cal_x_offsets[2]=0.
        cal_wavelengths = [self.parent.params['bluelambda'],self.parent.params['redlambda'], self.parent.params['irlambda']]     #Angstroms
        cal_energies = [self.parent.params['h'] * self.parent.params['c'] / (x * self.parent.params['ang2m']) for x in cal_wavelengths]             #eV
    except IndexError:
        cal_wavelengths = [self.parent.params['bluelambda'],self.parent.params['redlambda'], self.parent.params['irlambda']]     #Angstroms
        cal_energies = [self.parent.params['h'] * self.parent.params['c'] / (x * self.parent.params['ang2m']) for x in cal_wavelengths]             #eV
        cal_x_offsets = [0,0,0]
    driftFile.close()
    

    
    xArr = np.arange(-5000.,0.,0.1)
    if not np.all(polyfit==-1):
        modelArr=(parabola(p=polyfit,x=xArr,return_models=True))[0]
        self.axes.plot(xArr,modelArr,'k-')
    
    self.axes.plot(cal_x_offsets,cal_energies,'ko', label='Laser Cal')
    
    self.axes.plot(QE_x_offsets,QE_energies,'gx',label='QE wavelength fits')
    

    
    #Do linear fit of QE data to see if x-intercept is at 0
    QE_x_offsets = np.asarray(QE_x_offsets)
    QE_energies = np.asarray(QE_energies)
    energies = (QE_energies[np.where(QE_x_offsets<0)])
    sort_ind = np.argsort(energies)
    energies = energies[sort_ind]
    x_offsets = QE_x_offsets[np.where(QE_x_offsets<0)]
    x_offsets = x_offsets[sort_ind]
    if len(x_offsets)>1:
        parameter_guess = [0.0,(energies[0]-energies[-1])*1.0/(x_offsets[0]-x_offsets[-1]), 0.]
        parameter_lowerlimit=[None,None,0.]
        parameter_upperlimit=[None,None,0.]
        parameter_fit, redchi2gauss2, mpperr = fitData(x_offsets,energies,parameter_guess,parameter_lowerlimit,parameter_upperlimit,model='parabola')
        print 'x_intercept: '+str(-1.*parameter_fit[0]/parameter_fit[1])
    else:
        print 'x_intercept: Not enough data'
    #lin_fit=(parabola(p=parameter_fit,x=xArr,return_models=True))[0]
    #self.axes.plot(xArr,lin_fit,'g--')
    
    #Do parabola fit to QE data
    if len(x_offsets)>1:
        parameter_guess = [0.0,(energies[0]-energies[-1])*1.0/(x_offsets[0]-x_offsets[-1]), 0.]
        parameter_lowerlimit=[None,None,None]
        parameter_upperlimit=[None,None,None]
        parameter_fit, redchi2gauss2, mpperr = fitData(x_offsets,energies,parameter_guess,parameter_lowerlimit,parameter_upperlimit,model='parabola')
        print polyfit,' --> ',parameter_fit
        parab_fit=(parabola(p=parameter_fit,x=xArr,return_models=True))[0]
        self.axes.plot(xArr,parab_fit,'g--')
    
    #Plot Residuals
    if not np.all(polyfit==-1):
        cal_QE=(parabola(p=polyfit,x=x_offsets,return_models=True))[0]
        axes2=self.fig.add_subplot(gs[1],sharex=self.axes)
        #axes2 = self.fig.add_axes((.1,.1,.8,.2))
        axes2.plot(x_offsets,cal_QE-energies,'.')
        axes2.axhline(0,color='black')
    
    
    #Add labels etc
    self.axes.set_ylabel('Energy (eV)')
    self.axes.set_xlim(min(min(cal_x_offsets),min(QE_x_offsets))-200.,0.)
    self.axes.set_ylim(0.,5.)
    if not np.all(polyfit==-1):
        axes2.set_xlabel('Phase (ADC/DAC units)')
        axes2.set_ylabel('Calibration residuals (eV)')
        axes2.set_ylim(-3.,3.)
        #self.axes.set_xticklabels(self.axes.get_xticklabels(),visible=False)
        plt.setp(self.axes.get_xticklabels(),visible=False)

def plot_pixel_phaseHeights(self,row=None,col=None):

    self.row = row
    self.col = col
    self.laserCalFile = self.parent.wave_obs
    
    dataDict=self.laserCalFile.getTimedPacketList(self.row,self.col,timeSpacingCut=self.parent.params['danicas_cut'])
    peakHeights=np.asarray(dataDict['peakHeights'])*1.0
    ## less than 'min_amp' per second average count rate

    if dataDict['effIntTime']==0.0 or len(peakHeights)<=(dataDict['effIntTime']*self.parent.params['min_count_rate']):
        print 'Not enough data for phase histogram'
        raise ValueError
    baselines=np.asarray(dataDict['baselines'])*1.0
    peakHeights-=baselines
    biggest_photon = int(min(peakHeights))
    n_inbin,phase_bins=np.histogram(peakHeights,bins=np.abs(biggest_photon),range=(biggest_photon,0))
    phase_bins=(phase_bins+(phase_bins[1]-phase_bins[0])/2.0)[:-1]
    self.n_inbin = n_inbin
    self.phase_bins = phase_bins
    try:
        last_ind = np.where(n_inbin>self.parent.params['min_amp'])[0][-1]
    except IndexError:
        last_ind=len(n_inbin)-1
    ## Cut out all the zeros on the right


    self.axes.plot(phase_bins,n_inbin, 'b.',label="data")
    self.axes.set_xlim(phase_bins[(np.where(n_inbin >= 3))[0][0]],phase_bins[last_ind])

    n_inbin = n_inbin[:last_ind]
    phase_bins = phase_bins[:last_ind]

    #Check if this pixel is in drift file
    waveFN = FileName(obsFile=self.laserCalFile)
    driftFile=tables.openFile(waveFN.calDriftInfo(),mode='r')
    drift_row = driftFile.root.params_drift.driftparams.cols.pixelrow[:]    
    drift_col = driftFile.root.params_drift.driftparams.cols.pixelcol[:]    

    args=(np.where(np.multiply(drift_row==self.row, drift_col==self.col)))[0]
    if len(args==1):
        #pixel has a solution!
        fit_params = driftFile.root.params_drift.driftparams.cols.gaussparams[args[0]]
        #print fit_params
        model = driftFile.root.params_drift.driftparams.attrs.model_type
        model_list = {
            'parabola': parabola,
            'gaussian': gaussian,
            'fourgaussian': fourgaussian,
            'threegaussian_exp': threegaussian_exp,
            'threegaussian_exppow': threegaussian_exppow,
            'threegaussian_moyal': threegaussian_moyal,
            'threegaussian_power': threegaussian_power,
            'threegaussian_lorentzian': threegaussian_lorentzian}

        modelArr=model_list[model](p=fit_params,x=phase_bins,return_models=True)
        fullModel=np.zeros(len(n_inbin))
        for model_part in modelArr:
            fullModel = fullModel+model_part
            self.axes.plot(phase_bins,model_part, 'k--')
        self.axes.plot(phase_bins,fullModel, 'r-',label="model")

    self.axes.legend()
    self.axes.set_xlabel("phase")
    self.axes.set_ylabel("counts")
    self.axes.set_title('('+str(self.col)+', '+str(self.row)+') --> '+(self.laserCalFile.beamImage[self.row,self.col]).rsplit('/',1)[0])
    #res_blue = (self.drift.r_blue[self.row,self.col,self.closest_time_ind])
    #self.axes.text(0.15,0.95,"R ~ "+str(round(res_blue,3)),verticalalignment='center', horizontalalignment='center', transform=self.axes.transAxes)
    #print res_blue
    driftFile.close()

    return True



class waveCal_QE_PopUp(array_PopUp):
    def __init__(self,params,waveObs,QE_obs,QE_data,QE_param,parent=None):

        self.params=params
        self.wave_obs = waveObs
        self.QE_obs = QE_obs
        self.QE_data = QE_data
        self.QE_param = QE_param

        waveFN = FileName(obsFile=self.wave_obs)
        wvlcal_diag = waveCal_diagnostic(waveFN,params,save=False)
        wvlcal_diag.make_R_array()
        
        win_title='waveCal QE Popup'
        super(waveCal_QE_PopUp,self).__init__(wvlcal_diag.xyrarray,parent,win_title,title='Energy Resolution')


    def clickCanvas(self,event):
        if event.inaxes is self.axes:
            col = int(round(event.xdata))
            row = int(round(event.ydata))
            
            try:
                win_title='WaveCal Phase Heights for ('+str(col)+', '+str(row)+') --> '+(self.wave_obs.beamImage[row,col]).rsplit('/',1)[0]
                pop=pixel_PopUp(parent=self,win_title=win_title,plot_func=plot_pixel_phaseHeights,col=col,row=row)
                pop.show()
            except ValueError:
                #print '1'
                pass
            
            list_of_fits, list_of_phase, list_of_ninbin = fit_QE_gaussians(self,row=row,col=col)
            
            if not np.all(np.hstack(np.asarray(list_of_ninbin).flat)==0):
                win_title = 'QE Phase Heights for ('+str(col)+', '+str(row)+') --> '+(self.wave_obs.beamImage[row,col]).rsplit('/',1)[0]
                pop=pixel_PopUp(parent=self,win_title=win_title,plot_func=plot_QE_phaseHeights,list_of_fits=list_of_fits,list_of_phase=list_of_phase,list_of_ninbin=list_of_ninbin)
                pop.show()
            else:
                print "No QE data to plot"
            
            try:
                win_title = 'WaveCal Solution for ('+str(col)+', '+str(row)+') --> '+(self.wave_obs.beamImage[row,col]).rsplit('/',1)[0]
                pop=pixel_PopUp(parent=self,win_title=win_title,plot_func=plot_pixel_soln,col=col,row=row,list_of_fits=list_of_fits)
                pop.show()
            except ValueError:
                print "No QE nor waveCal data to plot"
                #print 'Unable to make plots for ('+str(col)+', '+str(row)+') --> '+(self.wave_obs.beamImage[row,col]).rsplit('/',1)[0]
                #raise
                pass

def pop_waveCal_QE(*args,**kwargs):
    #Waring: Does not play well with matplotlib state machine style plotting!
    block = kwargs.pop('block',False)
    def f(*args,**kwargs):
        form = waveCal_QE_PopUp(*args,**kwargs)
        #form = test_PopUp(*args,**kwargs)
        form.show()
    if block==True:
        f(*args,**kwargs)
        return None
    else:
        proc = Process(target=f,args=args,kwargs=kwargs)
        proc.start()
        return proc


if __name__ == "__main__":
    try:
        paramFile = sys.argv[1]
    except IndexError:
        paramFile=os.getenv('PYTHONPATH',default=os.path.expanduser('~')+'/ARCONS-pipeline').split(':')[0]+'/params/waveCal.dict'
        print "Loading parameters from: "+paramFile
    try:
        params = readDict(paramFile)
        params.readFromFile(paramFile)
    except:
        params = paramFile

    


    #waveObs = ObsFile('/Scratch/LABTEST/20140911/cal_20140911-190013.h5')
    #waveObs.loadBeammapFile('/Scratch/LABTEST/20140910/beammap_SCI6_B140731-Force_20140910.h5')
    #QE_obs = ObsFile('/Scratch/QEData/20140911/obs_20140911-194430.h5')
    #QE_obs.loadBeammapFile('/Scratch/LABTEST/20140910/beammap_SCI6_B140731-Force_20140910.h5')
    #QE_data = np.loadtxt('/Scratch/QEData/20140911/QE_20140911-194430.txt')
    #QE_param = np.loadtxt('/Scratch/QEData/20140911/QE_20140911-194430.param')
    waveObs = ObsFile('/Scratch/LABTEST/BOBAFETT/20140915/cal_20140915-231457.h5')
    #waveObs.loadBeammapFile('/Scratch/LABTEST/20140910/beammap_SCI6_B140731-Force_20140910.h5')
    QE_obs = ObsFile('/Scratch/QEData/20140915/obs_20140915-221013.h5')
    #QE_obs.loadBeammapFile('/Scratch/LABTEST/20140910/beammap_SCI6_B140731-Force_20140910.h5')
    QE_data = np.loadtxt('/Scratch/QEData/20140915/QE_20140915-221013.txt')
    QE_param = np.loadtxt('/Scratch/QEData/20140915/QE_20140915-221013.param')
    
    pop_waveCal_QE(params,waveObs,QE_obs,QE_data,QE_param)
    

