'''
Author: Alex Walter                             Date: March 18, 2014

This code analyzes individual wavecal solutions. You can use it by itself but usually it is automatically called by waveCal.py

Input:
    paramFile - same as wavecal parameter file
Output:
    - map of R at 400nm
    - map of number of lasers used in solution
    - histogram of R at 400 nm
    - histogram of various parameters
'''

import sys, os
import tables
from tables import *
import numpy as np

from util.readDict import readDict
from util.FileName import FileName
from fitFunctions import *
from utils import *
from waveCal import *


class waveCal_diagnostic():
    def __init__(self,calFN,params,save=True):
        self.calFN = calFN
        self.params=params

        self.outpath=None
        if save:
            self.outpath=self.calFN.calSoln().rsplit('/',1)[0]+os.sep
        ## open cal file
        print 'Opening file: '+self.calFN.calSoln()
        self.calsol_file = tables.openFile(self.calFN.calSoln(), mode="r")
        ## open drift file
        print 'Opening file: '+self.calFN.calDriftInfo()
        self.drift_file = tables.openFile(self.calFN.calDriftInfo(), mode="r")
        

    def __del__(self):
        """
        Closes any cal files that are open
        """
        try:
            self.calsol_file.close()
        except:
            pass
        try:
            self.drift_file.close()
        except:
            pass

    def make_R_array(self, laser='blue'):
        pixRow=self.calsol_file.root.wavecal.calsoln.cols.pixelrow[:]
        nRows=max(pixRow)+1
        pixCol=self.calsol_file.root.wavecal.calsoln.cols.pixelcol[:]
        nCols=max(pixCol)+1
        pixSigma=self.calsol_file.root.wavecal.calsoln.cols.sigma[:]
        pixRoach=self.calsol_file.root.wavecal.calsoln.cols.roach[:]
        pix_polyfit = self.calsol_file.root.wavecal.calsoln.cols.polyfit[:]
        drift_row = self.drift_file.root.params_drift.driftparams.cols.pixelrow[:]    
        drift_col = self.drift_file.root.params_drift.driftparams.cols.pixelcol[:]    
        drift_params = self.drift_file.root.params_drift.driftparams.cols.gaussparams[:]

        ## make array plot of R values    
        self.xyrarray = np.zeros((nRows, nCols))
        self.roacharray =  np.zeros((nRows, nCols))
        self.nlaserarray =  np.zeros((nRows, nCols))
        

        for k in range(len(pixSigma)):
            if pixSigma[k]>0:
                drift_ind = np.where((drift_row==pixRow[k]) * (drift_col==pixCol[k]))[0][0]
                peak_fit = drift_params[drift_ind]
                
                #Set up nlaserarray
                if peak_fit[8]>0:
                    self.nlaserarray[pixRow[k]][pixCol[k]]=3
                else:   #peak_fit[8]==0
                    self.nlaserarray[pixRow[k]][pixCol[k]]=2
                if peak_fit[8]<0:
                    print "shouldn't happen"
                
                if not ((laser=='IR' or laser=='ir') and self.nlaserarray[pixRow[k]][pixCol[k]]==2):
                    #Calculate R
                    if laser == 'blue':
                        wavelength = self.params['bluelambda']
                        peak_phase = peak_fit[1]
                        sigma = peak_fit[0]
                    elif laser=='red':
                        wavelength = self.params['redlambda']
                        peak_phase = peak_fit[1]+peak_fit[4]
                        sigma = peak_fit[3]
                    elif laser=='IR' or laser=='ir':
                        wavelength = self.params['irlambda']
                        peak_phase = peak_fit[1]+peak_fit[4]+peak_fit[7]
                        sigma = peak_fit[6]
                    else:
                        print str(laser)+" invalid. Only have 3 lasers implemented: 'blue', 'red', 'IR'"
                        raise ValueError
                    laser_energy=self.params['h'] * self.params['c'] / (wavelength * self.params['ang2m'])
                    
                    #Calculate sigma in energy space
                    laser_sig_high = laser_energy - (parabola(pix_polyfit[k],x=np.asarray([peak_phase+sigma]),return_models=True))[0][0]
                    laser_sig_low = (parabola(pix_polyfit[k],x=np.asarray([peak_phase-sigma]),return_models=True))[0][0] - laser_energy
                    if laser_sig_high>0 and laser_sig_low>0:
                        laser_sig = (laser_sig_high + laser_sig_low)/2.
                        
                        l_en = (parabola(pix_polyfit[k],x=np.asarray([peak_phase]),return_models=True))[0][0]
                        R=laser_energy/(self.params['fwhm2sig']*laser_sig)
                        #print laser_energy,l_en,pixSigma[k],laser_sig, sigma,R
                        if R<0:
                            print pixRow[k], pixCol[k]
                            print peak_fit
                            print pix_polyfit[k]
                            print laser_energy,l_en,pixSigma[k],laser_sig, sigma,R
                        
                        self.xyrarray[pixRow[k]][pixCol[k]]=laser_energy/(self.params['fwhm2sig']*laser_sig)
                
                

            #else:
            #    self.xyrarray[pixRow[k]][pixCol[k]]=0

            self.roacharray[pixRow[k]][pixCol[k]]=int(pixRoach[k])

    def plot_nlaser_array(self):
        if self.outpath==None:
            plotArray(self.nlaserarray, showMe=True, cbar=True,plotTitle='Number of Lasers for Fit')
        else:
            fname=self.outpath+self.params['figdir']+'/calsol_' + self.calFN.tstamp +'_nlaserPlot.png'
            plotArray(self.nlaserarray, showMe=False, cbar=True, plotFileName=fname,plotTitle='Number of Lasers for Fit')
            print '\tSaving n laser plot to: '+fname
            plt.close('all')

    def plot_R_array(self,laser='blue'):
        if laser=='blue':
            wavelength = str(self.params['bluelambda'])
        elif laser=='red':
            wavelength = str(self.params['redlambda'])
        elif laser=='IR' or laser=='ir':
            wavelength = str(self.params['irlambda'])
        else:
            print str(laser)+" invalid. Only have 3 lasers implemented: 'blue', 'red', 'IR'"
            raise ValueError
            
        if self.outpath==None:
            plotArray( self.xyrarray, showMe=True, cbar=True,plotTitle='Energy Resolution at '+wavelength+' Angstroms')
        else:
            fname=self.outpath+self.params['figdir']+'/calsol_' + self.calFN.tstamp +'_arrayPlot_'+laser+'.png'
            plotArray( self.xyrarray, showMe=False, cbar=True, plotFileName=fname,plotTitle='Energy Resolution at '+wavelength+' Angstroms')#,normMax=6.5)
            print '\tSaving R array to: '+fname
            plt.close('all')

    def plot_R_hist(self,laser='blue'):
        if laser=='blue':
            wavelength = str(self.params['bluelambda'])
        elif laser=='red':
            wavelength = str(self.params['redlambda'])
        elif laser=='IR' or laser=='ir':
            wavelength = str(self.params['irlambda'])
        else:
            print str(laser)+" invalid. Only have 3 lasers implemented: 'blue', 'red', 'IR'"
            raise ValueError
            
        plt.figure()
        plt.subplot(111)
        colormap = mpl.cm.gist_ncar
        #print np.unique(self.roacharray)
        plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, len(np.unique(self.roacharray)))])
        for k in np.unique(self.roacharray):
            roach_R_values = np.copy(self.xyrarray.flatten()[np.where(self.roacharray.flatten()==k)])
            n_res, resbins = np.histogram(roach_R_values, 11*5, range = (1,12))
            resbins=(resbins+(resbins[1]-resbins[0])/2.0)[:-1]
            plt.plot(resbins,n_res,label='roach '+str(int(k)))
        plt.legend(loc='upper left')
        plt.xlim(1.,9.)
        plt.xlabel('Energy Resolution at '+wavelength+' Angstroms')
        if self.outpath==None:
            plt.show()
        else:
            fname=self.outpath+self.params['figdir']+'/calsol_' + self.calFN.tstamp +'_R_Estimates_'+laser+'.png'
            print '\tSaving R histogram to: '+fname
            plt.savefig(fname)
            plt.close('all')

    def plot_parabolafit_hist(self):
        pix_polyfit = self.calsol_file.root.wavecal.calsoln.cols.polyfit[:]
        plt.figure()
        for k in range(len(pix_polyfit[0])):
            plt.subplot(1,len(pix_polyfit[0]),1+k)
            poly_arr = pix_polyfit[:,k]
            poly_arr = poly_arr[np.where(poly_arr!=-1.)]
            poly_num, poly_bins = np.histogram(poly_arr,200)
            poly_bins=(poly_bins+(poly_bins[1]-poly_bins[0])/2.0)[:-1]
            plt.plot(poly_bins,poly_num)
        if self.outpath==None:
            plt.show()
        else:
            pass
        #plt.close('all')
        
    def plot_noiseEnergy_hist(self):
        """
        Finds the energy of the 4th gaussian fit to the noise tail
        """
        
        pix_polyfit = self.calsol_file.root.wavecal.calsoln.cols.polyfit[:]
        pixRow=self.calsol_file.root.wavecal.calsoln.cols.pixelrow[:]
        pixCol=self.calsol_file.root.wavecal.calsoln.cols.pixelcol[:]
        pixSigma=self.calsol_file.root.wavecal.calsoln.cols.sigma[:]
        drift_params = self.drift_file.root.params_drift.driftparams.cols.gaussparams[:]
        drift_row = self.drift_file.root.params_drift.driftparams.cols.pixelrow[:]    
        drift_col = self.drift_file.root.params_drift.driftparams.cols.pixelcol[:]    
        
        noise_en = []
        noise_amp = []
        for k in range(len(pixSigma)):
            if pixSigma[k]>0:
                drift_ind = np.where((drift_row==pixRow[k]) * (drift_col==pixCol[k]))[0][0]
                
                print drift_ind
                noise_amp.append(drift_params[drift_ind,11])
                print drift_params[drift_ind,10]
                en = (parabola(pix_polyfit[k],x=np.asarray([drift_params[drift_ind,10]]),return_models=True))[0][0]
                print en
                noise_en.append(en)
                
        plt.figure()
        ax=plt.gca()
        n_en, bin_en = np.histogram(noise_en,100)
        bin_en = (bin_en+(bin_en[1]-bin_en[0])/2.0)[:-1]
        ax.plot(bin_en,n_en)
        ax.set_xlabel("Noise peak energy (eV)")
        ax.set_ylabel("# pixels")
        
        plt.figure()
        ax=plt.gca()
        n_amp, bin_amp = np.histogram(noise_amp,100)
        bin_amp = (bin_amp+(bin_amp[1]-bin_amp[0])/2.0)[:-1] 
        ax.plot(bin_amp,n_amp)
        ax.set_xlabel("Noise peak amp (photons/sec)")
        ax.set_ylabel("# pixels")
        
        plt.show()
        
        

    def plot_parameterfit_hist(self):
        """
        Assumes 3 gaussians and 3 parameter noise fit
        """
        drift_params = self.drift_file.root.params_drift.driftparams.cols.gaussparams[:]
        IR_amp = drift_params[:,8]
        plt.figure()
        for k in range(len(drift_params[0])):
            axL=plt.subplot(4,3,k+1)
            par_fit = drift_params[:,k]
            if k==11:
                #par_fit = par_fit[np.where(par_fit<10.**-12.)]
                #IR_amp = IR_amp[np.where(par_fit<10.**-12.)]
                par_fit = np.log(par_fit)
            par_fit_3 = par_fit[np.where(IR_amp>0)]
            par_fit_2 = par_fit[np.where(IR_amp==0)]
            if k in [6,7,8]:
                par_fit_2 = []
            #print len(drift_params[:,0])
            if k==9:
                par_num_3, par_bins_3 = np.histogram(par_fit_3,100,range=(2,3))
            else:
                par_num_3, par_bins_3 = np.histogram(par_fit_3,100)
            par_bins_3=(par_bins_3+(par_bins_3[1]-par_bins_3[0])/2.0)[:-1]
            axL.plot(par_bins_3,par_num_3,'b')
            #plt.ylabel('# pixels with 3 peaks')
            #plt.xlabel('Parameter Value')

            axR=axL.twinx()
            #axR.yaxis.tick_right()
            #axR.yaxis.set_label_position("right")
            if k==9:
                par_num_2, par_bins_2 = np.histogram(par_fit_2,100,range=(2,3))
            else:
                par_num_2, par_bins_2 = np.histogram(par_fit_2,100)
            par_bins_2=(par_bins_2+(par_bins_2[1]-par_bins_2[0])/2.0)[:-1]
            axR.plot(par_bins_2,par_num_2,'r')
            #plt.ylabel('# pixels with 2 peaks')

        if self.outpath==None:
            plt.show()
        else:
            pass
        #plt.close('all')

    def plot_sigmas(self):
        drift_params = self.drift_file.root.params_drift.driftparams.cols.gaussparams[:]
        blue_sigma = drift_params[:,0]
        red_sigma = drift_params[:,3]
        IR_sigma = drift_params[:,6]
        IR_amp = drift_params[:,8]
        #IR_sigma=IR_sigma[np.where(IR_amp>0)]
        plt.figure()
        plt.plot(blue_sigma,blue_sigma,'-b',label='blue')
        plt.plot(blue_sigma[np.where(IR_amp>0)],red_sigma[np.where(IR_amp>0)],'or',label='red 3')
        avg=np.mean(red_sigma[np.where(IR_amp>0)]/blue_sigma[np.where(IR_amp>0)])
        std=np.std(red_sigma[np.where(IR_amp>0)]/blue_sigma[np.where(IR_amp>0)])
        print "Average red/blue sigma_3 ratio: "+str(avg)+' +/- '+str(std)
        plt.plot(blue_sigma[np.where(IR_amp==0)],red_sigma[np.where(IR_amp==0)],'om',label='red 2')
        avg=np.mean(red_sigma[np.where(IR_amp==0)]/blue_sigma[np.where(IR_amp==0)])
        std=np.std(red_sigma[np.where(IR_amp==0)]/blue_sigma[np.where(IR_amp==0)])
        print "Average red/blue sigma_2 ratio: "+str(avg)+' +/- '+str(std)
        plt.plot(blue_sigma[np.where(IR_amp>0)],IR_sigma[np.where(IR_amp>0)],'ok',label='IR')
        avg=np.mean(IR_sigma[np.where(IR_amp>0)]/blue_sigma[np.where(IR_amp>0)])
        std=np.std(IR_sigma[np.where(IR_amp>0)]/blue_sigma[np.where(IR_amp>0)])
        print "Average IR/blue sigma ratio: "+str(avg)+' +/- '+str(std)
        plt.legend()
        plt.xlabel("Blue sigma")
        plt.ylabel("Laser sigma")
        if self.outpath==None:
            plt.show()
        else:
            pass
        #plt.close('all')

    def plot_amps(self):
        drift_params = self.drift_file.root.params_drift.driftparams.cols.gaussparams[:]
        blue_amp = drift_params[:,2]
        red_amp = drift_params[:,5]
        IR_amp = drift_params[:,8]
        #IR_amp = drift_params[:,8]
        #IR_sigma=IR_sigma[np.where(IR_amp>0)]
        plt.figure()
        plt.plot(blue_amp,blue_amp,'-b',label='blue')
        plt.plot(blue_amp[np.where(IR_amp>0)],red_amp[np.where(IR_amp>0)],'or',label='red 3')
        avg=np.mean(red_amp[np.where(IR_amp>0)]/blue_amp[np.where(IR_amp>0)])
        std=np.std(red_amp[np.where(IR_amp>0)]/blue_amp[np.where(IR_amp>0)])
        print "Average red/blue amp_3 ratio: "+str(avg)+' +/- '+str(std)
        plt.plot(blue_amp[np.where(IR_amp==0)],red_amp[np.where(IR_amp==0)],'om',label='red 2')
        avg=np.mean(red_amp[np.where(IR_amp==0)]/blue_amp[np.where(IR_amp==0)])
        std=np.std(red_amp[np.where(IR_amp==0)]/blue_amp[np.where(IR_amp==0)])
        print "Average red/blue amp_2 ratio: "+str(avg)+' +/- '+str(std)
        plt.plot(blue_amp[np.where(IR_amp>0)],IR_amp[np.where(IR_amp>0)],'ok',label='IR')
        avg=np.mean(IR_amp[np.where(IR_amp>0)]/blue_amp[np.where(IR_amp>0)])
        std=np.std(IR_amp[np.where(IR_amp>0)]/blue_amp[np.where(IR_amp>0)])
        print "Average IR/blue amp ratio: "+str(avg)+' +/- '+str(std)
        plt.legend()
        plt.xlabel("Blue amplitude")
        plt.ylabel("Laser amplitude")
        if self.outpath==None:
            plt.show()
        else:
            pass
        #plt.close('all')


if __name__ == '__main__':
    try:
        paramFile = sys.argv[1]
    except IndexError:
        paramFile=os.getenv('PYTHONPATH',default=os.path.expanduser('~')+'/ARCONS-pipeline').split(':')[0]+'/params/waveCal.dict'
        #paramFile = '/home/abwalter/ARCONS-pipeline/params/waveCal.dict'
        print "Loading parameters from: "+paramFile
    calFNs, params = getCalFileNames(paramFile)
    #calFNs, params = getCalFileNames(paramFile,wavecal='obs_',timeMaskDir='/ChargeDrift')

    for calFN in calFNs:
        try:
            diag_obj=waveCal_diagnostic(calFN,params,save=False)
            diag_obj.plot_noiseEnergy_hist()
            #laser='blue'
            #diag_obj.make_R_array(laser=laser)
            #diag_obj.plot_R_array(laser=laser)
            #diag_obj.plot_nlaser_array()
            #diag_obj.plot_R_hist(laser=laser)
            #diag_obj.plot_parameterfit_hist()
            #diag_obj.plot_parameterfit_hist()
            #diag_obj.plot_sigmas()
            #diag_obj.plot_amps()
            plt.show()
            plt.close('all')
            del diag_obj
        except:
            raise
            pass




