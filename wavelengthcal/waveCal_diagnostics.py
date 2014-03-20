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



def getCalFileNames(paramFile):
    params = readDict(paramFile)
    params.readFromFile(paramFile)

    #Load in cal files from param file
    run = params['run']
    sunsetDate = params['sunsetDate']
    tStampList = params['calTimeStamps']
    if tStampList == None:
        #search outpath for files
        intermDir=params['intermdir']
        outdir=params['outdir']
        if intermDir is None or intermDir is '':
            intermDir = os.getenv('INTERM_DIR', default="/Scratch")+os.sep
        if outdir is None or outdir is '':
            outdir = 'waveCalSolnFiles/'
        path=intermDir+outdir+run+os.sep+sunsetDate+os.sep
        tStampList = [f.split('.')[0].split('_')[1] for f in os.listdir(path) if f.endswith('.h5') and f.startswith('calsol_')]
    calFNs = [FileName(run=run, date=sunsetDate,tstamp=tStamp) for tStamp in tStampList]
    return calFNs, params


class waveCal_diagnostic():
    def __init__(self,calFN,params,save=True):
        self.calFN = calFN
        self.params=params
        intermDir=params['intermdir']
        outdir=params['outdir']
        if intermDir is None or intermDir is '':
            intermDir = os.getenv('INTERM_DIR', default="/Scratch")+os.sep
        if outdir is None or outdir is '':
            outdir = 'waveCalSolnFiles/'
        path=intermDir+outdir+calFN.run+os.sep+calFN.date+os.sep
        self.outpath=None
        if save:
            self.outpath=path
        ## open cal file
        cal_file_name = path+"calsol_" + calFN.tstamp + '.h5'
        print 'Opening file: '+cal_file_name
        self.calsol_file = tables.openFile(cal_file_name, mode="r")
        ## open drift file
        drift_file_name = path+params['driftdir']+'calsol_' + calFN.tstamp + '_drift.h5'
        print 'Opening file: '+drift_file_name
        self.drift_file = tables.openFile(drift_file_name, mode="r")

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

    def make_R_array(self):
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
        #blue_energy=params['h'] * params['c'] / (params['bluelambda'] * params['ang2m'])
        for k in range(len(pixSigma)):
            if pixSigma[k]>0:
                drift_ind = np.where((drift_row==pixRow[k]) * (drift_col==pixCol[k]))[0][0]
                peak_fit = drift_params[drift_ind]
                blue_energy = (parabola(pix_polyfit[k],x=np.asarray([peak_fit[1]]),return_models=True))[0][0]
                self.xyrarray[pixRow[k]][pixCol[k]]=blue_energy/(self.params['fwhm2sig']*pixSigma[k])
                if peak_fit[8]>0:
                    self.nlaserarray[pixRow[k]][pixCol[k]]=3
                else:
                    self.nlaserarray[pixRow[k]][pixCol[k]]=2
                if peak_fit[8]<0:
                    print "shouldn't happen"
            #else:
            #    self.xyrarray[pixRow[k]][pixCol[k]]=0

            self.roacharray[pixRow[k]][pixCol[k]]=int(pixRoach[k])

    def plot_nlaser_array(self):
        if self.outpath==None:
            plotArray(self.nlaserarray, showMe=True, cbar=True,plotTitle='Number of Lasers for Fit')
        else:
            fname=self.outpath+self.params['figdir']+'calsol_' + self.calFN.tstamp +'_nlaserPlot.png'
            plotArray(self.nlaserarray, showMe=False, cbar=True, plotFileName=fname,plotTitle='Number of Lasers for Fit')
            print '\tSaving n laser plot to: '+fname
            plt.close('all')

    def plot_R_array(self):
        if self.outpath==None:
            plotArray( self.xyrarray, showMe=True, cbar=True,plotTitle='Energy Resolution at 400nm')
        else:
            fname=self.outpath+self.params['figdir']+'calsol_' + self.calFN.tstamp +'_arrayPlot.png'
            plotArray( self.xyrarray, showMe=False, cbar=True, plotFileName=fname,plotTitle='Energy Resolution at 400nm')
            print '\tSaving R array to: '+fname
            plt.close('all')

    def plot_R_hist(self):
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
        plt.xlabel('Energy Resolution at 400nm')
        if self.outpath==None:
            plt.show()
        else:
            fname=self.outpath+self.params['figdir']+'calsol_' + self.calFN.tstamp +'_R_Estimates.png'
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
    paramFile = sys.argv[1]
    calFNs, params = getCalFileNames(paramFile)

    for calFN in calFNs:
        try:
            diag_obj=waveCal_diagnostic(calFN,params,save=False)
            #diag_obj.make_R_array()
            #diag_obj.plot_R_array()
            #diag_obj.plot_nlaser_array()
            #diag_obj.plot_R_hist()
            diag_obj.plot_parabolafit_hist()
            #diag_obj.plot_parameterfit_hist()
            #diag_obj.plot_sigmas()
            #diag_obj.plot_amps()
            plt.show()
            plt.close('all')
            del diag_obj
        except:
            pass




