'''
Author: Alex Walter                             Date: March 18, 2014

This code analyzes the wavecal solutions over time using the _drift.h5 files.
You can use it by itself but usually it is automatically called by master_waveCal.py

Input:
    paramFile - same as wavecal parameter file
Output:
    - Array map showing number of solutions for each pixel
    - PDF file of plots showing location of blue peak with error bars equal to sigma of blue peak for each pixel and each wavecal soln
'''



import sys, os
import tables
from tables import *
import numpy as np
import matplotlib.pylab as plt
from matplotlib.dates import strpdate2num
from matplotlib.backends.backend_pdf import PdfPages

from util.readDict import readDict
from util.FileName import FileName
from fitFunctions import *
from utils import *
from waveCal import fitData,print_guesses,getCalFileNames
from waveCal_diagnostics import waveCal_diagnostic


class drift_object:
    def __init__(self,driftFNs,params):

        self.driftFNs = driftFNs
        self.params=params
        
        intermDir=params['intermDir']
        outdir=params['outdir']
        if intermDir is None or intermDir is '':
            intermDir = os.getenv('MKID_PROC_PATH', default="/Scratch")
        if outdir is None or outdir is '':
            outdir = '/waveCalSolnFiles'
        else:
            self.calSoln_str=[intermDir+outdir+'/calsol_'+driftFN.tstamp+'.h5' for driftFN in self.driftFNs]
        #self.outpath=intermDir+outdir+driftFNs[0].run+os.sep+'drift_figs/'
        self.outpath=intermDir+outdir+os.sep+driftFNs[0].run+os.sep+'master_cals/figs/'


        temp=tables.openFile(driftFNs[0].cal(), mode='r')
        self.beammap = temp.root.beammap.beamimage.read()
        temp.close()

        self.timeArray=[strpdate2num("%Y%m%d-%H%M%S")(driftFN.tstamp) for driftFN in self.driftFNs]
        if len(self.timeArray)!=len(np.unique(self.timeArray)):
            print "We're reading in a duplicate file!"
            raise RuntimeError


    def populate_sig_range_data(self):
        print "Collecting Sigma and Range Data for "+str(self.driftFNs[0].run)+'...'

        self.sigma=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))
        self.range_min=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))
        self.range_max=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))

        for i in range(len(self.driftFNs)):
            try:
                calFile=tables.openFile(self.driftFNs[i].calSoln(),mode='r')
                cal_row = calFile.root.wavecal.calsoln.cols.pixelrow[:]    
                cal_col = calFile.root.wavecal.calsoln.cols.pixelcol[:]    
                cal_range = calFile.root.wavecal.calsoln.cols.solnrange[:]
                cal_sigma = calFile.root.wavecal.calsoln.cols.sigma[:]
                for p in range(len(cal_row)):
                    self.range_min[cal_row[p],cal_col[p],i]=(cal_range[p])[0]
                    self.range_max[cal_row[p],cal_col[p],i]=(cal_range[p])[1]
                    self.sigma[cal_row[p],cal_col[p],i]=cal_sigma[p]
                calFile.close()
            except:
                print '\tUnable to open: '+self.driftFNs[i].calSoln()
                #raise
        print "\tDone."

    def populate_gaussparam_data(self):
        print "Collecting Gaussian Param Data for "+str(self.driftFNs[0].run)+'...'
        self.gauss_params = np.zeros((len(self.driftFNs),self.beammap.shape[0],self.beammap.shape[1],12))
        
        for i in range(len(self.driftFNs)):
            try:
                driftFile=tables.openFile(self.driftFNs[i].calDriftInfo(),mode='r')
                drift_row = driftFile.root.params_drift.driftparams.cols.pixelrow[:]    
                drift_col = driftFile.root.params_drift.driftparams.cols.pixelcol[:]    
                drift_params = driftFile.root.params_drift.driftparams.cols.gaussparams[:]
                for p in range(len(drift_row)):
                    #print len(self.gauss_params[i,drift_row[p],drift_col[p],:])
                    #print len(drift_params)
                    self.gauss_params[i,drift_row[p],drift_col[p],:] = drift_params[p]
                driftFile.close()
            except:
                print '\tUnable to open: '+self.driftFNs[i].calDriftInfo()
                raise
        print "\tDone."
        #return self.gauss_params

    def populate_peak_data(self):
        print "Collecting Drift Data for "+str(self.driftFNs[0].run)+'...'



        self.blue_xOffset=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))
        self.blue_sigma=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))
        self.red_xOffset=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))
        self.red_sigma=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))
        self.IR_xOffset=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))
        self.IR_sigma=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))

        for i in range(len(self.driftFNs)):
            try:
                driftFile=tables.openFile(self.driftFNs[i].calDriftInfo(),mode='r')
                drift_row = driftFile.root.params_drift.driftparams.cols.pixelrow[:]    
                drift_col = driftFile.root.params_drift.driftparams.cols.pixelcol[:]    
                drift_params = driftFile.root.params_drift.driftparams.cols.gaussparams[:]
                for p in range(len(drift_row)):
                    self.blue_xOffset[drift_row[p],drift_col[p],i]=(drift_params[p])[1]
                    self.blue_sigma[drift_row[p],drift_col[p],i]=(drift_params[p])[0]
                    self.red_xOffset[drift_row[p],drift_col[p],i]=(drift_params[p])[4]
                    self.red_sigma[drift_row[p],drift_col[p],i]=(drift_params[p])[3]
                    mask= (drift_params[p])[8] > 0
                    self.IR_xOffset[drift_row[p],drift_col[p],i]=(drift_params[p])[7]*mask
                    self.IR_sigma[drift_row[p],drift_col[p],i]=(drift_params[p])[6]*mask
                driftFile.close()
            except:
                print '\tUnable to open: '+self.driftFNs[i].calDriftInfo()
                #raise
        print "\tDone."

    def populate_cal_data(self):
        print "Collecting Cal Data for "+str(self.driftFNs[0].run)+'...'


        self.parab_const=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))
        self.parab_lin=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))
        self.parab_quad=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))

        for i in range(len(self.driftFNs)):
            try:
                calFile=tables.openFile(self.driftFNs[i].calSoln(),mode='r')
                cal_row = calFile.root.wavecal.calsoln.cols.pixelrow[:]    
                cal_col = calFile.root.wavecal.calsoln.cols.pixelcol[:]    
                cal_params = calFile.root.wavecal.calsoln.cols.polyfit[:]
                for p in range(len(cal_row)):
                    self.parab_const[cal_row[p],cal_col[p],i]=(cal_params[p])[0]
                    self.parab_lin[cal_row[p],cal_col[p],i]=(cal_params[p])[1]
                    self.parab_quad[cal_row[p],cal_col[p],i]=(cal_params[p])[2]
                calFile.close()
            except:
                print '\tUnable to open: '+self.driftFNs[i].calSoln()
                #raise
        print "\tDone."

    def populate_R_data(self):
        print "Collecting R Data for "+str(self.driftFNs[0].run)+'...'
        self.r_blue=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))
        for i in range(len(self.driftFNs)):
            try:
                calFile=tables.openFile(self.driftFNs[i].calSoln(),mode='r')
                cal_row = calFile.root.wavecal.calsoln.cols.pixelrow[:]    
                cal_col = calFile.root.wavecal.calsoln.cols.pixelcol[:]    
                cal_params = calFile.root.wavecal.calsoln.cols.polyfit[:]
                cal_sigma = calFile.root.wavecal.calsoln.cols.sigma[:]
            except:
                print '\tUnable to open: '+self.driftFNs[i].calSoln()
                return
            try:
                driftFile=tables.openFile(self.driftFNs[i].calDriftInfo(),mode='r')
                drift_row = driftFile.root.params_drift.driftparams.cols.pixelrow[:]    
                drift_col = driftFile.root.params_drift.driftparams.cols.pixelcol[:]    
                drift_params = driftFile.root.params_drift.driftparams.cols.gaussparams[:]
            except:
                print '\tUnable to open: '+self.driftFNs[i].calDriftInfo()
                return
            for k in range(len(cal_sigma)):
                if cal_sigma[k]>0:
                    drift_ind = np.where((drift_row==cal_row[k]) * (drift_col==cal_col[k]))[0][0]
                    peak_fit = drift_params[drift_ind]
                    blue_energy = (parabola(cal_params[k],x=np.asarray([peak_fit[1]]),return_models=True))[0][0]
                    self.r_blue[cal_row[k],cal_col[k],i]=blue_energy/(self.params['fwhm2sig']*cal_sigma[k])
            calFile.close()
            driftFile.close()

        print "\tDone."

    def populate_master_sig_range_data(self):
        if hasattr(self, 'master_sigma'):
            return
        if not hasattr(self, 'masterFNs'):
            self.masterFNs,p = getCalFileNames(self.params,'mastercal_','_drift.h5')
        if len(self.masterFNs)==0:
            print "No master cal files found!"
            return

        print "Collecting Master Sigma and Range Data for "+str(self.driftFNs[0].run)+'...'
        self.master_sigma=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.masterFNs)))
        self.master_range_min=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.masterFNs)))
        self.master_range_max=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.masterFNs)))

        for i in range(len(self.masterFNs)):
            try:
                calFile=tables.openFile(self.masterFNs[i].mastercalSoln(),mode='r')
                cal_row = calFile.root.wavecal.calsoln.cols.pixelrow[:]    
                cal_col = calFile.root.wavecal.calsoln.cols.pixelcol[:]    
                cal_range = calFile.root.wavecal.calsoln.cols.solnrange[:]
                cal_sigma = calFile.root.wavecal.calsoln.cols.sigma[:]
                for p in range(len(cal_row)):
                    self.master_range_min[cal_row[p],cal_col[p],i]=(cal_range[p])[0]
                    self.master_range_max[cal_row[p],cal_col[p],i]=(cal_range[p])[1]
                    self.master_sigma[cal_row[p],cal_col[p],i]=cal_sigma[p]
                calFile.close()
            except:
                print '\tUnable to open: '+self.masterFNs[i].mastercalSoln()
                #raise
        print "\tDone."

    def populate_master_cal_data(self):
        if hasattr(self, 'master_parab_const'):
            return
        if not hasattr(self, 'masterFNs'):
            self.masterFNs,p = getCalFileNames(self.params,'mastercal_','_drift.h5')
        if len(self.masterFNs)==0:
            print "No master cal files found!"
            return

        print "Collecting Master Cal Data for "+str(self.driftFNs[0].run)+'...'
        self.master_parab_const=-1.0*np.ones((self.beammap.shape[0],self.beammap.shape[1],len(self.masterFNs)))
        self.master_parab_lin=-1.0*np.ones((self.beammap.shape[0],self.beammap.shape[1],len(self.masterFNs)))
        self.master_parab_quad=-1.0*np.ones((self.beammap.shape[0],self.beammap.shape[1],len(self.masterFNs)))

        for i in range(len(self.masterFNs)):
            try:
                calFile=tables.openFile(self.masterFNs[i].mastercalSoln(),mode='r')
                cal_row = calFile.root.wavecal.calsoln.cols.pixelrow[:]    
                cal_col = calFile.root.wavecal.calsoln.cols.pixelcol[:]    
                cal_params = calFile.root.wavecal.calsoln.cols.polyfit[:]
                for p in range(len(cal_row)):
                    self.master_parab_const[cal_row[p],cal_col[p],i]=(cal_params[p])[0]
                    self.master_parab_lin[cal_row[p],cal_col[p],i]=(cal_params[p])[1]
                    self.master_parab_quad[cal_row[p],cal_col[p],i]=(cal_params[p])[2]
                calFile.close()
            except:
                print '\tUnable to open: '+self.masterFNs[i].mastercalSoln()
                #raise
        print "\tDone."

    def populate_master_peak_data(self):
        if hasattr(self, 'master_blue'):
            return
        if not hasattr(self, 'masterFNs'):
            self.masterFNs,p = getCalFileNames(self.params,'mastercal_','_drift.h5')
        if len(self.masterFNs)==0:
            print "No master cal files found!"
            return

        print "Collecting master peak data for "+str(self.driftFNs[0].run)+'...'
        self.master_blue = np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.masterFNs)))
        self.master_red = np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.masterFNs)))
        self.master_IR = np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.masterFNs)))
        for i in range(len(self.masterFNs)):
            try:
                driftFile=tables.openFile(self.masterFNs[i].mastercalDriftInfo(),mode='r')
                drift_row = driftFile.root.params_drift.driftparams.cols.pixelrow[:]    
                drift_col = driftFile.root.params_drift.driftparams.cols.pixelcol[:]    
                drift_params = driftFile.root.params_drift.driftparams.cols.avglaserphase[:]
                for p in range(len(drift_row)):
                    self.master_blue[drift_row[p],drift_col[p],i]=(drift_params[p])[0]
                    self.master_red[drift_row[p],drift_col[p],i]=(drift_params[p])[1]
                    self.master_IR[drift_row[p],drift_col[p],i]=(drift_params[p])[2]
                driftFile.close()
            except:
                print '\tUnable to open: '+self.masterFNs[i].mastercalDriftInfo()
                #raise
        print "\tDone."

    def populate_master_times_data(self):
        if not hasattr(self, 'masterFNs'):
            self.masterFNs,p = getCalFileNames(self.params,'mastercal_','_drift.h5')
        if len(self.masterFNs)==0:
            print "No master cal files found!"
            return

        print "Collecting master wvlcal start and end times for "+str(self.driftFNs[0].run)+'...'
        self.master_start_time = np.zeros(len(self.masterFNs))
        self.master_end_time = np.zeros(len(self.masterFNs))
        for i in range(len(self.masterFNs)):
            try:
                driftFile=tables.openFile(self.masterFNs[i].mastercalDriftInfo(),mode='r')
                drift_times = driftFile.root.params_drift.drifttimes.cols.wavecal[:]
                drift_num_times = [strpdate2num("%Y%m%d-%H%M%S")(fname.split('_')[1].split('.')[0]) for fname in drift_times]
                #print drift_times
                self.master_start_time[i] = np.amin(drift_num_times)
                self.master_end_time[i] = np.amax(drift_num_times)
                driftFile.close()
            except:
                
                print '\tUnable to open: '+self.masterFNs[i].mastercalDriftInfo()

        #print self.master_start_time
        #print self.master_end_time
        print "\tDone."

    def plot_laser_xOffset(self,minSol=2,save=True):
        print "Making xOffset Plots..."
        if save:
            try:
                os.mkdir(self.outpath)
            except:
                pass
            pp = PdfPages(self.outpath+'/driftAna_plots.pdf')
        mpl.rcParams['font.size'] = 6
        n_plots_per_page = 4
        plotcounter = 0

        for i in range(self.beammap.shape[0]):
            for j in range(self.beammap.shape[1]):
                nSols=len(np.where(self.blue_sigma[i,j,:] > 0.0)[0])
                if (nSols < minSol):
                    continue
                if (plotcounter % n_plots_per_page == 0):
                    fig1=plt.figure(figsize=(8.25, 10), dpi=100)
                plt.subplot(4,1,plotcounter%4+1)
                plt.errorbar(self.timeArray,self.blue_xOffset[i,j,:],yerr=self.blue_sigma[i,j,:],fmt='bo', markersize=2., capsize=0, elinewidth=0.4, ecolor='b')
                plt.errorbar(self.timeArray,self.red_xOffset[i,j,:],yerr=self.red_sigma[i,j,:],fmt='ro', markersize=2., capsize=0, elinewidth=0.4, ecolor='r')
                plt.errorbar(self.timeArray,self.IR_xOffset[i,j,:],yerr=self.IR_sigma[i,j,:],fmt='ko', markersize=2., capsize=0, elinewidth=0.4, ecolor='k')
                plt.title('('+str(i)+', '+str(j)+') --> '+self.beammap[i,j])
                fig1.autofmt_xdate()
                
                if (plotcounter +1) % n_plots_per_page == 0:
                    if not save:
                        plt.show()
                    else:
                        pp.savefig(fig1)
                        #print 'Saved page: '+str(plotcounter)
                plotcounter+=1
        if plotcounter % n_plots_per_page != 0:
            if not save:
                plt.show()
            else:
                pp.savefig(fig1)

        if save:
            pp.close()

        print "\tFound " +str(plotcounter)+" pixels with at least "+str(minSol)+" solutions"

    def populate_drift_fluctuations(self):
        """
        Finds the fluctuations in each pixel for each mastercal time period

        The following needs to be run before:
            self.populate_peak_data()
            self.populate_master_peak_data()
            self.populate_master_times_data()
        """
        if hasattr(self, 'drift_fluct'):
            return
        if not hasattr(self, 'masterFNs'):
            self.masterFNs,p = getCalFileNames(self.params,'mastercal_','_drift.h5')
        if len(self.masterFNs)==0:
            print "No master cal files found!"
            return

        self.drift_fluct = np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.masterFNs)))
        for i in range(len(self.masterFNs)):
        #for i in [3]:
            indices_of_wavecals = (np.where((np.asarray(self.timeArray) >= self.master_start_time[i]) * (np.asarray(self.timeArray) <= self.master_end_time[i])))[0]

            flucts = [np.abs(self.blue_xOffset[:,:,j] - self.master_blue[:,:,i]) for j in indices_of_wavecals]
            for j in range(len(flucts)):
                flucts[j][np.where(self.master_blue[:,:,i] != 0)] = flucts[j][np.where(self.master_blue[:,:,i] != 0)]/np.abs(self.master_blue[:,:,i])[np.where(self.master_blue[:,:,i] != 0)]
                flucts[j][np.where(self.master_blue[:,:,i] == 0)] = 1.0
            self.drift_fluct[:,:,i] = np.asarray(np.ma.average(np.asarray(flucts),axis=0,weights=np.asarray(flucts)<1.0))


    def plot_fluct_map(self,save=True):
        if not save:
            showMe=True
            pltfilename=None
        else:
            try:
                os.mkdir(self.outpath)
            except:
                pass
            showMe=False
            
        for i in range(len(self.masterFNs)):
            pltfilename=self.outpath+'/mastercal_'+self.masterFNs[i].tstamp+'_fluctMap.png'
            plotArray(self.drift_fluct[:,:,i], colormap=mpl.cm.jet, showMe=showMe, normMin=0.0,normMax=0.05,cbar=True,plotFileName=pltfilename, plotTitle='Fluctuations in Blue Phase')


    def hist_fluct(self,save=True):
        if not save:
            showMe=True
            pltfilename=None
        else:
            try:
                os.mkdir(self.outpath)
            except:
                pass
            showMe=False
            
        for i in range(len(self.masterFNs)):
        #for i in [3]:
            pltfilename=self.outpath+'/mastercal_'+self.masterFNs[i].tstamp+'_fluctHist.png'
            #hist, bin_edges = np.histogram(self.drift_fluct[:,:,i],bins=200,range=(0.00001,self.drift_fluct[:,:,i].max()))
            hist, bin_edges = np.histogram(self.drift_fluct[:,:,i],bins=100,range=(0.00001,0.05))
            bins = bin_edges[:-1]+(bin_edges[1]-bin_edges[0])/2.0
            plt.plot(bins,hist,linestyle='steps-')
            plt.xlabel('Average pixel fluctuation from mean')
            plt.title('mastercal_'+self.masterFNs[i].tstamp+'.h5')
            if not save:
                plt.show()
            else:
                pass


    def make_numSols_array(self):
        self.numSols=np.zeros(self.beammap.shape)
        for i in range(len(self.numSols)):
            for k in range(len(self.numSols[i])):
                self.numSols[i,k]=len(np.where(self.blue_sigma[i,k,:] > 0.0)[0])

    def plot_numSols_map(self,save=True):
        print "Mapping number of Solutions in pixel"
        if not save:
            showMe=True
            pltfilename=None
        else:
            try:
                os.mkdir(self.outpath)
            except:
                pass
            showMe=False
            pltfilename=self.outpath+'/numSols_map.png'

        if not hasattr(self, 'numSols'):
            self.make_numSols_array()
        plotArray( self.numSols, colormap=mpl.cm.hot, showMe=showMe, cbar=True,plotFileName=pltfilename, plotTitle='Number of solutions out of %i' %(len(self.driftFNs)))
        print "\tDone."



if __name__ == '__main__':
    try:
        paramFile = sys.argv[1]
    except IndexError:
        paramFile=os.getenv('PYTHONPATH',default=os.path.expanduser('~')+'/ARCONS-pipeline/')+'params/waveCal.dict'
        print "Loading parameters from: "+paramFile
    driftFNs, params = getCalFileNames(paramFile,'calsol_','_drift.h5',getAll=True)
    
    try:
        drift = drift_object(driftFNs,params)
        #drift.populate_peak_data()
        #drift.plot_laser_xOffset(save=False)
        ##drift_ana.plot_fluct_map()
        #drift.plot_numSols_map(save=False)

        drift.populate_peak_data()
        drift.populate_master_peak_data()
        drift.populate_master_times_data()
        drift.populate_drift_fluctuations()
        drift.plot_fluct_map(save=False)
        drift.hist_fluct(save=False)


    except:
        print "ERROR!"
        raise


