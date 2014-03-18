'''
Author: Alex Walter                             Date: March 18, 2014

This code analyzes the wavecal solutions over time using the _drift.h5 files.

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
#from fitFunctions import *
from utils import *



def getDriftFileNames(paramFile):
    params = readDict(paramFile)
    params.readFromFile(paramFile)

    #Load in cal files from param file
    run = params['run']
    sunsetDate = params['sunsetDate']
    tStampList = params['calTimeStamps']
    walkPath = None

    if sunsetDate == None:
        intermDir=params['intermdir']
        outdir=params['outdir']
        if intermDir is None or intermDir is '':
            intermDir = os.getenv('INTERM_DIR', default="/Scratch")+os.sep
        if outdir is None or outdir is '':
            outdir = 'waveCalSolnFiles/'
        walkPath=intermDir+outdir+run+os.sep
    elif tStampList == None:
        intermDir=params['intermdir']
        outdir=params['outdir']
        if intermDir is None or intermDir is '':
            intermDir = os.getenv('INTERM_DIR', default="/Scratch")+os.sep
        if outdir is None or outdir is '':
            outdir = 'waveCalSolnFiles/'
        walkPath=intermDir+outdir+run+os.sep+sunsetDate+os.sep
    if walkPath != None:
        print 'Using all files from: '+walkPath
        calFNs = []
        for root,dirs,files in os.walk(walkPath):
            for f in files:
                if f.endswith('_drift.h5'):
                    d=(root.split(run)[-1]).split('/')[1]
                    t=f.split('_')[1]
                    calFNs.append(FileName(run=run, date=d,tstamp=t))
    else:
        calFNs = [FileName(run=run, date=sunsetDate,tstamp=tStamp) for tStamp in tStampList]
    return calFNs, params


class drift_ana_plots:
    def __init__(self,driftFNs,params,save=True):

        self.driftFNs = driftFNs
        self.params=params
        
        self.outpath=None
        if save:
            intermDir=params['intermdir']
            outdir=params['outdir']
            if intermDir is None or intermDir is '':
                intermDir = os.getenv('INTERM_DIR', default="/Scratch")+os.sep
            if outdir is None or outdir is '':
                outdir = 'waveCalSolnFiles/'
            self.outpath=intermDir+outdir+driftFNs[0].run+os.sep+'drift_figs/'
            try:
                os.mkdir(self.outpath)
            except:
                pass

        temp=tables.openFile(driftFNs[0].cal(), mode='r')
        self.beammap = temp.root.beammap.beamimage.read()
        temp.close()


    def populate_data(self):
        print "Collecting Data..."

        self.timeArray=[strpdate2num("%Y%m%d-%H%M%S")(driftFN.tstamp) for driftFN in self.driftFNs]

        self.blue_xOffset=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))
        self.blue_sigma=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))
        self.red_xOffset=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))
        self.red_sigma=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))
        self.IR_xOffset=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))
        self.IR_sigma=np.zeros((self.beammap.shape[0],self.beammap.shape[1],len(self.driftFNs)))

        for i in range(len(driftFNs)):
            try:
                driftFile=tables.openFile(driftFNs[i].calDriftInfo(),mode='r')
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
                print '\tUnable to open: '+driftFNs[i].calDriftInfo()
                #raise
        print "\tDone."

    def plot_laser_xOffset(self,minSol=2):
        print "Making xOffset Plots..."
        if self.outpath!=None:
            pp = PdfPages(self.outpath+'driftAna_plots.pdf')
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
                    if self.outpath==None:
                        plt.show()
                    else:
                        pp.savefig(fig1)
                        #print 'Saved page: '+str(plotcounter)
                plotcounter+=1
        if plotcounter % n_plots_per_page != 0:
            if self.outpath==None:
                plt.show()
            else:
                pp.savefig(fig1)

        if self.outpath!=None:
            pp.close()

        print "\tFound " +str(plotcounter)+" pixels with at least "+str(minSol)+" solutions"

    def plot_fluct_map(self,minSol=2):
        if self.outpath==None:
            showMe=True
            pltfilename=None
        else:
            showMe=False
            pltfilename=self.outpath+'fluct_map.png'
        fluctMap=np.zeros(self.beammap.shape)
        for i in range(len(fluctMap)):
            for k in range(len(fluctMap[i])):
                sol=np.where(self.blue_sigma[i,k,:] > 0.0)[0]
                if len(sol)>minSol:
                    #Add fluctuation value to map
                    pass

        #plotArray( fluctMap, colormap=mpl.cm.hot, showMe=showMe, cbar=True,plotFileName=pltfilename, plotTitle='Fluctuations in solutions')

    def hist_fluct(self,minSol=2):
        pass


    def plot_numSols_map(self):
        print "Mapping number of Solutions in pixel"
        if self.outpath==None:
            showMe=True
            pltfilename=None
        else:
            showMe=False
            pltfilename=self.outpath+'numSols_map.png'
        numSols=np.zeros(self.beammap.shape)
        for i in range(len(numSols)):
            for k in range(len(numSols[i])):
                numSols[i,k]=len(np.where(self.blue_sigma[i,k,:] > 0.0)[0])
        plotArray( numSols, colormap=mpl.cm.hot, showMe=showMe, cbar=True,plotFileName=pltfilename, plotTitle='Number of solutions out of %i' %(len(self.driftFNs)))
        print "\tDone."



if __name__ == '__main__':
    paramFile = sys.argv[1]
    driftFNs, params = getDriftFileNames(paramFile)
    
    try:
        drift_ana = drift_ana_plots(driftFNs,params,save=True)
        drift_ana.populate_data()
        drift_ana.plot_laser_xOffset()
        drift_ana.plot_fluct_map()
        drift_ana.plot_numSols_map()
    except:
        print "ERROR!"
        raise


