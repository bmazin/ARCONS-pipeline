'''
Author: Alex Walter                             Date: April 23, 2014

This code creates master cal files

Input:
    paramFile - same as wavecal parameter file
Output:
    - Array map showing 3 peak or 2 peak solutions for each pixel
    - PDF file of plots showing location of average blue peak for each pixel
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
from waveCal import fitData,print_guesses,getCalFileNames
from drift_diagnostics import *

class master_waveCal:
    """
        Class to combine multiple wavecals into a master wavecal
    """

    def __init__(self,drift_object,times_to_combine_str,save_plots=True):
        """
            Creates a drift_diagnositc object to load wavecal soln data into memory. Set's up outpath
        """
        self.drift = drift_object
        self.times_to_combine = [[strpdate2num("%Y%m%d-%H%M%S")(tstamp1),strpdate2num("%Y%m%d-%H%M%S")(tstamp2)] for [tstamp1,tstamp2] in times_to_combine_str]
        self.drift.mastercalFNs, self.params = getCalFileNames(self.drift.params,'mastercal_','_drift')
        self.save_plots=save_plots

        intermDir=self.params['intermDir']
        outdir=self.params['outdir']
        if intermDir is None or intermDir is '':
            intermDir = os.getenv('MKID_PROC_PATH', default="/Scratch")
        if outdir is None or outdir is '':
            outdir = '/waveCalSolnFiles'
        self.outpath=intermDir+outdir+os.sep+self.drift.driftFNs[0].run+os.sep+'master_cals'
        try:
            os.mkdir(self.outpath)
            os.mkdir(self.outpath+'/drift_study')
            if self.save_plots:
                os.mkdir(self.outpath+'/figs')
        except:
            pass


    def create_master_peak_data(self):
        """
            Calculates average location of the blue/red/IR peaks for the wavecals within a master cal. 
        """
        if not hasattr(self.drift, 'blue_xOffset'):
            self.drift.populate_peak_data()
        if not hasattr(self.drift, 'sigma'):
            self.drift.populate_sig_range_data()

        print "Calculating Master Cal Data for "+str(self.drift.driftFNs[0].run)+'...'

        self.master_sigma = -1.0*np.ones((self.drift.beammap.shape[0],self.drift.beammap.shape[1],len(self.times_to_combine)))
        self.master_range_min = -1.0*np.ones((self.drift.beammap.shape[0],self.drift.beammap.shape[1],len(self.times_to_combine)))
        self.master_range_max = -1.0*np.ones((self.drift.beammap.shape[0],self.drift.beammap.shape[1],len(self.times_to_combine)))
        self.master_blue = np.zeros((self.drift.beammap.shape[0],self.drift.beammap.shape[1],len(self.times_to_combine)))
        self.master_red = np.zeros((self.drift.beammap.shape[0],self.drift.beammap.shape[1],len(self.times_to_combine)))
        self.master_IR = np.zeros((self.drift.beammap.shape[0],self.drift.beammap.shape[1],len(self.times_to_combine)))


        for i in range(len(self.times_to_combine)):
        #for i in [0]:
            indices_of_wavecals = (np.where((np.asarray(self.drift.timeArray) >= self.times_to_combine[i][0]) * (np.asarray(self.drift.timeArray) <= self.times_to_combine[i][1])))[0]
            #print indices_of_wavecals

            std_blue = np.zeros(self.drift.beammap.shape)
            std_red = np.zeros(self.drift.beammap.shape)
            std_IR = np.zeros(self.drift.beammap.shape)
            

            for row in range(self.drift.beammap.shape[0]):
                for col in range(self.drift.beammap.shape[1]):
                    weights = self.drift.sigma[row,col,indices_of_wavecals] > 0.0
                    if np.sum(weights) <= 0.0:
                        #No solutions in this pixel
                        #print "No solutions for ("+str(col)+', '+str(row)+')'
                        continue

                    self.master_range_max[row,col,i] = max(self.drift.range_max[row,col,indices_of_wavecals][(np.where(weights))[0]])
                    self.master_range_min[row,col,i] = min(self.drift.range_min[row,col,indices_of_wavecals][(np.where(weights))[0]])
                    self.master_sigma[row,col,i] = np.average(self.drift.sigma[row,col,indices_of_wavecals], weights=weights)
                    weights = self.drift.blue_xOffset[row,col,indices_of_wavecals] < 0.0
                    #print weights
                    #if np.sum(weights) > 0.0:
                    self.master_blue[row,col,i] = np.average(self.drift.blue_xOffset[row,col,indices_of_wavecals], weights=weights)
                    std_blue[row,col] = np.std(self.drift.blue_xOffset[row,col,indices_of_wavecals][(np.where(weights))[0]])

                    if len(indices_of_wavecals) >=5:
                        weights = weights * (self.drift.blue_xOffset[row,col,indices_of_wavecals] > 2.5*std_blue[row,col]) * (self.drift.blue_xOffset[row,col,indices_of_wavecals] < -2.5*std_blue[row,col])
                        if np.sum(weights) > 0.0:       #This should always be true...
                            self.master_blue[row,col,i] = np.average(self.drift.blue_xOffset[row,col,indices_of_wavecals], weights=weights)
                            self.master_sigma[row,col,i] = np.average(self.drift.sigma[row,col,indices_of_wavecals], weights=weights)
                            self.master_range_max[row,col,i] = max(self.drift.range_max[row,col,indices_of_wavecals][(np.where(weights))[0]])


                    weights = self.drift.red_xOffset[row,col,indices_of_wavecals] < 0.0
                    #if np.sum(weights) > 0.0:
                    self.master_red[row,col,i] = np.average(self.drift.red_xOffset[row,col,indices_of_wavecals], weights=weights)
                    std_red[row,col] = np.std(self.drift.red_xOffset[row,col,indices_of_wavecals][(np.where(weights))[0]])
                    if len(indices_of_wavecals) >=5:
                        weights = weights * (self.drift.red_xOffset[row,col,indices_of_wavecals] > 2.5*std_red[row,col]) * (self.drift.red_xOffset[row,col,indices_of_wavecals] < -2.5*std_red[row,col])
                        if np.sum(weights) > 0.0:       #This should always be true...
                            self.master_red[row,col,i] = np.average(self.drift.red_xOffset[row,col,indices_of_wavecals], weights=weights)
                            self.master_range_min[row,col,i] = min(self.drift.range_min[row,col,indices_of_wavecals][(np.where(weights))[0]])

                    weights = self.drift.IR_xOffset[row,col,indices_of_wavecals] < 0.0
                    if np.sum(weights) > 0.0:
                        self.master_IR[row,col,i] = np.average(self.drift.IR_xOffset[row,col,indices_of_wavecals], weights=weights)
                        std_IR[row,col] = np.std(self.drift.IR_xOffset[row,col,indices_of_wavecals][(np.where(weights))[0]])
                        if len(indices_of_wavecals) >=5:
                            weights = weights * (self.drift.IR_xOffset[row,col,indices_of_wavecals] > 2.5*std_IR[row,col]) * (self.drift.IR_xOffset[row,col,indices_of_wavecals] < -2.5*std_IR[row,col])
                            if np.sum(weights) > 0.0:       #This should always be true...
                                self.master_IR[row,col,i] = np.average(self.drift.IR_xOffset[row,col,indices_of_wavecals], weights=weights)
                                self.master_range_min[row,col,i] = min(self.drift.range_min[row,col,indices_of_wavecals][(np.where(weights))[0]])

        print '\tDone'

    def create_master_cal_data(self):
        print "Fitting Parabolas to Master Cal Data for "+str(self.drift.driftFNs[0].run)+'...'

        self.master_parab_const=-1.0*np.ones((self.drift.beammap.shape[0],self.drift.beammap.shape[1],len(self.times_to_combine)))
        self.master_parab_lin=-1.0*np.ones((self.drift.beammap.shape[0],self.drift.beammap.shape[1],len(self.times_to_combine)))
        self.master_parab_quad=-1.0*np.ones((self.drift.beammap.shape[0],self.drift.beammap.shape[1],len(self.times_to_combine)))

        for i in range(len(self.times_to_combine)):
            n_sol = 0
            for row in range(self.drift.beammap.shape[0]):
                for col in range(self.drift.beammap.shape[1]):

                    if self.master_sigma[row,col,i] == -1:
                        continue

                    wavelengths = [self.params['bluelambda'],self.params['redlambda'], self.params['irlambda']]     #Angstroms
                    energies = [self.params['h'] * self.params['c'] / (x * self.params['ang2m']) for x in wavelengths]             #eV
                    laser_amps=np.asarray([self.master_blue[row,col,i],self.master_red[row,col,i],self.master_IR[row,col,i]])

                    # parameter_guess [constant, linear term (slope), quadratic term (perturbation to straight line)]
                    parameter_guess = [0.0,(energies[0]-energies[1])*1.0/(laser_amps[0]-laser_amps[1]), -10.**-5.]
                    parameter_lowerlimit=[None]*3
                    parameter_upperlimit=[None]*3
                    if self.master_IR[row,col,i]==0.0:      #No IR solution
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
                    
                    parameter_fit, redchi2gauss2, mpperr = fitData(laser_amps,energies,parameter_guess,parameter_lowerlimit,parameter_upperlimit,model='parabola',make_plot=False,verbose=False)

                    num_param_fails = np.sum(parameter_fit==parameter_guess)+np.sum(parameter_fit==parameter_lowerlimit) + np.sum(parameter_fit==parameter_upperlimit)
                    if parameter_fit[2]==0.0:       # no quadratic term, just fitting line
                        num_param_fails = 0
                    if mpperr==None or num_param_fails>0.0:
                        print "Unable to fit parabola: ("+str(col)+', '+str(row)+')'
                        print_guesses(parameter_guess,[None]*3,[None]*3,parameter_fit)
                        print "energies: "+str(energies)
                        print "laser_amps: "+str(laser_amps)
                        raise RuntimeError("Parabola Fit Failure. mpperr: "+str(mpperr))

                    self.master_parab_const[row,col,i] = parameter_fit[0]
                    self.master_parab_lin[row,col,i] = parameter_fit[1]
                    self.master_parab_quad[row,col,i] = parameter_fit[2]
                    n_sol+=1
            print '\tFound '+str(n_sol)+' parabola solutions for time interval '+str(i)

        print '\tDone'




    def write_master(self):
        self.write_master_cals()
        self.write_master_drifts()
        if self.save_plots:
            self.make_nlaser_plots()
            self.make_xOffsets_plots()

    def make_nlaser_plots(self):
        print "Saving nlaser plots..."
        for i in range(len(self.times_to_combine)):
            index_of_first_wavecal=np.argsort(self.drift.timeArray)[np.argmax(np.asarray(self.drift.timeArray)[np.argsort(self.drift.timeArray)] >= self.times_to_combine[i][0])]
            pltfilename = self.outpath + '/figs/mastercal_'+self.drift.driftFNs[index_of_first_wavecal].tstamp + '_nlaserPlot.png'
            self.nlaser = 1.0*(self.master_blue[:,:,i] < 0.0) + 1.0*(self.master_red[:,:,i] < 0.0) +1.0*(self.master_IR[:,:,i] < 0.0)
            plotArray( self.nlaser, showMe=False, cbar=True,plotFileName=pltfilename, plotTitle='Number of Lasers for Fit')
            print '\tSaving n laser plot to: '+pltfilename
        plt.close('all')
        print "\tDone"

    def make_xOffsets_plots(self):
        print "Saving xOffset plots..."
        for i in range(len(self.times_to_combine)):
        #for i in [1]:
            indices_of_wavecals = (np.where((np.asarray(self.drift.timeArray) >= self.times_to_combine[i][0]) * (np.asarray(self.drift.timeArray) <= self.times_to_combine[i][1])))[0]
            midtime= np.average([np.amax(np.asarray(self.drift.timeArray)[indices_of_wavecals]), np.amin(np.asarray(self.drift.timeArray)[indices_of_wavecals])])
            index_of_first_wavecal=np.argsort(self.drift.timeArray)[np.argmax(np.asarray(self.drift.timeArray)[np.argsort(self.drift.timeArray)] >= self.times_to_combine[i][0])]
            pltfilename = self.outpath + '/figs/mastercal_'+self.drift.driftFNs[index_of_first_wavecal].tstamp + '_xOffsets.pdf'
            pp = PdfPages(pltfilename)
            mpl.rcParams['font.size'] = 4
            n_plots_x = 3
            n_plots_y = 4
            n_plots_per_page = n_plots_x*n_plots_y
            plotcounter = 0

            fig_pdf = plt.figure(figsize=(8.25, 10), dpi=100)
                

            for row in range(self.drift.beammap.shape[0]):
                for col in range(self.drift.beammap.shape[1]):
                #for col in [23]:
                    if self.master_blue[row,col,i] >= 0.0:
                        continue
                    plt.subplot(n_plots_y,n_plots_x,plotcounter%n_plots_per_page+1)
                    plt.errorbar(np.asarray(self.drift.timeArray)[indices_of_wavecals],self.drift.blue_xOffset[row,col,indices_of_wavecals],yerr=self.drift.blue_sigma[row,col,indices_of_wavecals],fmt='bo', markersize=2., capsize=0, elinewidth=0.4, ecolor='b')
                    plt.errorbar(np.asarray(self.drift.timeArray)[indices_of_wavecals],self.drift.red_xOffset[row,col,indices_of_wavecals],yerr=self.drift.red_sigma[row,col,indices_of_wavecals],fmt='ro', markersize=2., capsize=0, elinewidth=0.4, ecolor='r')
                    plt.errorbar(np.asarray(self.drift.timeArray)[indices_of_wavecals],self.drift.IR_xOffset[row,col,indices_of_wavecals],yerr=self.drift.IR_sigma[row,col,indices_of_wavecals],fmt='ko', markersize=2., capsize=0, elinewidth=0.4, ecolor='k')

                    plt.errorbar([midtime],[self.master_blue[row,col,i]], xerr = [np.amax(np.asarray(self.drift.timeArray)[indices_of_wavecals])-midtime],fmt='go')
                    plt.errorbar([midtime],[self.master_red[row,col,i]], xerr = [np.amax(np.asarray(self.drift.timeArray)[indices_of_wavecals])-midtime],fmt='go')
                    plt.errorbar([midtime],[self.master_IR[row,col,i]], xerr = [np.amax(np.asarray(self.drift.timeArray)[indices_of_wavecals])-midtime],fmt='go')

                    plt.title('('+str(col)+', '+str(row)+') --> '+self.drift.beammap[row][col])

                    if ((plotcounter +1) % n_plots_per_page == 0):
                        pp.savefig(fig_pdf)
                        fig_pdf = plt.figure(figsize=(8.25, 10), dpi=100)
                    plotcounter+=1



            if plotcounter % n_plots_per_page!=0:
                try:
                    pp.savefig(fig_pdf)
                except:
                    pass
            plt.close('all')
            pp.close()
            print "\tSaved to: "+pltfilename
            print '\t#: '+str(plotcounter)

        print '\tDone'


    def write_master_cals(self):
        print "Writing Master Cal Files for "+str(self.drift.driftFNs[0].run)+'...'
        WaveCalSoln_Description = {
            "roach"     : UInt16Col(),      # ROACH board number
            "pixelnum"  : UInt16Col(),      # pixel number on the roach
            "pixelrow"  : UInt16Col(),      # physical x location - from beam map
            "pixelcol"  : UInt16Col(),      # physical y location 
            "polyfit"   : Float64Col(3),    # polynomial to convert from phase amplitude to wavelength float 64 precision
            "sigma"     : Float64Col(),     # 1 sigma (Gaussian width) in eV, for blue peak
            "solnrange" : Float32Col(2),    # start and stop wavelengths for the fit in Angstroms
            "wave_flag" : UInt16Col()}      # flag to indicate if pixel is good (0), unallocated (1), dead (2), or failed during wave cal fitting (2+) 


        for i in range(len(self.times_to_combine)):

            index_of_first_wavecal=np.argsort(self.drift.timeArray)[np.argmax(np.asarray(self.drift.timeArray)[np.argsort(self.drift.timeArray)] >= self.times_to_combine[i][0])]

            master_cal_file_name = self.outpath+"/mastercal_" + self.drift.driftFNs[index_of_first_wavecal].tstamp + '.h5'
            try:
                waveCalFile = tables.openFile(master_cal_file_name,mode='w')
            except:
                pass
            # Create a group that branches from the root node, save cal data table in this group
            calgroup = waveCalFile.createGroup(waveCalFile.root, 'wavecal', 'Table of calibration parameters for each pixel')

            # Create a table instance under group calgroup with node name "calsoln"
            # the WaveCalSoln class declared before is the description parameter to define the columns of the table
            caltable = waveCalFile.createTable(calgroup, 'calsoln', WaveCalSoln_Description, title='Wavelength Cal Table')

            for row in range(self.drift.beammap.shape[0]):
                for col in range(self.drift.beammap.shape[1]):
                    trow = caltable.row
                    trow['pixelrow'] = row
                    trow['pixelcol'] = col
                    pixel_name=self.drift.beammap[row][col]
                    trow['roach'] = int(pixel_name.split('/')[1][1:])
                    trow['pixelnum'] = int(pixel_name.split('/')[2][1:])
                    trow['wave_flag'] = int(max(-1.0*self.master_sigma[row,col,i],0))
                    trow['polyfit'] = [self.master_parab_const[row,col,i],self.master_parab_lin[row,col,i],self.master_parab_quad[row,col,i]]
                    trow['sigma'] = self.master_sigma[row,col,i]
                    trow['solnrange'] = [self.master_range_min[row,col,i],self.master_range_max[row,col,i]]
                    #trow['fitmodel'] = self.params['model_type']
                    trow.append()


            print "\tWrote to: "+master_cal_file_name

            # flush the table's I/O buffer to write the data to disk
            caltable.flush()
            # close the file, flush all remaining buffers
            waveCalFile.close()

        print '\tDone'

    def write_master_drifts(self):
        print "Writing Master Drift Files for "+str(self.drift.driftFNs[0].run)+'...'
        DriftTab_Description = {
            "pixelrow"      : UInt16Col(),                  # physical x location - from beam map
            "pixelcol"      : UInt16Col(),                  # physical y location 
            "avglaserphase" : Float64Col(3)}                # average phase height location of lasers
            #"perrors"       : Float64Col(num_params)}       # the errors on the fits

        DriftTimeTab_Description = {
            "wavecal"   : StringCol(22)}                # List the names of the wavecals used


        for i in range(len(self.times_to_combine)):

            index_of_first_wavecal=np.argsort(self.drift.timeArray)[np.argmax(np.asarray(self.drift.timeArray)[np.argsort(self.drift.timeArray)] >= self.times_to_combine[i][0])]

            master_drift_file_name = self.outpath+"/drift_study/mastercal_" + self.drift.driftFNs[index_of_first_wavecal].tstamp + '_drift.h5'

            try:
                driftCalFile = tables.openFile(master_drift_file_name,mode='w')
            except:
                pass
            
            driftgroup = driftCalFile.createGroup(driftCalFile.root, 'params_drift', 'Table of parameters for drift study')
            drifttable = driftCalFile.createTable(driftgroup, 'driftparams', DriftTab_Description, title='Drift Params Table')

            #for pix in range(len(self.drift_pixRowArr)):
            for row in range(self.drift.beammap.shape[0]):
                for col in range(self.drift.beammap.shape[1]):
                    trow = drifttable.row
                    trow['pixelrow'] = row
                    trow['pixelcol'] = col
                    trow['avglaserphase'] = [self.master_blue[row,col,i],self.master_red[row,col,i],self.master_IR[row,col,i]]
                    #trow['perrors'] = 
                    trow.append()

            # flush the table's I/O buffer to write the data to disk
            drifttable.flush()

            indices_of_wavecals = (np.where((np.asarray(self.drift.timeArray) >= self.times_to_combine[i][0]) * (np.asarray(self.drift.timeArray) <= self.times_to_combine[i][1])))[0]
            drifttimetable = driftCalFile.createTable(driftgroup, 'drifttimes', DriftTimeTab_Description, title='List of wavecals averaged')
            for k in range(len(indices_of_wavecals)):
                trow = drifttimetable.row
                trow['wavecal'] = "cal_" + self.drift.driftFNs[indices_of_wavecals[k]].tstamp + '.h5'
                trow.append()

            drifttimetable.flush()

            print "\tWrote to: "+master_drift_file_name

            # close the file, flush all remaining buffers
            driftCalFile.close()

        print '\tDone'




if __name__ == "__main__":
    try:
        paramFile = sys.argv[1]
    except IndexError:
        paramFile=os.getenv('PYTHONPATH',default=os.path.expanduser('~')+'/ARCONS-pipeline/')+'params/waveCal.dict'
        #paramFile = '/home/abwalter/ARCONS-pipeline/params/waveCal.dict'
        print "Loading parameters from: "+paramFile
    calFNs, params = getCalFileNames(paramFile,'calsol_','_drift.h5',getAll=True)   #Automatically grabs all cal files
    drift_object = drift_object(calFNs,params)
    drift_object.populate_peak_data()


    #PAL2012
    #times_to_combine_str = [['20121206-030038', '20121206-124822'],
    #                        ['20121207-023126', '20121207-135957'],
    #                        ['20121208-024053', '20121208-031249'],
    #                        ['20121208-031718', '20121208-133819'], 
    #                        ['20121209-021608', '20121209-133420'], 
    #                        ['20121211-020440', '20121211-135845'],
    #                        ['20121212-023030', '20121212-133822']]
    times_to_combine_str = [['20121206-030039', '20121206-034635'],         # Dec 5
                            ['20121206-044321', '20121206-115709'],         # Dec 5, realigned laser and retuned array
                            ['20121207-023127', '20121207-062840'],         # Dec 6
                            ['20121207-064749', '20121207-132325'],         # Dec 6, retuned
                            ['20121208-024054', '20121208-031248'],         # Dec 7, optimal filters from day 1, 2
                            ['20121208-031719', '20121208-063741'],         # Dec 7, template FIR filters
                            ['20121208-070505', '20121208-133818'],         # Dec 7, retuned
                            ['20121209-021609', '20121209-060704'],         # Dec 8
                            ['20121209-062404', '20121209-131132'],         # Dec 8, retuned
                            ['20121211-020441', '20121211-072632'],         # Dec 10
                            ['20121211-074031', '20121211-135844'],         # Dec 10, retuned
                            ['20121212-023031', '20121212-063518'],         # Dec 11
                            ['20121212-065247', '20121212-133821']]         # Dec 11, retuned
    #times_to_combine_str = [['20121206-030039', '20121212-133821']]

    #PAL2013
    #times_to_combine_str = [['20131209-093153','20131209-132225']]

    master = master_waveCal(drift_object, times_to_combine_str)
    master.create_master_peak_data()
    master.create_master_cal_data()
    master.write_master()
    drift_object.plot_laser_xOffset()
    drift_object.plot_numSols_map()
    #drift.populate_drift_fluctuations()
    #drift.plot_fluct_map()
    #drift.hist_fluct()







