from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import QtGui
import matplotlib.pyplot as plt
import numpy as np
import sys
from multiprocessing import Process
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib
from functools import partial

from util.ObsFile import ObsFile
from util.readDict import readDict
from util.FileName import FileName
from fitFunctions import *
from drift_diagnostics import *
from waveCal import fitData, waveCal

class PopUp(QMainWindow):
    def __init__(self, parent=None,plotFunc=None,title='',separateProcess=False, image=None,showMe=True):
        self.parent = parent
        if self.parent == None:
            self.app = QApplication([])
        super(PopUp,self).__init__(parent)
        self.setWindowTitle(title)
        self.plotFunc = plotFunc
        self.create_main_frame(title)
        self.create_status_bar()
        if plotFunc != None:
            plotFunc(fig=self.fig,axes=self.axes)
        if showMe == True:
            self.show()

    def draw(self):
        self.fig.canvas.draw()

    def create_main_frame(self,title):
        self.main_frame = QWidget()
      # Create the mpl Figure and FigCanvas objects. 
        self.dpi = 100
        self.fig = Figure((5, 5), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        self.axes = self.fig.add_subplot(111)
        #self.axes.set_title(title)

        # Create the navigation toolbar, tied to the canvas
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

    def create_status_bar(self):
        self.status_text = QLabel("")
        self.statusBar().addWidget(self.status_text, 1)


    def plotArray(self,image,normNSigma=3,title='',**kwargs):
        self.image = image
        self.handleMatshow = self.axes.matshow(image,cmap=matplotlib.cm.gnuplot2,origin='lower',vmax=np.mean(image)+normNSigma*np.std(image),**kwargs)
        self.fig.cbar = self.fig.colorbar(self.handleMatshow)
        self.axes.set_title(title)
        cid = self.fig.canvas.mpl_connect('scroll_event', partial(onscroll_cbar, self.fig))
        cid = self.fig.canvas.mpl_connect('motion_notify_event', self.hoverCanvas)
        cid = self.fig.canvas.mpl_connect('button_press_event', partial(onclick_cbar, self.fig))
        cid = self.fig.canvas.mpl_connect('button_press_event', self.clickCanvas)
        self.draw()

    def show(self):
        super(PopUp,self).show()
        if self.parent == None:
            self.app.exec_()

    def clickCanvas(self,event):
        if event.inaxes is self.axes:
            col = int(round(event.xdata))
            row = int(round(event.ydata))
            #print row,col
            if self.parent == None:
                try:
                    self.pop_x_offsets(row, col)
                except:
                    print "Failed: ("+str(col)+', '+str(row)+')'
                    raise
            else:
                print "Clicked: ("+str(col)+', '+str(row)+')'


    def pop_x_offsets(self,row,col):
        title = '('+str(col)+', '+str(row)+') --> '+(self.drift.beammap[row,col]).rsplit('/',1)[0]
        pop = PopUp(parent=self,title=title)
        pop.plot_x_offsets(row,col)
        pop.pixel = (row,col)

    def plot_x_offsets(self,i,j):
        self.axes.errorbar(self.parent.drift.timeArray,self.parent.drift.blue_xOffset[i,j,:],yerr=self.parent.drift.blue_sigma[i,j,:],fmt='bo', markersize=2., capsize=0, elinewidth=0.4, ecolor='b')
        self.axes.errorbar(self.parent.drift.timeArray,self.parent.drift.red_xOffset[i,j,:],yerr=self.parent.drift.red_sigma[i,j,:],fmt='ro', markersize=2., capsize=0, elinewidth=0.4, ecolor='r')
        self.axes.errorbar(self.parent.drift.timeArray,self.parent.drift.IR_xOffset[i,j,:],yerr=self.parent.drift.IR_sigma[i,j,:],fmt='ko', markersize=2., capsize=0, elinewidth=0.4, ecolor='k')

        cid = self.fig.canvas.mpl_connect('button_press_event', self.clickCanvas_x_offset)

        

    def clickCanvas_x_offset(self,event):
        if event.inaxes is self.axes and self.mpl_toolbar._active is None:
            if self.parent!=None and self.parent.parent==None:
                closest_time_ind = np.argmin(np.abs(self.parent.drift.timeArray - event.xdata))
                try:
                    self.pop_wavecal_Soln(self.pixel,closest_time_ind)
                except:
                    print "Clicked: ("+str(event.xdata)+', '+str(event.ydata)+')'
                    print "Failed: ("+str(self.parent.drift.timeArray[closest_time_ind])+', '+str(event.ydata)+')'
                    raise

    def pop_wavecal_Soln(self,pixel,calFN_ind):
        waveCal_obj = waveCal(self.parent.drift.driftFNs[calFN_ind],self.parent.drift.params,save_pdf=False,verbose=False,debug=False)
        #print (waveCal_obj.calFN.cal()).rsplit('/',1)[-1]
        title = (waveCal_obj.calFN.cal()).rsplit('/',1)[-1] +' '+ str(pixel) +  ' --> '+(self.parent.drift.beammap[pixel[0],pixel[1]]).rsplit('/',1)[0]
        pop = PopUp(parent=self,title=title)

        pop.pixel = pixel
        pop.waveCal_obj = waveCal_obj
        pop.driftFile=tables.openFile(waveCal_obj.calFN.calDriftInfo(),mode='r')
        
        try:
            pop.plot_pixel_data()
        except IndexError:
            print "There isn't any data to plot"
            

    def plot_pixel_data(self):
        dataDict=self.waveCal_obj.laserCalFile.getTimedPacketList(self.pixel[0],self.pixel[1],timeSpacingCut=self.waveCal_obj.params['danicas_cut'])
        peakHeights=np.asarray(dataDict['peakHeights'])*1.0
        ## less than 'min_amp' per second average count rate
        if dataDict['effIntTime']==0.0 or len(peakHeights)<=(dataDict['effIntTime']*self.waveCal_obj.params['min_count_rate']):
            raise IndexError
        baselines=np.asarray(dataDict['baselines'])*1.0
        peakHeights-=baselines
        biggest_photon = int(min(peakHeights))
        n_inbin,phase_bins=np.histogram(peakHeights,bins=np.abs(biggest_photon),range=(biggest_photon,0))
        phase_bins=(phase_bins+(phase_bins[1]-phase_bins[0])/2.0)[:-1]
        try:
            last_ind = np.where(n_inbin>self.waveCal_obj.params['min_amp'])[0][-1]
        except IndexError:
            last_ind=len(n_inbin)-1
        ## Cut out all the zeros on the right
        n_inbin = n_inbin[:last_ind]
        phase_bins = phase_bins[:last_ind]

        self.axes.plot(phase_bins,n_inbin, 'b.',label="data")

        #Check if this pixel is in drift file
        drift_row = self.driftFile.root.params_drift.driftparams.cols.pixelrow[:]    
        drift_col = self.driftFile.root.params_drift.driftparams.cols.pixelcol[:]    

        args=(np.where(np.multiply(drift_row==self.pixel[0], drift_col==self.pixel[1])))[0]
        if len(args==1):
            #pixel has a solution!
            fit_params = self.driftFile.root.params_drift.driftparams.cols.gaussparams[args[0]]
            print fit_params
            model = self.driftFile.root.params_drift.driftparams.attrs.model_type
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




            

    def hoverCanvas(self,event):
        if event.inaxes is self.axes:
            col = int(round(event.xdata))
            row = int(round(event.ydata))
            if row < np.shape(self.image)[0] and col < np.shape(self.image)[1]:
                self.status_text.setText('({:d},{:d}) {}'.format(row,col,self.image[row,col]))
            

    def create_status_bar(self):
        self.status_text = QLabel("Awaiting orders.")
        self.statusBar().addWidget(self.status_text, 1)




        
def onscroll_cbar(fig, event):
    if event.inaxes is fig.cbar.ax:
        increment=0.05
        currentClim = fig.cbar.mappable.get_clim()
        if event.button == 'up':
            newClim = (currentClim[0],(1.+increment)*currentClim[1])
        if event.button == 'down':
            newClim = (currentClim[0],(1.-increment)*currentClim[1])
        fig.cbar.mappable.set_clim(newClim)
        fig.canvas.draw()

def onclick_cbar(fig,event):
    if event.inaxes is fig.cbar.ax:
        if event.button == 1:
            fig.oldClim = fig.cbar.mappable.get_clim()
            fig.cbar.mappable.set_clim(fig.oldClim[0],event.ydata*fig.oldClim[1])
            fig.canvas.draw()
        if event.button == 3:
            fig.oldClim = fig.cbar.mappable.get_clim()
            fig.cbar.mappable.set_clim(fig.oldClim[0],1/event.ydata*fig.oldClim[1])
            fig.canvas.draw()

def plotArray(*args,**kwargs):
    #Waring: Does not play well with matplotlib state machine style plotting!
    block = kwargs.pop('block',False)
    def f(*args,**kwargs):
        form = PopUp(showMe=False)
        drift = kwargs.pop('drift',None)
        title = kwargs.get('title','')
        form.drift = drift
        form.title = title
        form.plotArray(*args,**kwargs)
        form.show()
    if block==True:
        f(*args,**kwargs)
        return None
    else:
        proc = Process(target=f,args=args,kwargs=kwargs)
        proc.start()
        return proc

def pop(*args,**kwargs):
    #Waring: Does not play well with matplotlib state machine style plotting!
    block = kwargs.pop('block',False)
    def f(*args,**kwargs):
        form = PopUp(showMe=False,*args,**kwargs)
        form.show()
    if block==True:
        f(*args,**kwargs)
        return None
    else:
        proc = Process(target=f,args=args,kwargs=kwargs)
        proc.start()
        return proc

if __name__ == "__main__":
    paramFile = sys.argv[1]
    calFNs, params = getDriftFileNames(paramFile)

    drift = drift_ana(calFNs,params,save=False)
    drift.populate_data()
    drift.make_numSols_array()

    plotArray(drift.numSols,title='Drift Analysis',drift=drift)












