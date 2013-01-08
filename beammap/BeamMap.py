import numpy as np
from tables import *
import time
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages
import scipy.signal as signal
from scipy import optimize
import scipy.stats as stats

from PyQt4.QtGui import *
from PyQt4.QtGui import *
#from beammap_gui import Ui_beammap_gui


# Define the various classes and functions needed for the beam mapping
# Define BeamMapper class - this allows you to change the location of a peak if there is a mistake
class BeamMapper():    
    # Initialize variables needed within the class
    def __init__(self,xtime,ytime,xfilelength,yfilelength):
        self.crx_median = np.zeros((2024,xtime))
        self.cry_median = np.zeros((2024,ytime))
        self.crx = np.zeros(((xfilelength,2024,xtime)))
        self.cry = np.zeros(((yfilelength,2024,ytime)))
        self.flag = np.zeros(2024)
        self.peakpos = np.zeros((2,2024))
    # Try to find a peak position by manually selecting an approximate peak location
    def on_click(self,event):
        # If x sweep plot (top plot) is clicked
        if(event.y > 250):
            self.xvals=np.arange(len(self.crx_median[pixelno][:]))
            self.xpeakguess=event.xdata
            self.xfitstart=max([self.xpeakguess-20,0])
            self.xfitend=min([self.xpeakguess+20,len(self.xvals)])
            params = fitgaussian(self.crx_median[pixelno][self.xfitstart:self.xfitend],self.xvals[self.xfitstart:self.xfitend])
            self.xfit = gaussian(params,self.xvals)
            self.peakpos[0,self.pixelno]=params[0]
        # If y sweep plot (bottom plot) is clicked
        else:
            self.yvals=np.arange(len(self.cry_median[pixelno][:]))
            self.ypeakguess=event.xdata
            self.yfitstart=max([self.ypeakguess-20,0])
            self.yfitend=min([self.ypeakguess+20,len(self.yvals)])
            params = fitgaussian(self.cry_median[pixelno][self.yfitstart:self.yfitend],self.yvals[self.yfitstart:self.yfitend])
            self.yfit = gaussian(params,self.yvals)
            self.peakpos[1,self.pixelno]=params[0]
    # Connect to plot
    def connect(self):
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)

# Define a standard Gaussian distribution function
def gaussian(pars, x):
    center, width, height, back = pars
    width = float(width)
    return back + height*np.exp(-(((center-x)/width)**2)/2)

# Define an error function between data and a Gaussian
def errorfunction(params, data, x):
    errorfunction = data - gaussian(params,x)
    return errorfunction

# Find an optimal Guassian fit for the data, return parameters of that Gaussian
def fitgaussian(data,x):
    params=(x.mean(),2.*(x[1]-x[0]),data.max(), 0.)
    p, success = optimize.leastsq(errorfunction, params, args=(data, x))
    return p

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

class Sweep_Number(QWidget):
    
    def __init__(self):
        super(Sweep_Number, self).__init__()
        
        self.initUI()
        
    def initUI(self):      

        self.btnx = QPushButton('X Sweeps', self)
        self.btnx.move(20, 20)
        self.btnx.clicked.connect(self.showDialogX)

        self.btny = QPushButton('Y Sweeps', self)
        self.btny.move(20, 50)
        self.btny.clicked.connect(self.showDialogY)

        self.btn1 = QPushButton('Frequency Path', self)
        self.btn1.move(20, 80)
        self.btn1.clicked.connect(self.freq_dialog)

        self.btn1 = QPushButton('Save Path', self)
        self.btn1.move(20, 110)
        self.btn1.clicked.connect(self.save_dialog)
        
        self.lex = QLineEdit(self)
        self.lex.move(130, 20)

        self.ley = QLineEdit(self)
        self.ley.move(130, 50)

        self.le1 = QLineEdit(self)
        self.le1.setGeometry(130, 80, 400, 20)

        self.le2 = QLineEdit(self)
        self.le2.setGeometry(130, 110, 400, 20)

        self.sweep_count = [0,0]

        self.freqpath = ''
        self.savepath = ''

        self.okbtn = QPushButton('OK', self)
        self.okbtn.move(20, 150)
        self.okbtn.clicked.connect(self.close)
        
        self.setGeometry(300, 300, 600, 200)
        self.setWindowTitle('Inputs')
        self.show()

    def freq_dialog(self):
        text =QFileDialog.getExistingDirectory()
        self.freqpath=str(text)
        self.le1.setText(str(text))
        

    def save_dialog(self):
        text =QFileDialog.getExistingDirectory()
        self.savepath=str(text)
        self.le2.setText(str(text))
        
        
    def showDialogX(self):
        
        text, ok = QInputDialog.getText(self, 'Input Dialog', 
            'Number of X Sweeps:')
        
        if ok:
            self.lex.setText(str(text))
            self.sweep_count[0] = int(text)

    def showDialogY(self):
        
        text, ok = QInputDialog.getText(self, 'Input Dialog', 
            'Number of Y Sweeps:')
        
        if ok:
            self.ley.setText(str(text))
            self.sweep_count[1] = int(text)
        
def numberrun():
    
    app = QApplication(sys.argv)
    ex1 = Sweep_Number()
    app.exec_()
    return [ex1.sweep_count[0],ex1.sweep_count[1],ex1.freqpath,ex1.savepath]

class Sweep_DialogX(QWidget):    
    def __init__(self):
        super(Sweep_DialogX, self).__init__()       
        self.initUI()
        
    def initUI(self):

        self.btn = QPushButton('Add', self)
        self.btn.move(20, 20)
        self.btn.clicked.connect(self.file_dialog)

        self.okbtn = QPushButton('Clear', self)
        self.okbtn.move(20, 50)
        self.okbtn.clicked.connect(self.clear_sweeps)
        
        self.okbtn = QPushButton('OK', self)
        self.okbtn.move(20, 100)
        self.okbtn.clicked.connect(self.close)

        self.le = []
        for i in range(input_params[0]):
            self.le.append(QLineEdit(self))
            self.le[i].setGeometry(130, 22+22*i,400,20)
        
        self.x_array = []
        
        self.setGeometry(300, 300, 600, 150)
        self.setWindowTitle('Choose X sweep files')
        self.show()

    def clear_sweeps(self):
        for i in range(len(self.x_array)):
            self.le[i].setText('')
        self.x_array = []
        
            
    def file_dialog(self):
        text =QFileDialog.getOpenFileName(parent=None, caption=str("Choose Sweep File"), filter=str("H5 (*.h5)")) 
        self.x_array.append(str(text))
        for i in range(len(self.x_array)):
            self.le[i].setText(self.x_array[i])

class Sweep_DialogY(QWidget):    
    def __init__(self):
        super(Sweep_DialogY, self).__init__()       
        self.initUI()
        
    def initUI(self):

        self.btn = QPushButton('Add', self)
        self.btn.move(20, 20)
        self.btn.clicked.connect(self.file_dialog)

        self.okbtn = QPushButton('Clear', self)
        self.okbtn.move(20, 50)
        self.okbtn.clicked.connect(self.clear_sweeps)
        
        self.okbtn = QPushButton('OK', self)
        self.okbtn.move(20, 100)
        self.okbtn.clicked.connect(self.close)

        self.le = []
        for i in range(input_params[1]):
            self.le.append(QLineEdit(self))
            self.le[i].setGeometry(130, 22+22*i,400,20)
        
        self.y_array = []
        
        self.setGeometry(300, 300, 600, 150)
        self.setWindowTitle('Choose Y sweep files')
        self.show()

    def clear_sweeps(self):
        for i in range(len(self.y_array)):
            self.le[i].setText('')
        self.y_array = []
        
            
    def file_dialog(self):
        text =QFileDialog.getOpenFileName(parent=None, caption=str("Choose Sweep File"), filter=str("H5 (*.h5)")) 
        self.y_array.append(str(text))
        for i in range(len(self.y_array)):
            self.le[i].setText(self.y_array[i])

def xrun():
    
    appx = QApplication(sys.argv)
    ex = Sweep_DialogX()
    appx.exec_()
    return ex.x_array

def yrun():
    
    appy = QApplication(sys.argv)
    ey = Sweep_DialogY()
    appy.exec_()
    return ey.y_array

input_params = numberrun()
xsweep = xrun()
ysweep = yrun()
print xsweep
print ysweep

# Specify input/output directory and files
#path = '/Users/kids/desktop/Beammapping/'

#xsweep = []
#xsweep.append('obs_20121129-014040.h5')
#xsweep.append('obs_20121129-014802.h5')

#ysweep = []
#ysweep.append('obs_20121129-005856.h5')
#ysweep.append('obs_20121129-010903.h5')

number_of_roaches = 8
roach_pixel_count = np.zeros(number_of_roaches)
for roachno in xrange(0,number_of_roaches):
    roach_pixel_count[roachno] = file_len(input_params[2] + '/ps_freq%i.txt' %roachno)-1

# Load the input files
# X sweep data
h5file_x = []
ts_x = []
exptime_x = []
for i in range(len(xsweep)):
#    h5file_x.append(openFile(path + xsweep[i], mode = 'r'))
    h5file_x.append(openFile(xsweep[i], mode = 'r'))
    try:
        ts_x.append(int(h5file_x[i].root.header.header.col('unixtime')[0]))
    except KeyError:
        ts_x.append(int(h5file_x[i].root.header.header.col('ut')[0]))
    exptime_x.append(int(h5file_x[i].root.header.header.col('exptime')[0]))
# Y sweep data
h5file_y = []
ts_y = []
exptime_y = []
for i in range(len(ysweep)):
#    h5file_y.append(openFile(path + ysweep[i], mode = 'r'))
    h5file_y.append(openFile(ysweep[i], mode = 'r'))
    try:
        ts_y.append(int(h5file_y[i].root.header.header.col('unixtime')[0]))
    except KeyError:
        ts_y.append(int(h5file_y[i].root.header.header.col('ut')[0]))
    exptime_y.append(int(h5file_y[i].root.header.header.col('exptime')[0]))

# Print start and sweep durations, also pixel selection commands
for i in range(len(xsweep)):
    print 'Start Time %i = ' %i,ts_x[i], ts_y[i]
for i in range(len(xsweep)):

    print 'exptime_x %i =' %i,exptime_x[i],'and exptime_y %i =' %i, exptime_y[i]
print '"A" = Accept, "D" = Delete, "Q" = Quit, "X" = X Only, "Y" = Y Only'

# Create a BeamMapper instance
# Go through each of the pixels (originally (0,2024))
for roachno in xrange(0,number_of_roaches):
    map = BeamMapper(exptime_x[0],exptime_y[0],len(xsweep),len(ysweep))
    f=open(input_params[3] + '/r%i.pos' %roachno,'w')
    for pixelno in xrange(0,int(roach_pixel_count[roachno])):
        map.pixelno=pixelno
        map.flag[pixelno] = pixelno
    
        # Store the x data into crx
        pn = []
        data = np.empty(((len(xsweep),exptime_x[0])), dtype = object)
        for i in range(len(xsweep)):
            pn.append('/r%d/p%d/t%d' % ( roachno ,pixelno, ts_x[i]))       
        try:
            for i in range(len(xsweep)):
                data[i][:] = h5file_x[i].root._f_getChild(pn[i]).read()
            for j in xrange(0,exptime_x[0]):
                median_array = []
                for i in range(len(xsweep)):
                    median_array.append(len(data[i][j]))
                map.crx_median[pixelno][j] = np.median(median_array)
                for i in range(len(xsweep)):
                    map.crx[i][pixelno][j] = len(data[i][j])
        except:
            pass
     
        # Store the y data into cry
        pn = []
        data = np.empty(((len(ysweep),exptime_y[0])), dtype = object)
        for i in range(len(ysweep)):
            pn.append('/r%d/p%d/t%d' % ( roachno ,pixelno, ts_y[i]))
        try:
            for i in range(len(ysweep)):
                data[i][:] = h5file_y[i].root._f_getChild(pn[i]).read()
            for j in xrange(exptime_y[0]):
                median_array = []                
                for i in range(len(ysweep)):
                    median_array.append(len(data[i][j]))             
                map.cry_median[pixelno][j] = np.median(median_array)                
                for i in range(len(ysweep)):
                    map.cry[i][pixelno][j] = len(data[i][j])
        except:
            pass

        print 'roach %d, pixel %d' % (roachno, pixelno)
        map.fig = plt.figure()

        # do fit of x-position
        map.ax = map.fig.add_subplot(211)
        map.xvals=np.arange(len(map.crx_median[pixelno][:]))
        plt.title(str(pixelno))
        plt.plot(map.xvals,map.crx_median[pixelno][:])       
        map.xpeakguess=np.where(map.crx_median[pixelno][:] == map.crx_median[pixelno][:].max())[0][0]
        map.xfitstart=max([map.xpeakguess-20,0])
        map.xfitend=min([map.xpeakguess+20,len(map.xvals)])
        params_x = fitgaussian(map.crx_median[pixelno][map.xfitstart:map.xfitend],map.xvals[map.xfitstart:map.xfitend])
        print 'x: [center, width, height, back] =', params_x
        map.xfit = gaussian(params_x,map.xvals)
        plt.plot(map.xvals, map.xfit)
        for i in range(len(xsweep)):
            plt.plot(map.xvals,map.crx[i][pixelno][:],alpha = .2)
        map.ax.set_ylim(params_x[3]-30,params_x[2]+75)
        map.peakpos[0,pixelno]=params_x[0]
            
        # do fit of y-position
        map.ax = map.fig.add_subplot(212)
        map.yvals=np.arange(len(map.cry_median[pixelno][:]))
        plt.title(str(pixelno+1))
        plt.plot(map.yvals,map.cry_median[pixelno][:])     
        map.ypeakguess=np.where(map.cry_median[pixelno][:] == map.cry_median[pixelno][:].max())[0][0]
        map.yfitstart=max([map.ypeakguess-20,0])
        map.yfitend=min([map.ypeakguess+20,len(map.yvals)])
        params_y = fitgaussian(map.cry_median[pixelno][map.yfitstart:map.yfitend],map.yvals[map.yfitstart:map.yfitend])
        print 'y: [center, width, height, back] =', params_y
        map.yfit = gaussian(params_y,map.yvals)
        plt.plot(map.yvals, map.yfit)
        for i in range(len(ysweep)):
            plt.plot(map.yvals,map.cry[i][pixelno][:],alpha = .2)
        map.ax.set_ylim(params_y[3]-30,params_y[2]+75)
        map.peakpos[1,pixelno]=params_y[0]

        # AUTOMATICALLY DELETE OBVIOUSLY BAD PIXELS
        
#        if((params_x[0] == 9.5 and params_y[0] == 9.5) or (params_x[1] >= 100 and params_y[1] >= 100)):
#            map.peakpos[0,pixelno]=0
#            map.peakpos[1,pixelno]=0
#            map.flag[pixelno] = -1*pixelno
#            plt.close()
#            f=open(path+'Sweep_1/r%i.pos' %roachno,'a')
#            f.write(str(map.peakpos[0,pixelno])+'\t'+str(map.peakpos[1,pixelno])+'\t'+str(int(map.flag[pixelno]))+'\n')
#            f.close()
#        else:
        map.connect()
        while True:
            # Accept pixel with 'a'
            response=raw_input('>')
            if(response == 'a'):
                print 'Response: ', response
                plt.close()
                f=open(input_params[3]+'/r%i.pos' %roachno,'a')
                f.write(str(map.peakpos[0,pixelno])+'\t'+str(map.peakpos[1,pixelno])+'\t'+str(int(map.flag[pixelno]))+ '\n')
                f.close()
                break
            # Delete pixel with 'd'
            elif(response == 'd'):
                print 'Response: ', response
                map.peakpos[0,pixelno]=0
                map.peakpos[1,pixelno]=0
                map.flag[pixelno] = -1*pixelno
                plt.close()
                f=open(input_params[3]+'/r%i.pos' %roachno,'a')
                f.write(str(map.peakpos[0,pixelno])+'\t'+str(map.peakpos[1,pixelno])+'\t'+str(int(map.flag[pixelno]))+ '\n')
                f.close()
                break
            # Control whether only the x or y pixel is selected
            # x is okay
            elif(response == 'x'):
                print 'Response: ', response
                map.peakpos[1,pixelno]=0
                map.flag[pixelno] = -1*pixelno
                plt.close()
                f=open(input_params[3]+'/r%i.pos' %roachno,'a')
                f.write(str(map.peakpos[0,pixelno])+'\t'+str(map.peakpos[1,pixelno])+'\t'+str(int(map.flag[pixelno]))+ '\n')
                f.close()
                break
            # y is okay
            elif(response == 'y'):
                print 'Response: ', response
                map.peakpos[0,pixelno]=0
                map.flag[pixelno] = -1*pixelno
                plt.close()
                f=open(input_params[3]+'/r%i.pos' %roachno,'a')
                f.write(str(map.peakpos[0,pixelno])+'\t'+str(map.peakpos[1,pixelno])+'\t'+str(int(map.flag[pixelno]))+ '\n')
                f.close()
                break
            # Quit program with 'q'
            elif(response == 'q'):
                sys.exit(0)
            plt.show()
            
h5file_x.close()
h5file_y.close()
f.close()
