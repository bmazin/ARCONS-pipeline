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


    def plotArray(self,image,normNSigma=3,title='',showColorBar=True,**kwargs):
        self.image = image
        if not 'vmax' in kwargs:
            goodImage = image[np.isfinite(image)]
            kwargs['vmax'] = np.mean(goodImage)+normNSigma*np.std(goodImage)
        if not 'cmap' in kwargs:
            defaultCmap=matplotlib.cm.gnuplot2
            defaultCmap.set_bad('0.15')
            kwargs['cmap'] = defaultCmap
        if not 'origin' in kwargs:
            kwargs['origin'] = 'lower'
            
        self.handleMatshow = self.axes.matshow(image,**kwargs)
        if showColorBar:
            self.fig.cbar = self.fig.colorbar(self.handleMatshow)
            cid = self.fig.canvas.mpl_connect('scroll_event', self.onscroll_cbar)
            cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick_cbar)
        self.axes.set_title(title)
        cid = self.fig.canvas.mpl_connect('motion_notify_event', self.hoverCanvas)
        self.draw()

    def show(self):
        super(PopUp,self).show()
        if self.parent == None:
            self.app.exec_()

    def hoverCanvas(self,event):
        if event.inaxes is self.axes:
            col = int(round(event.xdata))
            row = int(round(event.ydata))
            if row < np.shape(self.image)[0] and col < np.shape(self.image)[1]:
                self.status_text.setText('({:d},{:d}) {}'.format(row,col,self.image[row,col]))
            

    def create_status_bar(self):
        self.status_text = QLabel("Awaiting orders.")
        self.statusBar().addWidget(self.status_text, 1)
        
    def onscroll_cbar(self, event):
        if event.inaxes is self.fig.cbar.ax:
            increment=0.05
            currentClim = self.fig.cbar.mappable.get_clim()
            if event.button == 'up':
                newClim = (currentClim[0],(1.+increment)*currentClim[1])
            if event.button == 'down':
                newClim = (currentClim[0],(1.-increment)*currentClim[1])
            self.fig.cbar.mappable.set_clim(newClim)
            self.fig.canvas.draw()

    def onclick_cbar(self,event):
        if event.inaxes is self.fig.cbar.ax:
            if event.button == 1:
                self.fig.oldClim = self.fig.cbar.mappable.get_clim()
                self.fig.cbar.mappable.set_clim(self.fig.oldClim[0],event.ydata*self.fig.oldClim[1])
                self.fig.canvas.draw()
            if event.button == 3:
                self.fig.oldClim = self.fig.cbar.mappable.get_clim()
                self.fig.cbar.mappable.set_clim(self.fig.oldClim[0],1/event.ydata*self.fig.oldClim[1])
                self.fig.canvas.draw()

def plotArray(*args,**kwargs):
    #Waring: Does not play well with matplotlib state machine style plotting!
    block = kwargs.pop('block',False)
    def f(*args,**kwargs):
        form = PopUp(showMe=False)
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

def main():
    print 'non-blocking PopUp A, close when done'
    plotArray(title='A',image=np.arange(12).reshape(4,3))
    plotArray(title='C',image=np.arange(12).reshape(3,4))
    print 'blocking PopUp, close when done'
    form = PopUp(showMe=False,title='B')
    form.plotArray(np.arange(9).reshape(3,3))
    form.show()
    print 'done'


if __name__ == "__main__":
    main()
