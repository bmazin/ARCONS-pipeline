from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import QtGui
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

class PopUp(QMainWindow):
    def __init__(self, parent=None,plotFunc=None,clickFunc=None,title=''):
        if parent==None:
            app = QApplication(sys.argv)
        QMainWindow.__init__(self, parent)
        self.setWindowTitle(title)
        self.plotFunc = plotFunc
        self.clickFunc = clickFunc
        self.create_main_frame()
        self.create_status_bar()
        if plotFunc != None:
            plotFunc(fig=self.fig,axes=self.axes0)
        self.show()
        if parent==None:
            app.exec_()

    def clickCanvas(self,event):
        self.clickFunc(fig=self.fig,axes=self.axes0,event=event)
        self.canvas.draw()

    def create_main_frame(self):
        self.main_frame = QWidget()
      # Create the mpl Figure and FigCanvas objects. 
        self.dpi = 100
        self.fig = Figure((5, 5), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        self.axes0 = self.fig.add_subplot(111)
        if self.clickFunc != None:
            cid=self.canvas.mpl_connect('button_press_event', self.clickCanvas)

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

def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()
