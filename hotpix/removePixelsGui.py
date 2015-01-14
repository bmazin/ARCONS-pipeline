'''
Author: Matt Strader        Date: January 13, 2015
'''

import sys, os
import numpy as np
from PyQt4 import QtCore
from PyQt4 import QtGui
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from functools import partial

from util.FileName import FileName
from util.ObsFile import ObsFile
from util.CalLookupFile import CalLookupFile
from quicklook.ObsFileViewer import ObsFileViewer
from hotpix.manuallyRemovePixels import removePixel

class HotPixGui(ObsFileViewer):
    
    def __init__(self, **kwargs):
        super(HotPixGui,self).__init__(**kwargs)
        self.addClickFunc(self.removeBadPixel)
        
    def removeBadPixel(self,row,col,reason='unknown'):
        reply = QtGui.QMessageBox.question(self, 'Confirm',
                'Mark pixel (x,y)=({},{}) with tag \'{}\'?'.format(col,row,reason), QtGui.QMessageBox.Yes | 
                QtGui.QMessageBox.No, QtGui.QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes:
            if not self.obs.hotPixFileName is None:
                timeMaskPath = self.obs.hotPixFileName
            else:
                raise AttributeError('obs file does not have a timeMask loaded')
            self.obs.hotPixFile.close()
            removePixel(timeMaskPath=timeMaskPath,pixelRow=row,
                        pixelCol=col,reason=reason)
            self.obsMethod('loadHotPixCalFile',timeMaskPath,reasons=['hot pixel','dead pixel','unknown'])
        
if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == '-h' or sys.argv[1] == '--help':
            print 'Usage: {} obsFilePath'
        else:
            obsPath = sys.argv[1]
    else:
        obsPath = None
    
    form = HotPixGui(obsPath=obsPath)
    form.show()
