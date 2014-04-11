#!/bin/python

from PyQt4.QtGui import *
from PyQt4.QtGui import *
from LoadImageStack_gui import Ui_LoadImageStack_gui

'''
Author: Paul Szypryt		Date: November 4, 2013
'''

class LoadImageStack(QMainWindow):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent=None)
        self.ui = Ui_LoadImageStack_gui()
        self.ui.setupUi(self)
