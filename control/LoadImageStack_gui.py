# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'LoadImageStack_gui.ui'
#
# Created: Mon Nov  4 17:40:32 2013
#      by: PyQt4 UI code generator 4.9.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_LoadImageStack_gui(object):
    def setupUi(self, LoadImageStack_gui):
        LoadImageStack_gui.setObjectName(_fromUtf8("LoadImageStack_gui"))
        LoadImageStack_gui.resize(800, 600)
        self.centralwidget = QtGui.QWidget(LoadImageStack_gui)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        LoadImageStack_gui.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(LoadImageStack_gui)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 25))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        LoadImageStack_gui.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(LoadImageStack_gui)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        LoadImageStack_gui.setStatusBar(self.statusbar)

        self.retranslateUi(LoadImageStack_gui)
        QtCore.QMetaObject.connectSlotsByName(LoadImageStack_gui)

    def retranslateUi(self, LoadImageStack_gui):
        LoadImageStack_gui.setWindowTitle(QtGui.QApplication.translate("LoadImageStack_gui", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))

