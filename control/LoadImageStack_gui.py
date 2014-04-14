# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'LoadImageStack_gui.ui'
#
# Created: Mon Apr 14 13:27:48 2014
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
        LoadImageStack_gui.resize(515, 636)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(LoadImageStack_gui.sizePolicy().hasHeightForWidth())
        LoadImageStack_gui.setSizePolicy(sizePolicy)
        self.plotDock = MPL_Widget(LoadImageStack_gui)
        self.plotDock.setGeometry(QtCore.QRect(10, 40, 491, 521))
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plotDock.sizePolicy().hasHeightForWidth())
        self.plotDock.setSizePolicy(sizePolicy)
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(255, 245, 248))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 245, 248))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 245, 248))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 245, 248))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 245, 248))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 245, 248))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Window, brush)
        self.plotDock.setPalette(palette)
        self.plotDock.setObjectName(_fromUtf8("plotDock"))
        self.rightButton = QtGui.QPushButton(LoadImageStack_gui)
        self.rightButton.setGeometry(QtCore.QRect(380, 600, 92, 27))
        self.rightButton.setAutoDefault(False)
        self.rightButton.setObjectName(_fromUtf8("rightButton"))
        self.leftButton = QtGui.QPushButton(LoadImageStack_gui)
        self.leftButton.setGeometry(QtCore.QRect(40, 600, 92, 27))
        self.leftButton.setAutoDefault(False)
        self.leftButton.setObjectName(_fromUtf8("leftButton"))
        self.frameNumberLine = QtGui.QLineEdit(LoadImageStack_gui)
        self.frameNumberLine.setGeometry(QtCore.QRect(210, 600, 51, 27))
        self.frameNumberLine.setObjectName(_fromUtf8("frameNumberLine"))
        self.maxLabel = QtGui.QLabel(LoadImageStack_gui)
        self.maxLabel.setGeometry(QtCore.QRect(270, 600, 41, 21))
        self.maxLabel.setObjectName(_fromUtf8("maxLabel"))
        self.label = QtGui.QLabel(LoadImageStack_gui)
        self.label.setGeometry(QtCore.QRect(210, 570, 121, 20))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setObjectName(_fromUtf8("label"))
        self.jdLabel = QtGui.QLabel(LoadImageStack_gui)
        self.jdLabel.setGeometry(QtCore.QRect(10, 20, 491, 20))
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setBold(True)
        font.setWeight(75)
        self.jdLabel.setFont(font)
        self.jdLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.jdLabel.setIndent(0)
        self.jdLabel.setObjectName(_fromUtf8("jdLabel"))

        self.retranslateUi(LoadImageStack_gui)
        QtCore.QMetaObject.connectSlotsByName(LoadImageStack_gui)

    def retranslateUi(self, LoadImageStack_gui):
        LoadImageStack_gui.setWindowTitle(QtGui.QApplication.translate("LoadImageStack_gui", "Image Stack Viewer", None, QtGui.QApplication.UnicodeUTF8))
        self.rightButton.setText(QtGui.QApplication.translate("LoadImageStack_gui", "--->", None, QtGui.QApplication.UnicodeUTF8))
        self.leftButton.setText(QtGui.QApplication.translate("LoadImageStack_gui", "<---", None, QtGui.QApplication.UnicodeUTF8))
        self.maxLabel.setText(QtGui.QApplication.translate("LoadImageStack_gui", "/ MAX", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("LoadImageStack_gui", "Frame Number", None, QtGui.QApplication.UnicodeUTF8))
        self.jdLabel.setText(QtGui.QApplication.translate("LoadImageStack_gui", "JD", None, QtGui.QApplication.UnicodeUTF8))

from mpl_pyqt4_widget import MPL_Widget
