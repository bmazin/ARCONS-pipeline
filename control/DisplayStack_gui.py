# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'DisplayStack_gui.ui'
#
# Created: Fri Apr 11 15:00:23 2014
#      by: PyQt4 UI code generator 4.9.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_DisplayStack_gui(object):
    def setupUi(self, DisplayStack_gui):
        DisplayStack_gui.setObjectName(_fromUtf8("DisplayStack_gui"))
        DisplayStack_gui.resize(1101, 797)
        self.centralwidget = QtGui.QWidget(DisplayStack_gui)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.integrationTimeLine = QtGui.QLineEdit(self.centralwidget)
        self.integrationTimeLine.setGeometry(QtCore.QRect(840, 300, 51, 20))
        self.integrationTimeLine.setObjectName(_fromUtf8("integrationTimeLine"))
        self.integrationTimeLabel = QtGui.QLabel(self.centralwidget)
        self.integrationTimeLabel.setGeometry(QtCore.QRect(700, 300, 141, 21))
        self.integrationTimeLabel.setObjectName(_fromUtf8("integrationTimeLabel"))
        self.lowerWavelengthCutoffLabel = QtGui.QLabel(self.centralwidget)
        self.lowerWavelengthCutoffLabel.setGeometry(QtCore.QRect(270, 550, 231, 21))
        self.lowerWavelengthCutoffLabel.setObjectName(_fromUtf8("lowerWavelengthCutoffLabel"))
        self.lowerWavelengthCutoffLine = QtGui.QLineEdit(self.centralwidget)
        self.lowerWavelengthCutoffLine.setEnabled(True)
        self.lowerWavelengthCutoffLine.setGeometry(QtCore.QRect(210, 550, 51, 20))
        self.lowerWavelengthCutoffLine.setAcceptDrops(True)
        self.lowerWavelengthCutoffLine.setFrame(True)
        self.lowerWavelengthCutoffLine.setObjectName(_fromUtf8("lowerWavelengthCutoffLine"))
        self.upperWavelengthCutoffLine = QtGui.QLineEdit(self.centralwidget)
        self.upperWavelengthCutoffLine.setGeometry(QtCore.QRect(210, 580, 51, 20))
        self.upperWavelengthCutoffLine.setObjectName(_fromUtf8("upperWavelengthCutoffLine"))
        self.upperWavelengthCutoffLabel = QtGui.QLabel(self.centralwidget)
        self.upperWavelengthCutoffLabel.setGeometry(QtCore.QRect(270, 580, 231, 21))
        self.upperWavelengthCutoffLabel.setObjectName(_fromUtf8("upperWavelengthCutoffLabel"))
        self.wavelengthCalibrationBox = QtGui.QCheckBox(self.centralwidget)
        self.wavelengthCalibrationBox.setGeometry(QtCore.QRect(210, 300, 231, 22))
        self.wavelengthCalibrationBox.setChecked(True)
        self.wavelengthCalibrationBox.setObjectName(_fromUtf8("wavelengthCalibrationBox"))
        self.flatCalibrationBox = QtGui.QCheckBox(self.centralwidget)
        self.flatCalibrationBox.setGeometry(QtCore.QRect(500, 300, 191, 22))
        self.flatCalibrationBox.setChecked(True)
        self.flatCalibrationBox.setObjectName(_fromUtf8("flatCalibrationBox"))
        self.hotPixelBox = QtGui.QCheckBox(self.centralwidget)
        self.hotPixelBox.setGeometry(QtCore.QRect(700, 360, 191, 22))
        self.hotPixelBox.setChecked(True)
        self.hotPixelBox.setObjectName(_fromUtf8("hotPixelBox"))
        self.line_2 = QtGui.QFrame(self.centralwidget)
        self.line_2.setGeometry(QtCore.QRect(180, -40, 20, 811))
        self.line_2.setFrameShape(QtGui.QFrame.VLine)
        self.line_2.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_2.setObjectName(_fromUtf8("line_2"))
        self.deadPixelBox = QtGui.QCheckBox(self.centralwidget)
        self.deadPixelBox.setGeometry(QtCore.QRect(500, 550, 171, 22))
        self.deadPixelBox.setChecked(True)
        self.deadPixelBox.setObjectName(_fromUtf8("deadPixelBox"))
        self.stackButton = QtGui.QPushButton(self.centralwidget)
        self.stackButton.setGeometry(QtCore.QRect(710, 410, 111, 27))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.stackButton.setFont(font)
        self.stackButton.setObjectName(_fromUtf8("stackButton"))
        self.label = QtGui.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(210, 260, 62, 17))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setObjectName(_fromUtf8("label"))
        self.timeAdjustmentBox = QtGui.QCheckBox(self.centralwidget)
        self.timeAdjustmentBox.setGeometry(QtCore.QRect(700, 330, 201, 22))
        self.timeAdjustmentBox.setChecked(True)
        self.timeAdjustmentBox.setObjectName(_fromUtf8("timeAdjustmentBox"))
        self.obsList = QtGui.QListWidget(self.centralwidget)
        self.obsList.setGeometry(QtCore.QRect(340, 40, 151, 192))
        self.obsList.setObjectName(_fromUtf8("obsList"))
        self.obsTimeLabel = QtGui.QLabel(self.centralwidget)
        self.obsTimeLabel.setGeometry(QtCore.QRect(340, 20, 151, 17))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.obsTimeLabel.setFont(font)
        self.obsTimeLabel.setObjectName(_fromUtf8("obsTimeLabel"))
        self.addButton = QtGui.QPushButton(self.centralwidget)
        self.addButton.setGeometry(QtCore.QRect(510, 80, 92, 27))
        self.addButton.setObjectName(_fromUtf8("addButton"))
        self.removeButton = QtGui.QPushButton(self.centralwidget)
        self.removeButton.setGeometry(QtCore.QRect(510, 120, 92, 27))
        self.removeButton.setObjectName(_fromUtf8("removeButton"))
        self.clearButton = QtGui.QPushButton(self.centralwidget)
        self.clearButton.setGeometry(QtCore.QRect(510, 160, 92, 27))
        self.clearButton.setObjectName(_fromUtf8("clearButton"))
        self.inputList = QtGui.QListWidget(self.centralwidget)
        self.inputList.setGeometry(QtCore.QRect(620, 40, 151, 191))
        self.inputList.setObjectName(_fromUtf8("inputList"))
        self.label_2 = QtGui.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(620, 20, 131, 16))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.wavelengthList = QtGui.QListWidget(self.centralwidget)
        self.wavelengthList.setGeometry(QtCore.QRect(210, 330, 151, 192))
        self.wavelengthList.setObjectName(_fromUtf8("wavelengthList"))
        self.flatList = QtGui.QListWidget(self.centralwidget)
        self.flatList.setGeometry(QtCore.QRect(500, 330, 151, 192))
        self.flatList.setObjectName(_fromUtf8("flatList"))
        self.sunsetList = QtGui.QListWidget(self.centralwidget)
        self.sunsetList.setGeometry(QtCore.QRect(220, 40, 81, 192))
        self.sunsetList.setObjectName(_fromUtf8("sunsetList"))
        self.sunsetLabel = QtGui.QLabel(self.centralwidget)
        self.sunsetLabel.setGeometry(QtCore.QRect(220, 20, 91, 17))
        self.sunsetLabel.setObjectName(_fromUtf8("sunsetLabel"))
        self.line = QtGui.QFrame(self.centralwidget)
        self.line.setGeometry(QtCore.QRect(190, 280, 491, 16))
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.line_3 = QtGui.QFrame(self.centralwidget)
        self.line_3.setGeometry(QtCore.QRect(470, 290, 20, 511))
        self.line_3.setFrameShape(QtGui.QFrame.VLine)
        self.line_3.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_3.setObjectName(_fromUtf8("line_3"))
        self.line_4 = QtGui.QFrame(self.centralwidget)
        self.line_4.setGeometry(QtCore.QRect(670, 290, 21, 481))
        self.line_4.setFrameShape(QtGui.QFrame.VLine)
        self.line_4.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_4.setObjectName(_fromUtf8("line_4"))
        self.loadStackButton = QtGui.QPushButton(self.centralwidget)
        self.loadStackButton.setGeometry(QtCore.QRect(710, 450, 111, 27))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.loadStackButton.setFont(font)
        self.loadStackButton.setObjectName(_fromUtf8("loadStackButton"))
        self.runList = QtGui.QListWidget(self.centralwidget)
        self.runList.setGeometry(QtCore.QRect(20, 40, 141, 192))
        self.runList.setObjectName(_fromUtf8("runList"))
        self.runLabel = QtGui.QLabel(self.centralwidget)
        self.runLabel.setGeometry(QtCore.QRect(20, 20, 91, 17))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.runLabel.setFont(font)
        self.runLabel.setObjectName(_fromUtf8("runLabel"))
        self.targetList = QtGui.QListWidget(self.centralwidget)
        self.targetList.setGeometry(QtCore.QRect(20, 330, 141, 281))
        self.targetList.setObjectName(_fromUtf8("targetList"))
        self.targetLabel = QtGui.QLabel(self.centralwidget)
        self.targetLabel.setGeometry(QtCore.QRect(20, 310, 91, 17))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.targetLabel.setFont(font)
        self.targetLabel.setObjectName(_fromUtf8("targetLabel"))
        self.targetButton = QtGui.QPushButton(self.centralwidget)
        self.targetButton.setGeometry(QtCore.QRect(30, 650, 111, 27))
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.targetButton.setFont(font)
        self.targetButton.setObjectName(_fromUtf8("targetButton"))
        DisplayStack_gui.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(DisplayStack_gui)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1101, 25))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        DisplayStack_gui.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(DisplayStack_gui)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        DisplayStack_gui.setStatusBar(self.statusbar)

        self.retranslateUi(DisplayStack_gui)
        QtCore.QMetaObject.connectSlotsByName(DisplayStack_gui)

    def retranslateUi(self, DisplayStack_gui):
        DisplayStack_gui.setWindowTitle(QtGui.QApplication.translate("DisplayStack_gui", "Display Stack Options", None, QtGui.QApplication.UnicodeUTF8))
        self.integrationTimeLine.setText(QtGui.QApplication.translate("DisplayStack_gui", "10", None, QtGui.QApplication.UnicodeUTF8))
        self.integrationTimeLabel.setText(QtGui.QApplication.translate("DisplayStack_gui", "Integration Time (s)", None, QtGui.QApplication.UnicodeUTF8))
        self.lowerWavelengthCutoffLabel.setText(QtGui.QApplication.translate("DisplayStack_gui", "Lower Wavelength Cutoff (Å)", None, QtGui.QApplication.UnicodeUTF8))
        self.lowerWavelengthCutoffLine.setText(QtGui.QApplication.translate("DisplayStack_gui", "3000", None, QtGui.QApplication.UnicodeUTF8))
        self.upperWavelengthCutoffLine.setText(QtGui.QApplication.translate("DisplayStack_gui", "8000", None, QtGui.QApplication.UnicodeUTF8))
        self.upperWavelengthCutoffLabel.setText(QtGui.QApplication.translate("DisplayStack_gui", "Upper Wavelength Cutoff (Å)", None, QtGui.QApplication.UnicodeUTF8))
        self.wavelengthCalibrationBox.setText(QtGui.QApplication.translate("DisplayStack_gui", "Wavelength Calibration", None, QtGui.QApplication.UnicodeUTF8))
        self.flatCalibrationBox.setText(QtGui.QApplication.translate("DisplayStack_gui", "Flat Calibration", None, QtGui.QApplication.UnicodeUTF8))
        self.hotPixelBox.setText(QtGui.QApplication.translate("DisplayStack_gui", "Hot Pixel Masking", None, QtGui.QApplication.UnicodeUTF8))
        self.deadPixelBox.setText(QtGui.QApplication.translate("DisplayStack_gui", "Remove Dead Pixels", None, QtGui.QApplication.UnicodeUTF8))
        self.stackButton.setText(QtGui.QApplication.translate("DisplayStack_gui", "Stack Images", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("DisplayStack_gui", "Options", None, QtGui.QApplication.UnicodeUTF8))
        self.timeAdjustmentBox.setText(QtGui.QApplication.translate("DisplayStack_gui", "Time Adjustment", None, QtGui.QApplication.UnicodeUTF8))
        self.obsTimeLabel.setText(QtGui.QApplication.translate("DisplayStack_gui", "Observations", None, QtGui.QApplication.UnicodeUTF8))
        self.addButton.setText(QtGui.QApplication.translate("DisplayStack_gui", "Add ->", None, QtGui.QApplication.UnicodeUTF8))
        self.removeButton.setText(QtGui.QApplication.translate("DisplayStack_gui", "<- Remove", None, QtGui.QApplication.UnicodeUTF8))
        self.clearButton.setText(QtGui.QApplication.translate("DisplayStack_gui", "Clear", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("DisplayStack_gui", "Files to process", None, QtGui.QApplication.UnicodeUTF8))
        self.sunsetLabel.setText(QtGui.QApplication.translate("DisplayStack_gui", "Sunset Date", None, QtGui.QApplication.UnicodeUTF8))
        self.loadStackButton.setText(QtGui.QApplication.translate("DisplayStack_gui", "Load Stack", None, QtGui.QApplication.UnicodeUTF8))
        self.runLabel.setText(QtGui.QApplication.translate("DisplayStack_gui", "Runs", None, QtGui.QApplication.UnicodeUTF8))
        self.targetLabel.setText(QtGui.QApplication.translate("DisplayStack_gui", "Targets", None, QtGui.QApplication.UnicodeUTF8))
        self.targetButton.setText(QtGui.QApplication.translate("DisplayStack_gui", "Load Target", None, QtGui.QApplication.UnicodeUTF8))

