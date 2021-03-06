# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'h5quicklook_v8.ui'
#
# Created: Thu Sep 06 16:20:23 2012
#      by: PyQt4 UI code generator 4.8.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_h5quicklook(object):
    def setupUi(self, h5quicklook):
        h5quicklook.setObjectName(_fromUtf8("h5quicklook"))
        h5quicklook.resize(850, 750)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8("Archon.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        h5quicklook.setWindowIcon(icon)
        self.centralwidget = QtGui.QWidget(h5quicklook)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.browse_button = QtGui.QPushButton(self.centralwidget)
        self.browse_button.setGeometry(QtCore.QRect(180, 10, 91, 31))
        self.browse_button.setObjectName(_fromUtf8("browse_button"))
        self.filename_label = QtGui.QLabel(self.centralwidget)
        self.filename_label.setGeometry(QtCore.QRect(80, 50, 261, 21))
        self.filename_label.setFrameShape(QtGui.QFrame.Box)
        self.filename_label.setText(_fromUtf8(""))
        self.filename_label.setObjectName(_fromUtf8("filename_label"))
        self.label_2 = QtGui.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(20, 50, 62, 17))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.darksub_checkbox = QtGui.QCheckBox(self.centralwidget)
        self.darksub_checkbox.setGeometry(QtCore.QRect(20, 190, 131, 21))
        self.darksub_checkbox.setObjectName(_fromUtf8("darksub_checkbox"))
        self.flat_checkbox = QtGui.QCheckBox(self.centralwidget)
        self.flat_checkbox.setGeometry(QtCore.QRect(20, 230, 111, 21))
        self.flat_checkbox.setObjectName(_fromUtf8("flat_checkbox"))
        self.displayimage_button = QtGui.QPushButton(self.centralwidget)
        self.displayimage_button.setGeometry(QtCore.QRect(40, 340, 121, 41))
        self.displayimage_button.setObjectName(_fromUtf8("displayimage_button"))
        self.saveimage_button = QtGui.QPushButton(self.centralwidget)
        self.saveimage_button.setGeometry(QtCore.QRect(190, 340, 111, 41))
        self.saveimage_button.setObjectName(_fromUtf8("saveimage_button"))
        self.label_3 = QtGui.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(20, 90, 101, 16))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.obstime_lcd = QtGui.QLCDNumber(self.centralwidget)
        self.obstime_lcd.setGeometry(QtCore.QRect(120, 80, 111, 31))
        self.obstime_lcd.setObjectName(_fromUtf8("obstime_lcd"))
        self.label_4 = QtGui.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(240, 90, 101, 17))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.label_5 = QtGui.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(20, 120, 251, 17))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.starttime_spinbox = QtGui.QSpinBox(self.centralwidget)
        self.starttime_spinbox.setGeometry(QtCore.QRect(20, 140, 57, 21))
        self.starttime_spinbox.setMaximum(9999)
        self.starttime_spinbox.setObjectName(_fromUtf8("starttime_spinbox"))
        self.label_6 = QtGui.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(90, 140, 16, 17))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.endtime_spinbox = QtGui.QSpinBox(self.centralwidget)
        self.endtime_spinbox.setGeometry(QtCore.QRect(120, 140, 57, 21))
        self.endtime_spinbox.setMaximum(9999)
        self.endtime_spinbox.setObjectName(_fromUtf8("endtime_spinbox"))
        self.label_7 = QtGui.QLabel(self.centralwidget)
        self.label_7.setGeometry(QtCore.QRect(190, 140, 62, 17))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.label_8 = QtGui.QLabel(self.centralwidget)
        self.label_8.setGeometry(QtCore.QRect(20, 170, 91, 17))
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.label_11 = QtGui.QLabel(self.centralwidget)
        self.label_11.setGeometry(QtCore.QRect(20, 16, 151, 21))
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.tv_image = QtGui.QGraphicsView(self.centralwidget)
        self.tv_image.setGeometry(QtCore.QRect(360, 10, 321, 321))
        self.tv_image.setMouseTracking(True)
        self.tv_image.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.tv_image.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.tv_image.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignTop)
        self.tv_image.setObjectName(_fromUtf8("tv_image"))
        self.histogram_plot = MPL_Widget(self.centralwidget)
        self.histogram_plot.setGeometry(QtCore.QRect(30, 410, 641, 291))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 245, 248))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
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
        self.histogram_plot.setPalette(palette)
        self.histogram_plot.setObjectName(_fromUtf8("histogram_plot"))
        self.label = QtGui.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(420, 380, 101, 21))
        self.label.setObjectName(_fromUtf8("label"))
        self.nbins = QtGui.QSpinBox(self.centralwidget)
        self.nbins.setGeometry(QtCore.QRect(520, 380, 57, 21))
        self.nbins.setMaximum(10000)
        self.nbins.setProperty(_fromUtf8("value"), 100)
        self.nbins.setObjectName(_fromUtf8("nbins"))
        self.plot_pushButton = QtGui.QPushButton(self.centralwidget)
        self.plot_pushButton.setGeometry(QtCore.QRect(580, 370, 101, 41))
        self.plot_pushButton.setObjectName(_fromUtf8("plot_pushButton"))
        self.label_13 = QtGui.QLabel(self.centralwidget)
        self.label_13.setGeometry(QtCore.QRect(730, 10, 101, 20))
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.header_info = QtGui.QTextEdit(self.centralwidget)
        self.header_info.setGeometry(QtCore.QRect(690, 30, 151, 601))
        self.header_info.setObjectName(_fromUtf8("header_info"))
        self.time_up = QtGui.QPushButton(self.centralwidget)
        self.time_up.setGeometry(QtCore.QRect(310, 100, 51, 40))
        self.time_up.setObjectName(_fromUtf8("time_up"))
        self.time_down = QtGui.QPushButton(self.centralwidget)
        self.time_down.setGeometry(QtCore.QRect(310, 130, 51, 41))
        self.time_down.setObjectName(_fromUtf8("time_down"))
        self.satpercent = QtGui.QDoubleSpinBox(self.centralwidget)
        self.satpercent.setGeometry(QtCore.QRect(20, 250, 62, 25))
        self.satpercent.setDecimals(2)
        self.satpercent.setMaximum(100.0)
        self.satpercent.setSingleStep(0.1)
        self.satpercent.setProperty(_fromUtf8("value"), 2.5)
        self.satpercent.setObjectName(_fromUtf8("satpercent"))
        self.label_14 = QtGui.QLabel(self.centralwidget)
        self.label_14.setGeometry(QtCore.QRect(90, 250, 121, 21))
        self.label_14.setObjectName(_fromUtf8("label_14"))
        self.bad_pix = QtGui.QCheckBox(self.centralwidget)
        self.bad_pix.setEnabled(False)
        self.bad_pix.setGeometry(QtCore.QRect(540, 340, 131, 21))
        self.bad_pix.setObjectName(_fromUtf8("bad_pix"))
        self.choosedark = QtGui.QPushButton(self.centralwidget)
        self.choosedark.setGeometry(QtCore.QRect(360, 330, 81, 41))
        self.choosedark.setObjectName(_fromUtf8("choosedark"))
        self.groupBox = QtGui.QGroupBox(self.centralwidget)
        self.groupBox.setGeometry(QtCore.QRect(180, 170, 171, 121))
        self.groupBox.setTitle(_fromUtf8(""))
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.drag_select_radioButton = QtGui.QRadioButton(self.groupBox)
        self.drag_select_radioButton.setGeometry(QtCore.QRect(0, 0, 141, 31))
        self.drag_select_radioButton.setAutoExclusive(False)
        self.drag_select_radioButton.setObjectName(_fromUtf8("drag_select_radioButton"))
        self.mode_buttonGroup = QtGui.QButtonGroup(h5quicklook)
        self.mode_buttonGroup.setObjectName(_fromUtf8("mode_buttonGroup"))
        self.mode_buttonGroup.addButton(self.drag_select_radioButton)
        self.rect_select_radioButton = QtGui.QRadioButton(self.groupBox)
        self.rect_select_radioButton.setGeometry(QtCore.QRect(0, 30, 161, 21))
        self.rect_select_radioButton.setChecked(True)
        self.rect_select_radioButton.setAutoExclusive(False)
        self.rect_select_radioButton.setObjectName(_fromUtf8("rect_select_radioButton"))
        self.mode_buttonGroup.addButton(self.rect_select_radioButton)
        self.rect_x_spinBox = QtGui.QSpinBox(self.groupBox)
        self.rect_x_spinBox.setGeometry(QtCore.QRect(30, 50, 57, 21))
        self.rect_x_spinBox.setMinimum(1)
        self.rect_x_spinBox.setMaximum(32)
        self.rect_x_spinBox.setProperty(_fromUtf8("value"), 1)
        self.rect_x_spinBox.setObjectName(_fromUtf8("rect_x_spinBox"))
        self.rect_y_spinBox = QtGui.QSpinBox(self.groupBox)
        self.rect_y_spinBox.setGeometry(QtCore.QRect(110, 50, 57, 21))
        self.rect_y_spinBox.setMinimum(1)
        self.rect_y_spinBox.setMaximum(32)
        self.rect_y_spinBox.setProperty(_fromUtf8("value"), 1)
        self.rect_y_spinBox.setObjectName(_fromUtf8("rect_y_spinBox"))
        self.label_23 = QtGui.QLabel(self.groupBox)
        self.label_23.setGeometry(QtCore.QRect(20, 50, 16, 21))
        self.label_23.setObjectName(_fromUtf8("label_23"))
        self.label_24 = QtGui.QLabel(self.groupBox)
        self.label_24.setGeometry(QtCore.QRect(100, 50, 16, 21))
        self.label_24.setObjectName(_fromUtf8("label_24"))
        self.circ_r_spinBox = QtGui.QSpinBox(self.groupBox)
        self.circ_r_spinBox.setGeometry(QtCore.QRect(30, 90, 57, 21))
        self.circ_r_spinBox.setMinimum(0)
        self.circ_r_spinBox.setMaximum(16)
        self.circ_r_spinBox.setProperty(_fromUtf8("value"), 0)
        self.circ_r_spinBox.setObjectName(_fromUtf8("circ_r_spinBox"))
        self.label_25 = QtGui.QLabel(self.groupBox)
        self.label_25.setGeometry(QtCore.QRect(20, 90, 16, 21))
        self.label_25.setObjectName(_fromUtf8("label_25"))
        self.circ_select_radioButton = QtGui.QRadioButton(self.groupBox)
        self.circ_select_radioButton.setGeometry(QtCore.QRect(0, 70, 151, 21))
        self.circ_select_radioButton.setAutoExclusive(False)
        self.circ_select_radioButton.setObjectName(_fromUtf8("circ_select_radioButton"))
        self.mode_buttonGroup.addButton(self.circ_select_radioButton)
        self.countlabel = QtGui.QLabel(self.centralwidget)
        self.countlabel.setGeometry(QtCore.QRect(690, 670, 151, 31))
        self.countlabel.setFrameShape(QtGui.QFrame.Box)
        self.countlabel.setObjectName(_fromUtf8("countlabel"))
        self.choosesky = QtGui.QPushButton(self.centralwidget)
        self.choosesky.setGeometry(QtCore.QRect(450, 330, 81, 41))
        self.choosesky.setObjectName(_fromUtf8("choosesky"))
        self.skysub_checkbox = QtGui.QCheckBox(self.centralwidget)
        self.skysub_checkbox.setGeometry(QtCore.QRect(20, 210, 131, 21))
        self.skysub_checkbox.setObjectName(_fromUtf8("skysub_checkbox"))
        self.pixelpath = QtGui.QLabel(self.centralwidget)
        self.pixelpath.setGeometry(QtCore.QRect(690, 640, 151, 31))
        self.pixelpath.setFrameShape(QtGui.QFrame.Box)
        self.pixelpath.setObjectName(_fromUtf8("pixelpath"))
        self.histogram_radio = QtGui.QRadioButton(self.centralwidget)
        self.histogram_radio.setGeometry(QtCore.QRect(320, 380, 102, 21))
        self.histogram_radio.setChecked(False)
        self.histogram_radio.setAutoExclusive(False)
        self.histogram_radio.setObjectName(_fromUtf8("histogram_radio"))
        self.plot_buttonGroup = QtGui.QButtonGroup(h5quicklook)
        self.plot_buttonGroup.setObjectName(_fromUtf8("plot_buttonGroup"))
        self.plot_buttonGroup.addButton(self.histogram_radio)
        self.timestream_radio = QtGui.QRadioButton(self.centralwidget)
        self.timestream_radio.setGeometry(QtCore.QRect(160, 380, 102, 21))
        self.timestream_radio.setChecked(True)
        self.timestream_radio.setAutoExclusive(False)
        self.timestream_radio.setObjectName(_fromUtf8("timestream_radio"))
        self.plot_buttonGroup.addButton(self.timestream_radio)
        self.label_12 = QtGui.QLabel(self.centralwidget)
        self.label_12.setGeometry(QtCore.QRect(10, 380, 171, 21))
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.label_15 = QtGui.QLabel(self.centralwidget)
        self.label_15.setGeometry(QtCore.QRect(280, 380, 21, 21))
        self.label_15.setObjectName(_fromUtf8("label_15"))
        self.increment = QtGui.QSpinBox(self.centralwidget)
        self.increment.setGeometry(QtCore.QRect(310, 80, 51, 25))
        self.increment.setMinimum(1)
        self.increment.setMaximum(100)
        self.increment.setProperty(_fromUtf8("value"), 10)
        self.increment.setObjectName(_fromUtf8("increment"))
        self.histogram_plot_2 = MPL_Widget(self.centralwidget)
        self.histogram_plot_2.setGeometry(QtCore.QRect(30, 720, 641, 291))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 245, 248))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
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
        self.histogram_plot_2.setPalette(palette)
        self.histogram_plot_2.setObjectName(_fromUtf8("histogram_plot_2"))
        self.topplot = QtGui.QComboBox(self.centralwidget)
        self.topplot.setGeometry(QtCore.QRect(690, 820, 131, 31))
        self.topplot.setMaxVisibleItems(3)
        self.topplot.setObjectName(_fromUtf8("topplot"))
        self.topplot.addItem(_fromUtf8(""))
        self.topplot.addItem(_fromUtf8(""))
        self.topplot.addItem(_fromUtf8(""))
        self.bottomplot = QtGui.QComboBox(self.centralwidget)
        self.bottomplot.setGeometry(QtCore.QRect(690, 920, 131, 31))
        self.bottomplot.setEditable(False)
        self.bottomplot.setMaxVisibleItems(3)
        self.bottomplot.setMaxCount(3)
        self.bottomplot.setObjectName(_fromUtf8("bottomplot"))
        self.bottomplot.addItem(_fromUtf8(""))
        self.bottomplot.addItem(_fromUtf8(""))
        self.bottomplot.addItem(_fromUtf8(""))
        self.checkBox = QtGui.QCheckBox(self.centralwidget)
        self.checkBox.setGeometry(QtCore.QRect(700, 860, 141, 21))
        self.checkBox.setObjectName(_fromUtf8("checkBox"))
        self.checkBox_2 = QtGui.QCheckBox(self.centralwidget)
        self.checkBox_2.setGeometry(QtCore.QRect(700, 960, 141, 21))
        self.checkBox_2.setObjectName(_fromUtf8("checkBox_2"))
        self.label_16 = QtGui.QLabel(self.centralwidget)
        self.label_16.setGeometry(QtCore.QRect(700, 900, 111, 20))
        font = QtGui.QFont()
        font.setWeight(75)
        font.setBold(True)
        self.label_16.setFont(font)
        self.label_16.setObjectName(_fromUtf8("label_16"))
        self.label_17 = QtGui.QLabel(self.centralwidget)
        self.label_17.setGeometry(QtCore.QRect(700, 800, 111, 20))
        font = QtGui.QFont()
        font.setWeight(75)
        font.setBold(True)
        self.label_17.setFont(font)
        self.label_17.setObjectName(_fromUtf8("label_17"))
        self.man_contrast = QtGui.QCheckBox(self.centralwidget)
        self.man_contrast.setGeometry(QtCore.QRect(20, 300, 131, 21))
        self.man_contrast.setObjectName(_fromUtf8("man_contrast"))
        self.label_18 = QtGui.QLabel(self.centralwidget)
        self.label_18.setGeometry(QtCore.QRect(70, 270, 62, 31))
        self.label_18.setObjectName(_fromUtf8("label_18"))
        self.vmin = QtGui.QSpinBox(self.centralwidget)
        self.vmin.setGeometry(QtCore.QRect(186, 300, 61, 25))
        self.vmin.setMaximum(1000000)
        self.vmin.setProperty(_fromUtf8("value"), 500)
        self.vmin.setObjectName(_fromUtf8("vmin"))
        self.vmax = QtGui.QSpinBox(self.centralwidget)
        self.vmax.setGeometry(QtCore.QRect(286, 300, 61, 25))
        self.vmax.setMaximum(1000000)
        self.vmax.setProperty(_fromUtf8("value"), 2000)
        self.vmax.setObjectName(_fromUtf8("vmax"))
        self.graphicsView = QtGui.QGraphicsView(self.centralwidget)
        self.graphicsView.setGeometry(QtCore.QRect(160, 300, 21, 21))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Light, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.AlternateBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Light, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.AlternateBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Light, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.AlternateBase, brush)
        self.graphicsView.setPalette(palette)
        self.graphicsView.setInteractive(False)
        self.graphicsView.setObjectName(_fromUtf8("graphicsView"))
        self.graphicsView_2 = QtGui.QGraphicsView(self.centralwidget)
        self.graphicsView_2.setGeometry(QtCore.QRect(260, 300, 21, 21))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Light, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.AlternateBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Light, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Inactive, QtGui.QPalette.AlternateBase, brush)
        brush = QtGui.QBrush(QtGui.QColor(69, 69, 69))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Light, brush)
        brush = QtGui.QBrush(QtGui.QColor(69, 69, 69))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Text, brush)
        brush = QtGui.QBrush(QtGui.QColor(106, 104, 100))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Base, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.Window, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Disabled, QtGui.QPalette.AlternateBase, brush)
        self.graphicsView_2.setPalette(palette)
        self.graphicsView_2.setInteractive(False)
        self.graphicsView_2.setObjectName(_fromUtf8("graphicsView_2"))
        self.label_9 = QtGui.QLabel(self.centralwidget)
        self.label_9.setGeometry(QtCore.QRect(690, 740, 131, 61))
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setWeight(75)
        font.setBold(True)
        self.label_9.setFont(font)
        self.label_9.setObjectName(_fromUtf8("label_9"))
        h5quicklook.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(h5quicklook)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 850, 22))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        h5quicklook.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(h5quicklook)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        h5quicklook.setStatusBar(self.statusbar)

        self.retranslateUi(h5quicklook)
        self.topplot.setCurrentIndex(0)
        self.bottomplot.setCurrentIndex(2)
        QtCore.QMetaObject.connectSlotsByName(h5quicklook)

    def retranslateUi(self, h5quicklook):
        h5quicklook.setWindowTitle(QtGui.QApplication.translate("h5quicklook", "ARCONS QuickLook", None, QtGui.QApplication.UnicodeUTF8))
        self.browse_button.setText(QtGui.QApplication.translate("h5quicklook", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("h5quicklook", "Loaded:", None, QtGui.QApplication.UnicodeUTF8))
        self.darksub_checkbox.setText(QtGui.QApplication.translate("h5quicklook", "Dark Subtraction?", None, QtGui.QApplication.UnicodeUTF8))
        self.flat_checkbox.setText(QtGui.QApplication.translate("h5quicklook", "Flat fielding?", None, QtGui.QApplication.UnicodeUTF8))
        self.displayimage_button.setText(QtGui.QApplication.translate("h5quicklook", "Display Image", None, QtGui.QApplication.UnicodeUTF8))
        self.saveimage_button.setText(QtGui.QApplication.translate("h5quicklook", "Save Image", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("h5quicklook", "Observation is", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("h5quicklook", "seconds.", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("h5quicklook", "Which chunk would you like to image?", None, QtGui.QApplication.UnicodeUTF8))
        self.label_6.setText(QtGui.QApplication.translate("h5quicklook", "to", None, QtGui.QApplication.UnicodeUTF8))
        self.label_7.setText(QtGui.QApplication.translate("h5quicklook", "seconds.", None, QtGui.QApplication.UnicodeUTF8))
        self.label_8.setText(QtGui.QApplication.translate("h5quicklook", "Do you want:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_11.setText(QtGui.QApplication.translate("h5quicklook", "Choose an observation:", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("h5quicklook", "with # of bins=", None, QtGui.QApplication.UnicodeUTF8))
        self.plot_pushButton.setText(QtGui.QApplication.translate("h5quicklook", "Plot", None, QtGui.QApplication.UnicodeUTF8))
        self.label_13.setText(QtGui.QApplication.translate("h5quicklook", "Header Info", None, QtGui.QApplication.UnicodeUTF8))
        self.time_up.setText(QtGui.QApplication.translate("h5quicklook", "t++", None, QtGui.QApplication.UnicodeUTF8))
        self.time_down.setText(QtGui.QApplication.translate("h5quicklook", "t--", None, QtGui.QApplication.UnicodeUTF8))
        self.satpercent.setToolTip(QtGui.QApplication.translate("h5quicklook", "Default is to have brightest 10% of pixels saturated \n"
"and scale grey color gradient over lower 90%.  \n"
"This prevents hot pixels from dominating color scaling.", None, QtGui.QApplication.UnicodeUTF8))
        self.label_14.setToolTip(QtGui.QApplication.translate("h5quicklook", "Default is to have brightest 10% of pixels saturated \n"
"and scale grey color gradient over lower 90%.  \n"
"This prevents hot pixels from dominating color scaling.", None, QtGui.QApplication.UnicodeUTF8))
        self.label_14.setText(QtGui.QApplication.translate("h5quicklook", "% Saturated", None, QtGui.QApplication.UnicodeUTF8))
        self.bad_pix.setText(QtGui.QApplication.translate("h5quicklook", "Mark bad pixels?", None, QtGui.QApplication.UnicodeUTF8))
        self.choosedark.setText(QtGui.QApplication.translate("h5quicklook", "Dark", None, QtGui.QApplication.UnicodeUTF8))
        self.drag_select_radioButton.setText(QtGui.QApplication.translate("h5quicklook", "Click && Drag", None, QtGui.QApplication.UnicodeUTF8))
        self.rect_select_radioButton.setText(QtGui.QApplication.translate("h5quicklook", "Single Click Rectangle", None, QtGui.QApplication.UnicodeUTF8))
        self.label_23.setText(QtGui.QApplication.translate("h5quicklook", "x", None, QtGui.QApplication.UnicodeUTF8))
        self.label_24.setText(QtGui.QApplication.translate("h5quicklook", "y", None, QtGui.QApplication.UnicodeUTF8))
        self.label_25.setText(QtGui.QApplication.translate("h5quicklook", "r", None, QtGui.QApplication.UnicodeUTF8))
        self.circ_select_radioButton.setText(QtGui.QApplication.translate("h5quicklook", "Single Click Circle", None, QtGui.QApplication.UnicodeUTF8))
        self.countlabel.setText(QtGui.QApplication.translate("h5quicklook", "Hist. Counts", None, QtGui.QApplication.UnicodeUTF8))
        self.choosesky.setText(QtGui.QApplication.translate("h5quicklook", "Sky", None, QtGui.QApplication.UnicodeUTF8))
        self.skysub_checkbox.setText(QtGui.QApplication.translate("h5quicklook", "Sky Subtraction?", None, QtGui.QApplication.UnicodeUTF8))
        self.pixelpath.setText(QtGui.QApplication.translate("h5quicklook", "Beammap Path", None, QtGui.QApplication.UnicodeUTF8))
        self.histogram_radio.setText(QtGui.QApplication.translate("h5quicklook", "Histogram", None, QtGui.QApplication.UnicodeUTF8))
        self.timestream_radio.setText(QtGui.QApplication.translate("h5quicklook", "Timestream", None, QtGui.QApplication.UnicodeUTF8))
        self.label_12.setText(QtGui.QApplication.translate("h5quicklook", "Select a pixel and plot:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_15.setText(QtGui.QApplication.translate("h5quicklook", "or", None, QtGui.QApplication.UnicodeUTF8))
        self.topplot.setItemText(0, QtGui.QApplication.translate("h5quicklook", "Peak Height", None, QtGui.QApplication.UnicodeUTF8))
        self.topplot.setItemText(1, QtGui.QApplication.translate("h5quicklook", "Parabola Fit", None, QtGui.QApplication.UnicodeUTF8))
        self.topplot.setItemText(2, QtGui.QApplication.translate("h5quicklook", "Baseline", None, QtGui.QApplication.UnicodeUTF8))
        self.bottomplot.setItemText(0, QtGui.QApplication.translate("h5quicklook", "Peak Height", None, QtGui.QApplication.UnicodeUTF8))
        self.bottomplot.setItemText(1, QtGui.QApplication.translate("h5quicklook", "Parabola Fit", None, QtGui.QApplication.UnicodeUTF8))
        self.bottomplot.setItemText(2, QtGui.QApplication.translate("h5quicklook", "Baseline", None, QtGui.QApplication.UnicodeUTF8))
        self.checkBox.setText(QtGui.QApplication.translate("h5quicklook", "Minus Baseline", None, QtGui.QApplication.UnicodeUTF8))
        self.checkBox_2.setText(QtGui.QApplication.translate("h5quicklook", "Minus Baseline", None, QtGui.QApplication.UnicodeUTF8))
        self.label_16.setText(QtGui.QApplication.translate("h5quicklook", "Plot 2", None, QtGui.QApplication.UnicodeUTF8))
        self.label_17.setText(QtGui.QApplication.translate("h5quicklook", "Plot 1", None, QtGui.QApplication.UnicodeUTF8))
        self.man_contrast.setText(QtGui.QApplication.translate("h5quicklook", "Manual Contrast", None, QtGui.QApplication.UnicodeUTF8))
        self.label_18.setText(QtGui.QApplication.translate("h5quicklook", "or:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_9.setText(QtGui.QApplication.translate("h5quicklook", "Histogram\n"
"Options", None, QtGui.QApplication.UnicodeUTF8))

from mpl_pyqt4_widget import MPL_Widget
