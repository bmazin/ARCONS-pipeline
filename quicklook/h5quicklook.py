#-----------------------------------
# h5quicklook_v10.py
#
# Written by Seth Meeker 07/16/11
#
# standalone program for generating images and pulse height histograms from h5 observation files
# inherits most of its functionality from the arcons GUI's quicklook imaging
# version 6 updated to subtract dark file
# version 7 updated to reimplement different aperture selections and plot histograms with dark counts subtracted
# version 9 updated to standardize orientation of images
# version 10 updated to display images from arbitrarily sized arrays
#
# KNOWN BUG: does not properly rescale histogram when multiple spectra are plotted.  Doesn't matter though
# since it makes no sense to plot multiple pixels while the spectra are uncalibrated.
#
# ----------------------------------

#import standard python libraries
import sys
import time
import struct
import os
from os.path import isfile
#import installed libraries
from matplotlib import pylab
from matplotlib import pyplot as plt
from numpy import *
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from tables import *
#import my functions
from h5quicklook_v8_gui import Ui_h5quicklook

c = 3.0E17 #[nm/s]
h = 4.13567E-15 #[ev*s]

class StartQt4(QMainWindow):
	def __init__(self,parent=None):
		QWidget.__init__(self, parent)
		self.ui = Ui_h5quicklook()
		self.ui.setupUi(self)
		
		#set some important permanent parameters (array size, energy range, directories)
		#self.nxpix = 32
		#self.nypix = 32
		#self.total_pix = self.nxpix*self.nypix
		
		#self.imagex = self.ui.tv_image.width()
		#self.imagey = self.ui.tv_image.height()
		
		self.imagex = 308
		self.imagey = 322

		self.position = self.ui.tv_image.pos()

		self.Emin = 0.92 #eV corresponds to 1350 nm
		self.Emax = 3.18 #eV corresponds to 390 nm (changed from original 4.2 ev)
		self.bintype = "wavelength"
		
		self.pix_select_mode = "rect"
		
		#store list of pixels we want in our histogram
		self.histogram_pixel = []
		self.plot_type = "bars" #change to "points" for histogram plotted as point scatter, anything else for bars
		
		#read in bad pixels from file so they can be "X"d out in the quicklook
		self.bad_pix = []
		
		self.flatfile = ""
		self.darkfile = ""
		self.skyfile = ""

		#use mouse to select pixel from tv_image, also triggers a new spectrum to be displayed
		self.ui.tv_image.mousePressEvent = self.start_pixel_select
		self.ui.tv_image.mouseReleaseEvent = self.end_pixel_select
		#signal from image thread to make spectrum
		#QObject.connect(self.image_thread,SIGNAL("new_spectrum"),self.display_spectra)
		
		self.ui.histogram_plot.canvas.xtitle = "Pulse Heights"
		self.ui.histogram_plot.canvas.format_labels()
		QObject.connect(self.ui.plot_pushButton, SIGNAL("clicked()"), self.plot_histogram)
		
		#button signals
		QObject.connect(self.ui.mode_buttonGroup, SIGNAL("buttonClicked(int)"), self.set_pix_select_mode)
		QObject.connect(self.ui.browse_button, SIGNAL("clicked()"), self.file_dialog)
		QObject.connect(self.ui.displayimage_button, SIGNAL("clicked()"), self.make_image)
		QObject.connect(self.ui.saveimage_button,SIGNAL("clicked()"), self.save_image)
		QObject.connect(self.ui.time_up, SIGNAL("clicked()"), self.time_up)
		QObject.connect(self.ui.time_down, SIGNAL("clicked()"), self.time_down)
		QObject.connect(self.ui.choosedark,SIGNAL("clicked()"), self.choosedark)
		QObject.connect(self.ui.choosesky,SIGNAL("clicked()"), self.choosesky)
		#QObject.connect(self.ui.satpercent, SIGNAL("valueChanged(double)"), self.make_image)
		
	#open dialog to select data path
	def file_dialog(self):
		self.newdatafile = QFileDialog.getOpenFileName(parent=None, caption=QString(str("Choose Obs File")),directory = ".",filter=QString(str("H5 (*.h5)")))
		if len(self.newdatafile)!=0:
			self.datafile = self.newdatafile
			splitname = self.datafile.split("/")
			justname = splitname[-1]
			self.ui.filename_label.setText(str(justname))
			self.blankscene = QGraphicsScene()
			self.blankscene.clear()
			self.ui.tv_image.setScene(self.blankscene)
			self.ui.tv_image.show()
			self.ui.histogram_plot.canvas.ax.clear()
			self.ui.histogram_plot.canvas.format_labels()
			self.ui.histogram_plot.canvas.draw()
			self.display_obs_time()
			self.display_header()
		
	def choosedark(self):
		self.darkfile = QFileDialog.getOpenFileName(parent=None, caption=QString(str("Choose Dark File")),directory = ".",filter=QString(str("H5 (*.h5)")))
		self.darkfile = str(self.darkfile)
		print "loading darkfile",str(self.darkfile)
		
	def choosesky(self):
		self.skyfile = QFileDialog.getOpenFileName(parent=None, caption=QString(str("Choose Sky File")),directory = ".",filter=QString(str("H5 (*.h5)")))
		self.skyfile = str(self.skyfile)
		print "loading skyfile",str(self.skyfile)
		
	def save_image(self):
		self.savefile = QFileDialog.getSaveFileName(parent=None, caption=QString(str("Save File")), directory = ".")
		os.rename(self.imfile, self.savefile)		
		
	def display_obs_time(self):
		
		try:
			h5file = openFile(str(self.datafile), mode = "r")
			htable = h5file.root.header.header.read()
			obstime = htable["exptime"]
			self.ui.obstime_lcd.display(obstime)
			self.ui.endtime_spinbox.setValue(obstime)
			h5file.close()
		except:
			print "Unable to load Header, checking beammap"
			h5file = openFile(str(self.datafile), mode = "r")
			bmap = h5file.root.beammap.beamimage.read()
			#bmap = rot90(bmap,2)
			h5file.close()
			
			self.nxpix = shape(bmap)[1]
			self.nypix = shape(bmap)[0]
			
			self.scalex=float(self.imagex/float(self.nxpix))
			self.scaley=float(self.imagey/float(self.nypix))
			#if self.nxpix > self.nypix:
			#	self.scaley = self.scaley * float(self.nypix/float(self.nxpix))
			#else:
			#	self.scalex = self.scalex * float(self.nxpix/float(self.nypix))
			
			self.bad_pix = []
			tempbmap = reshape(bmap,self.nxpix*self.nypix)
			for i in range(self.nxpix*self.nypix):
				if len(concatenate(h5file.root._f_getChild(tempbmap[i])[:])) == 0:
					self.bad_pix.append(i)
					
			for i in range(self.nxpix):
				for j in range(self.nypix):
					try:
						photons= h5file.root._f_getChild(bmap[i][j]).read()
						obstime = len(photons)
						self.ui.obstime_lcd.display(obstime)
						self.ui.endtime_spinbox.setValue(obstime)
						return
					except NoSuchNodeError:
						continue
			print "unable to find any pixels with data. Check beammap is present"

	def display_header(self):
		h5file = openFile(str(self.datafile), mode = "r")
		try:
			header = h5file.root.header.header
			self.ui.header_info.clear()	
			#titles = ['LOfreq', 'airmass', 'alt', 'az', 'beammap', 'calfile', 'datadir', 'dec', 'epoch', 'equinox', 'exptime', 'focus', 'inatten', 'instrument', 'jd', 'localtime', 'lst', 'obsalt', 'obslat', 'obslong', 'outatten','parallactic', 'platescl', 'ra', 'target', 'telescope', 'timezone', 'ut']
			titles = header.colnames
			info = header[0]
			for i in range(len(titles)):
				self.ui.header_info.append(titles[i] + ":\n" + str(info[i]) + "\n")
		except NoSuchNodeError:
			self.ui.header_info.clear()
			self.ui.header_info.append('No header info')
		h5file.close()
	
	def set_dark_sub(self):
		#Use radio button to set if we want sky subtraction to be on or off
		if self.ui.darksub_checkbox.isChecked():
			self.dark_subtraction = True
		else:
			self.dark_subtraction = False
			
	def set_sky_sub(self):
		#Use radio button to set if we want sky subtraction to be on or off
		if self.ui.skysub_checkbox.isChecked():
			self.sky_subtraction = True
		else:
			self.sky_subtraction = False
			
	def set_flat_sub(self):
		pass
			
	def set_image_color(self):
		if self.ui.color_checkbox.isChecked():
			self.color = True
		else:
			self.color = False
			
	def unpack_file(self, file):
		f = open(file, 'rb')
		inputdata = f.read()
		f.close()
		numbers = struct.unpack((self.nxpix*self.nypix)*'10I', inputdata)
		darray = reshape(numbers,((self.nxpix*self.nypix),10))
		return darray		
	
	def make_image(self):
		self.set_dark_sub()
		#self.set_image_color()
		self.set_flat_sub()
		self.set_sky_sub()
		
		if self.ui.darksub_checkbox.isChecked():
			self.dark_subtraction = True
		if self.ui.skysub_checkbox.isChecked():
			self.sky_subtraction = True
		
		nphot=0
		
		ti = self.ui.starttime_spinbox.value()
		tf = self.ui.endtime_spinbox.value()
		
		if ti>tf:
			copyti = ti
			ti=tf
			tf=copyti
			print "WARNING: selected ti > tf. They were switched for you."
		
		h5file = openFile(str(self.datafile), 'r')
		bmap = h5file.root.beammap.beamimage.read()
		#bmap = rot90(bmap,2)
		self.nxpix = shape(bmap)[1]
		self.nypix = shape(bmap)[0]

		self.scalex=float(self.imagex/float(self.nxpix))
		self.scaley=float(self.imagey/float(self.nypix))

		#if self.nxpix > self.nypix:
		#	self.scaley = self.scaley * float(self.nypix/float(self.nxpix))
		#else:
		#	self.scalex = self.scalex * float(self.nxpix/float(self.nypix))
		
		all_photons = []
		for j in range(self.nxpix*self.nypix):
			all_photons.append([])
		
		
		if self.dark_subtraction == True:
			darkrate = zeros((self.nypix,self.nxpix))
			darkh5 = openFile(str(self.darkfile), 'r')
			darkbmap = darkh5.root.beammap.beamimage.read()
			#darkbmap = rot90(darkbmap,2)
			
		if self.sky_subtraction == True:
			skyrate = zeros((self.nypix,self.nxpix))
			skyh5 = openFile(str(self.skyfile), 'r')
			skybmap = skyh5.root.beammap.beamimage.read()
			#skybmap = rot90(skybmap,2)

		totalhist = zeros((self.ui.nbins.value()))

		h5file = openFile(str(self.datafile), 'r')
		bmap = h5file.root.beammap.beamimage.read()
		#bmap = rot90(bmap,2)
		
		#print "(0,0)", str(bmap[0][0])
		#print "(0,31)", str(bmap[0][31])
		#print "(31,0)", str(bmap[31][0])
		#print "(31,31)", str(bmap[31][31])
		
		counts = zeros((self.nypix,self.nxpix))
		
		bins = range(self.ui.nbins.value()+1)
		
		for i in xrange(self.nypix):
			for j in xrange(self.nxpix):
				#if i*32+j in self.histogram_pixel:
					if bmap[i][j] == '':
						counts[i][j]=0
						subtracted = zeros((self.ui.nbins.value()))
						continue
					try:
						#etime = len(h5file.root._f_getChild(bmap[i][j]))
						photons= concatenate(h5file.root._f_getChild(bmap[i][j])[ti:tf])
						#obsheights = right_shift(photons,32)%4096
						parabheights= right_shift(photons,44)%4096
						npparabheights = array(parabheights, dtype=float)
						baseline= right_shift(photons,20)%4096
						npbaseline = array(baseline, dtype=float)
						obsheights = npbaseline-npparabheights
						#obsheights = baseline						
						#obsheights = parabheights

						if len(obsheights)==0:
							counts[i][j]=0
							continue
						else:
							obshist,bins = histogram(obsheights,bins=self.ui.nbins.value(),range=(obsheights.min(),obsheights.max()))
						#print bins
						
						#print bins
						if self.dark_subtraction == True:
							dtime = float(len(darkh5.root._f_getChild(darkbmap[i][j])))
							photons= concatenate(darkh5.root._f_getChild(darkbmap[i][j])[:])
							darkheights= right_shift(photons,32)%4096
							
							if len(darkheights)==0:
								subtracted = obshist
								darkhist = zeros((self.ui.nbins.value()))
							else:
								darkhist,bins = histogram(darkheights,bins=self.ui.nbins.value(),range=(obsheights.min(),obsheights.max()))
								subtracted = obshist-(darkhist*(tf-ti)/float(dtime))
								#print "dark subtracted..."
						else:
							subtracted = obshist
						
						for p in range(len(subtracted)):
							if subtracted[p]<0:
								subtracted[p]=0
						
						if self.sky_subtraction == True:
							if self.dark_subtraction != True:
								stime = float(len(skyh5.root._f_getChild(skybmap[i][j])))
								photons= concatenate(skyh5.root._f_getChild(skybmap[i][j])[:])
								skyheights= right_shift(photons,32)%4096
								
								if len(skyheights)==0:
									pass
								else:
									skyhist,bins = histogram(skyheights,bins=self.ui.nbins.value(),range=(obsheights.min(),obsheights.max()))
									skysubtracted = skyhist
									#print skysubtracted
									
									for p in range(len(skysubtracted)):
										if skysubtracted[p] <0:
											skysubtracted[p]=0
									
									subtracted -= (skysubtracted*(tf-ti)/float(stime))
									#print "sky subtracted..."
							else:
								stime = float(len(skyh5.root._f_getChild(skybmap[i][j])))
								photons= concatenate(skyh5.root._f_getChild(skybmap[i][j])[:])
								skyheights= right_shift(photons,32)%4096
								
								if len(skyheights)==0:
									pass
								else:
									skyhist,bins = histogram(skyheights,bins=self.ui.nbins.value(),range=(obsheights.min(),obsheights.max()))
									skysubtracted = skyhist-(darkhist*(stime)/float(dtime))
									#print skysubtracted
									
									for p in range(len(skysubtracted)):
										if skysubtracted[p] <0:
											skysubtracted[p]=0
									
									subtracted -= (skysubtracted*(tf-ti)/float(stime))
									#print "sky subtracted..."
						else:
							pass
						
						for p in range(len(subtracted)):
							if subtracted[p]<0:
								subtracted[p]=0
						counts[i][j] = sum(subtracted)
						
					except NoSuchNodeError:
						counts[i][j]=0
						subtracted = zeros((self.ui.nbins.value()))
					
					for p in range(len(subtracted)):
						if subtracted[p]<0:
							subtracted[p]=0
						
					if counts[i][j]<0:
						counts[i][j]=0
			
					#totalhist += subtracted
					nphot += counts[i][j]
	
		photon_count = counts
		im = photon_count
		
		#histogram equalization on image
		#nbr_bins = 256
		#imhist,bins = histogram(im, nbr_bins, normed=True)
		#cdf = imhist.cumsum() #cumulative dist func
		#cdf = 255*cdf/cdf[-1] #normalize
		#photon_count = interp(im, bins[:-1],cdf)
		
		#photon_count = im.reshape(32,32)
		
		#if self.ui.darksub_checkbox.isChecked() == True and self.darkfile!="":
			#print "subtracting dark rate"
			#photon_count -= ((self.tf-self.ti)*(self.darkrate))
		
		#for n in range(32):
			#for m in range(32):
				#if photon_count[n][m] < 0:
					#photon_count[n][m] = 0
		
		photon_count = flipud(photon_count)
		
		indices = sort(reshape(photon_count,self.nxpix*self.nypix))
		#print indices
		brightest = int(self.ui.satpercent.value()*(self.nxpix*self.nypix))
		#print brightest
		self.vmax = indices[(-1*brightest)]
		print self.vmax
		fig = plt.figure(figsize=(0.01*self.nxpix,0.01*self.nypix), dpi=100, frameon=False)
		im = plt.figimage(photon_count, cmap='gray', vmax = self.vmax)
			
		self.imfile = "TV_frame.png"
		plt.savefig(self.imfile, pad_inches=0)
		print "done making image."
			
		self.display_image()
			
		if len(self.histogram_pixel) != 0:
			self.plot_histogram()
		
	def time_up(self):
		ti = self.ui.starttime_spinbox.value()
		tf = self.ui.endtime_spinbox.value()
		ti += 1
		tf += 1
		self.ui.starttime_spinbox.setValue(ti)
		self.ui.endtime_spinbox.setValue(tf)
		self.make_image()
		
	def time_down(self):
		ti = self.ui.starttime_spinbox.value()
		tf = self.ui.endtime_spinbox.value()
		ti -= 1
		tf -= 1
		self.ui.starttime_spinbox.setValue(ti)
		self.ui.endtime_spinbox.setValue(tf)
		self.make_image()
		
	def makepixmap(self, imagefile, scalex=1, scaley=1):
		'''Given an image file this function converts them to pixmaps to be displayed by QT gui'''
		qtimage = QImage(imagefile)
		width, height = qtimage.size().width(), qtimage.size().height()
		qtimage = qtimage.scaled(width*scalex,height*scaley)
		pix = QPixmap.fromImage(qtimage)
		return pix

	def display_image(self):
		#search directory for image
		self.imagefile = "./TV_frame.png"
		
		self.ui.tv_image.setGeometry(self.position.x(), self.position.y()-8, self.scalex*(self.nxpix)+4, self.scaley*(self.nypix)+4)
		#print self.nxpix
		#print self.nypix
		#convert to pixmap
		if isfile(self.imagefile):
			pix = self.makepixmap(self.imagefile, scalex=self.scalex, scaley=self.scaley)
			#display pixmap
			self.scene = QGraphicsScene()
			self.scene.addPixmap(pix)
			
			### BAD PIXELS AND SPECTRA SELECTION NOT IMPLEMENTED ###
			if self.ui.bad_pix.isChecked() == True:
				for bp in self.bad_pix:
					x = bp%(self.nxpix)
					y = (self.nypix-1) - bp/(self.nxpix)
					self.scene.addLine(self.scalex*x,self.scaley*y,self.scale*x+(self.scalex),self.scaley*y+(self.scaley),Qt.red)
					self.scene.addLine(self.scalex*x,self.scaley*y+(self.scaley),self.scalex*x+(self.scalex),self.scaley*y,Qt.red)
			if self.histogram_pixel != []:
				for p in self.histogram_pixel:
					x = p%(self.nxpix)
					y = (self.nypix-1) - p/self.nxpix					
					self.scene.addRect(self.scalex*(x),self.scaley*(y),(self.scalex),(self.scaley), Qt.blue)
					
			self.ui.tv_image.setScene(self.scene)
			self.ui.tv_image.show()
			#os.remove(str(imagefile)) #want to keep this around for saving purposes
		else: 
			self.blankscene = QGraphicsScene()
			self.blankscene.clear()
			self.ui.tv_image.setScene(self.blankscene)
			self.ui.tv_image.show()
	
	#def set_plot_mode(self):
		#change between spectra plots and signal to noise ratio plots
		#if self.ui.plot_snr_radioButton.isChecked():
			#self.image_thread.set_plot_mode("snr")
			#self.ui.spectra_plot.canvas.ytitle = "Signal to Noise"
		#else:
			#self.image_thread.set_plot_mode("spectra")
			#self.ui.spectra_plot.canvas.ytitle = "Counts"
		#pass
	
	def set_pix_select_mode(self):
		if self.ui.drag_select_radioButton.isChecked():
			self.pix_select_mode = "drag"
		elif self.ui.rect_select_radioButton.isChecked():
			self.pix_select_mode = "rect"
		elif self.ui.circ_select_radioButton.isChecked():
			self.pix_select_mode = "circ"
			
	def start_pixel_select(self,event):
		#Mouse press returns x,y position of first pixel to be used in spectra
		self.startrawx,self.startrawy = event.pos().x(), event.pos().y()
		self.startpx = int(self.startrawx/self.scalex)
		self.startpy = int((self.nypix) - self.startrawy/self.scaley)
		self.startpix = self.nxpix*self.startpy+self.startpx
	
	def end_pixel_select(self,event):
		#Mouse release returns x,y position of last pixel to be used in spectra
		self.endrawx,self.endrawy = event.pos().x(), event.pos().y()
		self.endpx = int(self.endrawx/self.scalex)
		self.endpy = int((self.nypix) - self.endrawy/self.scaley)
		self.endpix = self.nxpix*self.endpy+self.endpx
		self.pixel_list()
		
	def pixel_list(self):
		#if click and drag selection is on, add new pixels to the list of all pixels being plotted
		if self.pix_select_mode == "drag":
			if self.startpix != self.endpix:
				#get all pix in box
				allx = range(min(self.startpx,self.endpx),max(self.startpx, self.endpx)+1)
				ally = range(min(self.startpy,self.endpy),max(self.startpy, self.endpy)+1)
				pix = []
				for x in allx:
					for y in ally:
						pix.append(y*self.nxpix+x)
						#pix.append((31-y)*32+x)
			else:
				pix = [self.startpix]
		elif self.pix_select_mode == "rect":
			#get all pix in box
			length = self.ui.rect_x_spinBox.value()
			height = self.ui.rect_y_spinBox.value()
			allx = range(self.startpx-int(ceil(length/2.0)-1),self.startpx+int(floor(length/2.0))+1)
			ally = range(self.startpy-int(ceil(height/2.0)-1),self.startpy+int(floor(height/2.0))+1)
			self.histogram_pixel = [] #clear spectrum array
			pix=[]
			#self.histogram_pixel.append(self.startpix)
			for x in allx:
				for y in ally:
					pix.append(y*self.nxpix+x)
					#pix.append((31-y)*32+x)
		elif self.pix_select_mode == "circ":
			r = self.ui.circ_r_spinBox.value()
			length = 2*r
			height = length
			allx = range(self.startpx-int(ceil(length/2.0)),self.startpx+int(floor(length/2.0))+1)
			ally = range(self.startpy-int(ceil(height/2.0)),self.startpy+int(floor(height/2.0))+1)
			self.histogram_pixel = [] #clear spectrum array
			pix = []
			for x in allx:
				for y in ally:
					if (abs(x-self.startpx))**2+(abs(y-self.startpy))**2 <= (r)**2:
						pix.append(y*self.nxpix+x)
						#pix.append((31-y)*32+x)
		for element in pix:
				#check for repeated pixels (clicked for deselection) and out of bounds pixels, remove from total array
				if element not in self.histogram_pixel:
					if element >= 0 and element <= (self.nxpix*self.nypix-1):# and element not in self.bad_pix:
						self.histogram_pixel.append(element)
				else:
					spot = self.histogram_pixel.index(element)
					del(self.histogram_pixel[spot])
		#if self.observing == False:
			#self.image_thread.update_spectrum(self.bindir)
		self.display_image()
		self.plot_histogram()
	
	def plot_histogram(self):		
		self.set_dark_sub()
		self.set_sky_sub()
		self.set_flat_sub()
		
		if self.histogram_pixel != []:
			
			nphot=0
			
			ti = self.ui.starttime_spinbox.value()
			tf = self.ui.endtime_spinbox.value()
			
			if ti>tf:
				copyti = ti
				ti=tf
				tf=copyti
				print "WARNING: selected ti > tf. They were switched for you."
			
			h5file = openFile(str(self.datafile), 'r')
			bmap = h5file.root.beammap.beamimage.read()
			#bmap = rot90(bmap,2)
			
			all_photons = []
			for j in range(self.nxpix*self.nypix):
				all_photons.append([])
			
			if self.dark_subtraction == True:
				darkrate = zeros((self.nypix,self.nxpix))

				darkh5 = openFile(str(self.darkfile), 'r')
				darkbmap = darkh5.root.beammap.beamimage.read()
				#darkbmap = rot90(darkbmap,2)
				
			if self.sky_subtraction == True:
				skyrate = zeros((self.nypix,self.nxpix))

				skyh5 = openFile(str(self.skyfile), 'r')
				skybmap = skyh5.root.beammap.beamimage.read()
				#skybmap = rot90(skybmap,2)

			totalhist = zeros((self.ui.nbins.value()))
	
			counts = zeros((self.nypix,self.nxpix))
			
			bins = range(self.ui.nbins.value()+1)
			
			for i in xrange(self.nypix):
				for j in xrange(self.nxpix):
					if i*self.nxpix+j in self.histogram_pixel:
						
						self.ui.pixelpath.setText(str(bmap[i][j]))
						
						if bmap[i][j] == '':
							counts[i][j]=0
							subtracted = zeros((self.ui.nbins.value()))
							continue
						try:
							#etime = len(h5file.root._f_getChild(bmap[i][j]))
							photons= concatenate(h5file.root._f_getChild(bmap[i][j])[ti:tf])
							#obsheights= right_shift(photons,32)%4096
							parabheights= right_shift(photons,44)%4096
							npparabheights = array(parabheights, dtype=float)
							baseline= right_shift(photons,20)%4096
							npbaseline = array(baseline, dtype=float)
							obsheights = npbaseline-npparabheights
							#obsheights = baseline
							#obsheights = parabheights

							if len(obsheights)==0:
								counts[i][j]=0
								continue
							else:
								obshist,bins = histogram(obsheights,bins=self.ui.nbins.value(),range=(obsheights.min(),obsheights.max()))
							#print bins
							if self.dark_subtraction == True:
								dtime = float(len(darkh5.root._f_getChild(darkbmap[i][j])))
								photons= concatenate(darkh5.root._f_getChild(darkbmap[i][j])[:])
								darkheights= right_shift(photons,32)%4096
							
								if len(darkheights)==0:
									subtracted = obshist
									darkhist = zeros((self.ui.nbins.value()))
								else:
									darkhist,bins = histogram(darkheights,bins=self.ui.nbins.value(),range=(obsheights.min(),obsheights.max()))
									subtracted = obshist-(darkhist*(tf-ti)/float(dtime))
							else:
								subtracted=obshist
								
							for p in range(len(subtracted)):
								if subtracted[p]<0:
									subtracted[p]=0
								
							if self.sky_subtraction == True:
								if self.dark_subtraction != True:
									stime = float(len(skyh5.root._f_getChild(skybmap[i][j])))
									photons = concatenate(skyh5.root._f_getChild(skybmap[i][j])[:])
									skyheights = right_shift(photons,32)%4096
									
									if len(skyheights)==0:
										pass
									else:
										skyhist,bins = histogram(skyheights,bins=self.ui.nbins.value(),range=(obsheights.min(),obsheights.max()))
										skysubtracted = skyhist
										
										for p in range(len(skysubtracted)):
											if skysubtracted[p] <0:
												skysubtracted[p]=0
										
										subtracted = subtracted-(skysubtracted*(tf-ti)/float(stime))
								else:
									stime = float(len(skyh5.root._f_getChild(skybmap[i][j])))
									photons = concatenate(skyh5.root._f_getChild(skybmap[i][j])[:])
									skyheights = right_shift(photons,32)%4096
									
									if len(skyheights)==0:
										pass
									else:
										skyhist,bins = histogram(skyheights,bins=self.ui.nbins.value(),range=(obsheights.min(),obsheights.max()))
										skysubtracted = skyhist-(darkhist*(stime)/float(dtime))
										
										for p in range(len(skysubtracted)):
											if skysubtracted[p] <0:
												skysubtracted[p]=0
										
										subtracted = subtracted-(skysubtracted*(tf-ti)/float(stime))
							else:
								pass
							
							counts[i][j] = sum(subtracted)
							if counts[i][j]<0:
								counts[i][j]=0
				
						except NoSuchNodeError:
							counts[i][j]=0
							subtracted = zeros((self.ui.nbins.value()))
						
						for p in range(len(subtracted)):
							if subtracted[p]<0:
								subtracted[p]=0
						
						totalhist += subtracted
						nphot += counts[i][j]
			
			print "plotting histogram of ", nphot, " pulse heights"
			self.ui.countlabel.setText(str(nphot))
			nbins = self.ui.nbins.value()
			self.ui.histogram_plot.canvas.ax.clear()
			#if self.plot_type == "point":
				#photon_hist = histogram(new_photons, bins = nbins, range = (min(new_photons), max(new_photons)))
				#self.ui.histogram_plot.canvas.ax.plot(photon_hist[1][1:],photon_hist[0],'o')
			#else:
				#self.ui.histogram_plot.canvas.ax.hist(new_photons, bins = nbins, range = (min(new_photons),max(new_photons)), histtype='bar')
			self.ui.histogram_plot.canvas.ax.bar(bins[:-1],totalhist, width=(bins[1]-bins[0]), bottom=0)
			
			#self.ui.histogram_plot.canvas.format_labels()
			self.ui.histogram_plot.canvas.draw()
			print "done"
		#except ValueError:
			#self.ui.histogram_plot.canvas.ax.clear()
			#self.ui.histogram_plot.canvas.format_labels()
			#self.ui.histogram_plot.canvas.draw()
		#except NoSuchNodeError:
			#print "file has no data"
	
	#def display_spectra(self, pixel, bins, counts, SNR):
		#plots spectra directly from bin file data passed from image thread
		#self.ui.spectra_plot.canvas.ax.clear()
		#self.ui.pixel_no_lcd.display(pixel)
		#self.ui.integrated_snr_lcd.display(SNR)
		#self.ui.spectra_plot.canvas.ax.plot(bins,counts,'o')
		#self.ui.spectra_plot.canvas.format_labels()
		#self.ui.spectra_plot.canvas.draw()
		
	#Data generation in dataSim.py and organization in dataBin.py, and image generation is in own thread that stays in this gui application.	
	#image_Worker retrieves and opens organized files, saves to local arrays, storing as cumulative arrays 
	#if we are in an observation, and makes quicklook images and spectra from these local arrays.
	'''
	def subtract_sky(self, meds):
		#Sky subtraction works by taking the median of counts over all pixels in every energy slice
		#and subtracting that number from each pixel's count total.  Since there are generally 8 to 9 times more
		#sky pixels than object pixels this median should come from a sky pixel and can be used as 
		#a good representation of sky counts at that energy
		#subtract each bins median values from each pixel in that bin
		self.C0[:] = [x-int(meds[0]) for x in self.C0]
		self.C1[:] = [x-int(meds[1]) for x in self.C1]
		self.C2[:] = [x-int(meds[2]) for x in self.C2]
		self.C3[:] = [x-int(meds[3]) for x in self.C3]
		self.C4[:] = [x-int(meds[4]) for x in self.C4]
		self.C5[:] = [x-int(meds[5]) for x in self.C5]
		self.C6[:] = [x-int(meds[6]) for x in self.C6]
		self.C7[:] = [x-int(meds[7]) for x in self.C7]
		self.C8[:] = [x-int(meds[8]) for x in self.C8]
		self.C9[:] = [x-int(meds[9]) for x in self.C9]
		'''
	def calc_mean_energy(self):
		for p in range(self.nxpix*self.nypix):
			self.me[p] = ((self.C0[p]*self.E0 + self.C1[p]*self.E1 + self.C2[p]*self.E2 + self.C3[p]*self.E3 \
						+ self.C4[p]*self.E4 + self.C5[p]*self.E5 + self.C6[p]*self.E6 + self.C7[p]*self.E7 \
						+ self.C8[p]*self.E8 + self.C9[p]*self.E9) /self.pc[p])
		if self.bintype == "wavelength":
			self.me[:] = [h*c/e for e in self.me]
			'''
	def calculate_SNR(self, totalcounts, medians, npix):
		total_signal = 0
		total_n = 0
		self.integrated_SNR = 0
		self.SNR = [0,0,0,0,0,0,0,0,0,0]
		for i in range(len(medians)):
			if self.sky_subtraction == True:
				signal = totalcounts[i]
				noise = npix*medians[i]
			else:
				signal = totalcounts[i]-npix*medians[i]
				noise = npix*medians[i]
			if signal < 0:
				signal = 0
			if noise == 0:
				noise = 1 #pretty much only happens when a bin has no photons in it at all
			total_signal += signal
			total_n += noise
			self.SNR[i] = signal/(sqrt(noise))
		self.integrated_SNR = total_signal/(sqrt(total_n))
	'''

	def closeEvent(self, event=None):
		if isfile(self.imagefile):
			os.remove(str(self.imagefile))
				
if __name__ == "__main__":
	app = QApplication(sys.argv)
	myapp = StartQt4()
	myapp.show()
	app.exec_()
