#!/bin/python

import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt

class ObsFile:
    def __init__(self,fileName):
        self.loadFile(fileName)

    def __del__(self):
        self.file.close()

    def loadFile(self,fileName):
        self.fileName = fileName
        #make the full file name by joining the input name to the MKID_DATA_DIR (or . if the environment variable is not defined)
        dataDir = os.getenv('MKID_DATA_DIR','/')
        self.fullFileName = os.path.join(dataDir,self.fileName)
        if (not os.path.exists(self.fullFileName)):
            print 'file does not exist: ',self.fullFileName
            sys.exit(1)

        #open the hdf5 file
        self.file = tables.openFile(self.fullFileName,mode='r')

        #get the header
        header = self.file.root.header.header
        titles = header.colnames
        info = header[0] #header is a table with one row

        #get the beam image.
        try:
            self.beamImage = self.file.getNode('/beammap/beamimage')
        except:
            print 'Can\'t access beamimage'
            sys.exit(2)

        beamShape = self.beamImage.shape
        self.nRow = beamShape[0]
        self.nCol = beamShape[1]

    def __iter__(self):
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                pixelLabel = self.beamImage[iRow][iCol]
                pixelData = self.file.getNode('/'+pixelLabel)
                yield pixelData

    def getPixel(self,iRow,iCol):
        pixelLabel = self.beamImage[iRow][iCol]
        pixelData = self.file.getNode('/'+pixelLabel)
        return pixelData


