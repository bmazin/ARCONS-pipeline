#!/bin/python

import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
from util.confirm import confirm
from util.loadRawFile import ObsFile
from util import utils
from headers import ArconsHeaders

def applyWvlCal(obsFileName,photonListFileName,wvlCalFileName):
    h = 4.135668e-15 #eV s
    c = 2.998e8 #m/s
    angstromPerMeter = 1e10
    tickDuration = 1e-6 # s

    obsFile = ObsFile(obsFileName)
    wvlCalTable = loadWvlCalFile(wvlCalFileName,obsFile.nRow,obsFile.nCol)
    plTable = createPhotonListFile(photonListFileName)

    for iRow in xrange(obsFile.nRow):
        for iCol in xrange(obsFile.nCol):
            for iSec,sec in enumerate(obsFile.getPixel(iRow,iCol)):
                for photonPacket in sec:
                    timestamp,parabolaPeak,baseline = utils.parsePhotonPacket(photonPacket)
                    pulseHeight = parabolaPeak - baseline
                    xOffset = wvlCalTable[iRow,iCol,0]
                    amplitude = wvlCalTable[iRow,iCol,2]
                    yOffset = wvlCalTable[iRow,iCol,1]
                    energy = amplitude*(pulseHeight-xOffset)**2+yOffset
                    wavelength = h*c*angstromPerMeter/energy
                    print wavelength
                    if wavelength > 0 and wavelength != np.inf:
                        newRow = plTable.row
                        newRow['Xpix'] = iCol
                        newRow['Ypix'] = iRow
                        newRow['ArrivalTime'] = iSec+timestamp*tickDuration
                        newRow['Wavelength'] = wavelength
                        newRow.append()
    plTable.flush()


def createPhotonListFile(photonListFileName):
    scratchDir = os.getenv('INTERM_PATH','/')
    fullPhotonListFileName = os.path.join(scratchDir,photonListFileName)
    if (os.path.exists(fullPhotonListFileName)):
        if confirm('Photon list file  %s exists. Overwrite?'%fullPhotonListFileName,resp=False) == False:
            exit(0)
    zlibFilter = tables.Filters(complevel=1,complib='zlib',fletcher32=False)
    plFile = tables.openFile(fullPhotonListFileName,mode='w')
    plGroup = plFile.createGroup('/','photons','Group containing photon list')
    plTable = plFile.createTable(plGroup,'photons',ArconsHeaders.PhotonList,'Photon List Data',filters=zlibFilter)
    return plTable

def loadWvlCalFile(wvlCalFileName,nRow,nCol):
    scratchDir = os.getenv('INTERM_PATH','/')
    fullWvlCalFileName = os.path.join(scratchDir,wvlCalFileName)
    if (not os.path.exists(fullWvlCalFileName)):
        print 'wavelength cal file does not exist: ',fullWvlCalFileName
        exit(1)
    wvlCalFile = tables.openFile(fullWvlCalFileName,mode='r')
    wvlCalData = wvlCalFile.root.wavecal.calsoln
    nCalCoeffs = 3
    wvlCalTable = np.zeros([nRow,nCol,nCalCoeffs])
    for calPixel in wvlCalData:
        if calPixel['wave_flag'] == 0:
            wvlCalTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['polyfit']
    wvlCalFile.close()
    return wvlCalTable



if __name__ == '__main__':
    applyWvlCal('obs_20120919-131142.h5','pl_20120919-131142.h5','calsol_20120917-072537.h5')

    
 
