#!/bin/python
'''
Author: Matt Strader        Date: March 6,2013

This program opens a series of observations of the crab pulsar.  It calculates the period, and calcluates the phase for every photon in the series of observations, It then plots a light curve, and a histogram of counts per period
'''
from flatcal.flatCal import FlatCal
from util.ObsFile import ObsFile
from util.FileName import FileName
from util.popup import PopUp
import matplotlib.pyplot as plt
import numpy as np
import datetime
import tables
import ephem
import matplotlib
import matplotlib.cm as cm
import os
import glob
import time
import linecache

class roachDelaysDescription(tables.IsDescription):
    """Description of roach delay structure for the output .h5 file."""
    obsFileName = tables.StringCol(40)  #To record associated obs file name
    roachDelays = tables.UInt8Col(8) #number of seconds that need to be added to timestamps from each roach for this obs

class firmwareDelayDescription(tables.IsDescription):
    """Description of firmware specific delay structure"""
    firmwareName = tables.StringCol(100)
    firmwareDelay = tables.Float64Col()

def main():

    run = 'LICK2012'
    year = '2012'
    initialPath = '/ScienceData'
    initialPath = os.path.join(initialPath,run)
    outPath = FileName(run=run).timeAdjustments()
    outFile = tables.openFile(outPath,'w')
    timeAdjustGroup = outFile.createGroup('/','timeAdjust','Times to add to timestamps')
    firmwareDelayTable = outFile.createTable(timeAdjustGroup,'firmwareDelay',firmwareDelayDescription,'Times to add to all timestamps taken with a firmware bof')
    newFirmwareEntry = firmwareDelayTable.row
    newFirmwareEntry['firmwareName']='chan_snap_v3_2012_Oct_30_1216.bof'
    newFirmwareEntry['firmwareDelay']=-41e-6 #s, subtract 41 us from timestamps
    newFirmwareEntry.append()
    firmwareDelayTable.flush()
    firmwareDelayTable.close()
    roachDelayTable = outFile.createTable(timeAdjustGroup,'roachDelays',roachDelaysDescription,'Times to add to each roach\'s timestamps')


    for sunsetDatePath in sorted(glob.glob(os.path.join(initialPath,year+'*'))):
        sunsetDate = os.path.basename(sunsetDatePath)
        for fullObsPath in sorted(glob.glob(os.path.join(sunsetDatePath,'obs*.h5'))):
            obsFileName = os.path.basename(fullObsPath)
            obsTStamp = obsFileName.split('.')[0].split('_')[1]
            print obsFileName
        

            obsFN = FileName(run=run,date=sunsetDate,tstamp=obsTStamp)
            pmLogFileName = obsFN.packetMasterLog()
            try:
                ob = ObsFile(fullObsPath)
            except:
                continue
                
            try:
                if os.path.getsize(pmLogFileName) <= 0:
                    continue
            except:
                print 'can\'t open Packet Master Log ',pmLogFileName
                continue

            try:
                f = open(pmLogFileName,'r')
            except:
                print 'can\'t open Packet Master Log ',pmLogFileName
                continue
            lastTstampLines = np.zeros(8)
            firstTstampLines = np.zeros(8)

            tstampLine = ''
            for line in f:
                if 'here' in line:
                    #skip over large log files with debug info
                    print 'skipping file with "here"'
                    continue
                if line.split(' ')[0] == 'bundle':
                    tstampLine = line
                    break

            if tstampLine == '':
                print 'skipping file without "bundle"'
                #didn't find lines with 'bundle' in them
                continue

            f.seek(0)
            for line in f:
                if line.split(' ')[0] == 'bundle':
                    try:
                        at = float(line.split('at')[1].split()[0])
                    except:
                        break
                    lastTstampLine = at
                    iRoach = int(line.split('roach')[1].split('took')[0].strip())
                    lastTstampLines[iRoach] = at
                    if firstTstampLines[iRoach] == 0:
                        firstTstampLines[iRoach] = at
                    
            packetReceivedUnixTimestamp = float((tstampLine.split('took')[1].split('total')[0].strip()))
            firstPacketDelay = packetReceivedUnixTimestamp-int(ob.getFromHeader('unixtime'))
                

            roachSecDelays =np.array(np.floor(lastTstampLines+firstPacketDelay-ob.getFromHeader('exptime')),dtype=np.int)
            if np.all(roachSecDelays >= 0):
                newEntry = roachDelayTable.row
                newEntry['obsFileName'] = os.path.basename(fullObsPath)
                newEntry['roachDelays'] = roachSecDelays
                newEntry.append()
                roachDelayTable.flush()
                print os.path.basename(fullObsPath),
                for i in range(8):
                    print roachSecDelays[i],
                print ''
    roachDelayTable.close()
    outFile.close()


if __name__=='__main__':
    main()

