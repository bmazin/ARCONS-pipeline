"""
Author:  Chris Stoughton   
Build up complete file names from:

mkidDataDir -- root of all raw. If not specified, looks for system variable
                MKID_DATA_DIR, otherwise '/ScienceData'.
intermDir -- root of all generated files. If not specified, looks for 
                sys. variable INTERM_DIR, otherwise '/Scratch')
run -- such as LICK201 or PAL2012
date -- in format yyyymmdd -- the locat year, month, and date of sunset
flavor -- obs or cal are the only ones I know about
tstamp -- in format yyyymmdd-hhmmss -- such as 20120920-123350

NOTES
Updated 4/25/2013, JvE - if root directory names not provided, looks for system variables.

"""

import os
class FileName:
    
    def __init__(self, run='', date='', tstamp='', \
                 mkidDataDir=None, intermDir=None):
        
        if mkidDataDir is None:
            mkidDataDir = os.getenv('MKID_DATA_DIR', default="/ScienceData")
        if intermDir is None:
            intermDir = os.getenv('INTERM_DIR', default="/Scratch")
        
        self.mkidDataDir = mkidDataDir
        self.intermDir = intermDir
        self.run = run
        self.date = date
        self.tstamp = tstamp

    def makeName(self, prefix="plot_", suffix="png"):
        return prefix+self.tstamp+suffix

    def obs(self):
        return self.mkidDataDir + os.sep + \
            self.run + os.sep + \
            self.date + os.sep + \
            "obs_" + self.tstamp + '.h5'

    def cal(self):
        return self.mkidDataDir + os.sep + \
            self.run + os.sep + \
            self.date + os.sep + \
            "cal_" + self.tstamp + '.h5'

    def flat(self):
        return self.mkidDataDir + os.sep + \
            self.run + os.sep + \
            self.date + os.sep + \
            "flat_" + self.tstamp + '.h5'

    def timeMask(self):
        return self.intermDir + os.sep + \
            'timeMasks' + os.sep + \
            self.date + os.sep + \
            "timeMask_" + self.tstamp + '.h5'

    def calSoln(self):
        return self.intermDir + os.sep + \
            'waveCalSolnFiles' + os.sep + \
            self.date + os.sep + \
            "calsol_" + self.tstamp + '.h5'

    def calDriftInfo(self):
        return self.intermDir + os.sep + \
            'waveCalSolnFiles' + os.sep + \
            self.date + os.sep + \
            'drift_study'+ os.sep+\
            "calsol_" + self.tstamp + '_drift.h5'

    def flatSoln(self):
        return self.intermDir + os.sep + \
            'flatCalSolnFiles' + os.sep + \
            self.date + os.sep + \
            "flatsol_" + self.date + '.h5'

    def photonList(self):
        return self.intermDir + os.sep + \
            'photonLists' + os.sep + \
            self.date + os.sep + \
            "photons_" + self.tstamp + '.h5'
                           
    def packetMasterLog(self):
        return self.mkidDataDir + os.sep + \
            'PacketMasterLogs' + os.sep + \
            "obs_" + self.tstamp + '.log'

    def timeAdjustments(self):
        return self.intermDir + os.sep + \
            'timeAdjust' + os.sep + \
            self.run + '.h5'

                            
