"""
Author:  Chris Stoughton   
Build up complete file names from:

MKDIR_DATA_DIR -- root of all raw
INTERM_DIR -- root of all generated files
run -- such as LICK201 or PAL2012
date -- in format yyyymmdd -- the locat year, month, and date of sunset
flavor -- obs or cal are the only ones I know about
tstamp -- in format yyyymmdd-hhmmss -- such as 20120920-123350

"""
import os
class FileName:
    
    def __init__(self, run, date, tstamp, \
                 mkidDataDir="/ScienceData", intermDir="/Scratch"):
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

                           
                            
