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
    
    def __init__(self, run, date, 
                 mkidDataDir="/ScienceData", intermDir="/Scratch"):
        self.mkidDataDir = mkidDataDir
        self.intermDir = intermDir
        self.run = run
        self.date = date

    def raw(self, flavor, tstamp):
        return self.mkidDataDir + os.sep + \
            self.run + os.sep + \
            self.date + os.sep + \
            flavor + "_" + tstamp + '.h5'

    def timeMask(self, tstamp):
        return self.intermDir + os.sep + \
            'timeMasks' + os.sep + \
            self.date + os.sep + \
            "timeMask_" + tstamp + '.h5'

    def calSoln(self, tstamp):
        return self.intermDir + os.sep + \
            'waveCalSolnFiles' + os.sep + \
            self.date + os.sep + \
            "calsol_" + tstamp + '.h5'

                           
                            
