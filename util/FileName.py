"""
Author:  Chris Stoughton   
Build up complete file names from:

mkidDataDir -- root of all raw. If not specified, looks for system variable
                MKID_DATA_DIR, otherwise '/ScienceData'.
intermDir -- root of all generated files. If not specified, looks for 
                sys. variable INTERM_DIR, otherwise '/Scratch')
run -- such as LICK201 or PAL2012
date -- in format yyyymmdd -- the local year, month, and date of sunset
flavor -- obs or cal are the only ones I know about
tstamp -- in format yyyymmdd-hhmmss -- such as 20120920-123350

NOTES
Updated 4/25/2013, JvE - if root directory names not provided, looks for system variables.
5/22/2013 - can now optionally supply an ObsFile instance on calling instead of run, date, and tstamp.

"""

import os
class FileName:
    
    def __init__(self, run='', date='', tstamp='', \
                 mkidDataDir=None, intermDir=None, obsFile=None):
        '''
        To create a FileName instance, supply:
            run - string, name of folder which contains all obs files for this run
            date - string, date (usually of sunset), name of folder in which obs files for this date are stored
            tstamp - string, timestamp of the obs. file.
            obsFile - instead of run, date, tstamp, supply an obsFile instance instead, and it will pull out the required parameters automatically.
            mkidDataDir - raw data directory (uses path pointed to by system variable 'MKID_DATA_DIR' if not specified.)
            intermDir - data reduction product directory (uses path pointed to by system variable 'INTERM_DIR' if not specified.
        '''
            
        if mkidDataDir is None:
            mkidDataDir = os.getenv('MKID_DATA_DIR', default="/ScienceData")
        if intermDir is None:
            intermDir = os.getenv('INTERM_DIR', default="/Scratch")
        

        self.mkidDataDir = mkidDataDir
        self.intermDir = intermDir
        if obsFile is None:
            self.run = run
            self.date = date
            self.tstamp = tstamp
        else:
            # Split full file path of ObsFile instance into 
            # individual directory names

            # Warning:  Chris S. changed obsFile.fullFile name to obsFile
            # here and then again a few lines down.  July 11, 2013
            dirs = os.path.dirname(os.path.normpath(obsFile)).split(os.sep)
            #Pull out the relevant bits

            self.run = dirs[-2]
            self.date = dirs[-1]
            self.tstamp = (os.path.basename(obsFile)
                           .partition('obs_')[2]
                           .rpartition('.h5')[0])

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

    def cosmicMask(self):
        return self.intermDir + os.sep + \
            'cosmicMasks' + os.sep + \
            "cosmicMask_" + self.tstamp + '.h5'

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

    def centroidList(self):
        if self.tstamp == '' or self.tstamp == None:
            return self.intermDir + os.sep + \
                'centroidListFiles' + os.sep + \
                self.date + os.sep + \
                "centroid_" + self.date + '.h5'
        else:
            return self.intermDir + os.sep + \
                'centroidListFiles' + os.sep + \
                self.date + os.sep + \
                "centroid_" + self.tstamp + '.h5'

    def fluxSoln(self):
        if self.tstamp == '' or self.tstamp == None:
            return self.intermDir + os.sep + \
                'fluxCalSolnFiles' + os.sep + \
                self.date + os.sep + \
                "fluxsol_" + self.date + '.h5'
        else:
            return self.intermDir + os.sep + \
                'fluxCalSolnFiles' + os.sep + \
                self.date + os.sep + \
                "fluxsol_" + self.tstamp + '.h5'


    def flatSoln(self):
        if self.tstamp == '' or self.tstamp == None:
            return self.intermDir + os.sep + \
                'flatCalSolnFiles' + os.sep + \
                self.date + os.sep + \
                "flatsol_" + self.date + '.h5'
        else:
            return self.intermDir + os.sep + \
                'flatCalSolnFiles' + os.sep + \
                self.date + os.sep + \
                "flatsol_" + self.tstamp + '.h5'

    def flatInfo(self):
        if self.tstamp == '' or self.tstamp == None:
            return self.intermDir + os.sep + \
                'flatCalSolnFiles' + os.sep + \
                self.date + os.sep + \
                "flatsol_" + self.date + '.npz'
        else:
            return self.intermDir + os.sep + \
                'flatCalSolnFiles' + os.sep + \
                self.date + os.sep + \
                "flatsol_" + self.tstamp + '.npz'


    def oldFlatInfo(self):
        if self.tstamp == '' or self.tstamp == None:
            return self.intermDir + os.sep + \
                'oldFlatCalSolnFiles' + os.sep + \
                self.date + os.sep + \
                "flatsol_" + self.date + '.npz'
        else:
            return self.intermDir + os.sep + \
                'oldFlatCalSolnFiles' + os.sep + \
                self.date + os.sep + \
                "flatsol_" + self.tstamp + '.npz'

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

    
    ##################################
    
    def getComponents(self):
        '''
        Return a tuple of run, date, and timestamp.
        Potentially useful if an obsFile instance was passed on creation instead of individual components
        '''
        return self.run, self.date, self.tstamp
                            
