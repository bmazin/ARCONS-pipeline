import os
import tables
import ephem
import math

class Overview:
    '''
    A simple class to generate summary information for ARCONS data 
    and calibration files. E.g. useage:
    
        from util import overview as ov
        ov.Overview('./').all()
        
        - lists info on all the ARCONS data files in the 
          current working directory.
    '''
    
    def __init__(self,directory):
        if os.path.exists(directory):
            self.dataDir = directory
        else:
            self.dataDir = os.getenv('MKID_DATA_DIR','.')
                
        self.calFiles = []
        self.obsFiles = []
        self.allFiles = []
        for (path,dirs,files) in os.walk(self.dataDir):
            for file in files:
                if file.endswith(".h5"):
                    fullFileName = os.path.join(path,file)
                    self.allFiles.append(fullFileName)
                    if file.startswith("cal"):
                        self.calFiles.append(fullFileName)
                    elif file.startswith("obs"):
                        self.obsFiles.append(fullFileName)
        self.obsFiles.sort()
        self.calFiles.sort()
        self.allFiles.sort(key=getSortKey)
    def obs(self):
        """Print summary information for observation files"""
        overviewFileList(self.obsFiles)
    def cal(self):
        """Print summary information for calibration files"""
        overviewFileList(self.calFiles)
    def all(self):
        """Print summary information for all files"""
        overviewFileList(self.allFiles)

def overviewFileList(fileList):
    for file in fileList:
        print getSummaryLine(file)

def getSortKey(file):
    return os.path.basename(file)[4:]

def getSummaryLine(fullFileName):
    fid = tables.openFile(fullFileName)
    header = fid.root.header.header
    titles = header.colnames
    info = header[0]
    raDeg = info[titles.index('ra')]
    filt = info[titles.index('filt')]
    raStr = str(ephem.hours(raDeg*math.pi/180.0))
    decDeg = info[titles.index('dec')]
    decStr = str(ephem.degrees(decDeg*math.pi/180.0))
    target = info[titles.index('target')]
    exptime = info[titles.index('exptime')]
    descrip = info[titles.index('description')]
    line = "%s%12s%12s%7s %6s %s " % (fullFileName, raStr, decStr, filt, exptime, target)
    if len(descrip) > 0:
        line += "\n       "+descrip
    fid.close()
    return line
