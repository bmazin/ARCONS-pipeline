'''
Author: Matt Strader        Date: January 08, 2015
'''

import numpy as np
import os
from headers.CalHeaders import CalLookup_Description
from util import utils
from util.FileName import FileName
from util.ObsFile import ObsFile
import tables
import glob
from matplotlib.dates import strpdate2num
import fnmatch

class CalLookupFile:
    def __init__(self,path=None,mode='r'):
        if path is None:
            path = os.environ['MKID_CAL_LOOKUP']
        self.path = path
        
        if not (mode in ['r','a']):
            raise TypeError('Mode must be in [\'r\',\'a\']. Write (\'w\') is not allowd. Received: {}'.format(mode))
        self.file = tables.openFile(self.path,mode=mode)
        self.table = self.file.root.lookup

    def lookupObs(self,obsTstamp):
        result = self.table.readWhere('obs_tstamp == obsTstamp')
        if len(result) > 1:
            print 'multiple matches! Returning first'
        if len(result) == 0:
            print 'no matches'
            return None
        #convert to dict
        resultDict = {name:result[0][name] for name in result[0].dtype.names}
        return resultDict

    def updateObs(self,obsTstamp=None,newParams={}):
    
        if obsTstamp is None:
            try:
                obsTstamp = newParams['obs_tstamp']
            except KeyError:
                raise TypeError('Either obsTstamp or a dict containing \'obs_tstamp\' must be given')
        else:
            newParams['obs_tstamp'] = obsTstamp
        
        existingEntry = self.table.readWhere('obs_tstamp == obsTstamp')
        if len(existingEntry) == 0:
            if not ('obs_run' in newParams.keys() and 'obs_date' in newParams.keys()):
                raise TypeError('For a new entry, \'obs_run\' and \'obs_date\' musb be in newParams')
            print 'adding obs {} to table'.format(obsTstamp)
            newRow = self.table.row
            for key in newParams.keys():
                newRow[key] = newParams[key]
            newRow.append()
        elif len(existingEntry) == 1:
            print 'entry found.  Updating entry'
            for oldRow in self.table.where('obs_tstamp == obsTstamp'):
                for key in newParams.keys():
                    oldRow[key] = newParams[key]
                oldRow.update()
            self.file.flush()
        
    def getComponents(self,obsTstamp,keyPrefix):
        if keyPrefix == 'calSoln':
            keyPrefix = 'waveSoln'
        entry = self.lookupObs(obsTstamp)
        run = entry[keyPrefix+'_run']
        date = entry[keyPrefix+'_date']
        tstamp = entry[keyPrefix+'_tstamp']
        
        return {'run':run,'date':date,'tstamp':tstamp}

    def makeFileName(self,obsTstamp,keyPrefix):
        entry = self.lookupObs(obsTstamp)
        if entry is None:
            print 'entry none'
            return None
        run = entry[keyPrefix+'_run']
        date = entry[keyPrefix+'_date']
        tstamp = entry[keyPrefix+'_tstamp']
        if run == '' and date == '' and tstamp == '':
            return None
        fn = FileName(run=run,date=date,tstamp=tstamp)
        return fn

    def obs(self,obsTstamp):
        fn = self.makeFileName(obsTstamp,'obs')
        if fn is None:
            return ''
        else:
            return fn.obs()

    def calSoln(self,obsTstamp):
        fn = self.makeFileName(obsTstamp,'waveSoln')
        entry = self.lookupObs(obsTstamp)
        if fn is None:
            return ''
        else:
            isMasterCal = entry['waveSoln_isMasterCal']
            if isMasterCal:
                print 'master cal'
                return fn.mastercalSoln()
            else:
                print 'calSoln'
                return fn.calSoln()

    def flatSoln(self,obsTstamp):
        fn = self.makeFileName(obsTstamp,'flatSoln')
        entry = self.lookupObs(obsTstamp)
        if fn is None:
            return ''
        else:
            print fn
            isIllumCal = entry['flatSoln_isIllumCal']
            if isIllumCal:
                return fn.illumSoln()
            else:
                return fn.flatSoln()

    def fluxSoln(self,obsTstamp):
        fn = self.makeFileName(obsTstamp,'fluxSoln')
        if fn is None:
            return ''
        else:
            return fn.fluxSoln()

    def timeMask(self,obsTstamp):
        fn = self.makeFileName(obsTstamp,'timeMask')
        if fn is None:
            return ''
        else:
            return fn.timeMask()
        
    def cosmicMask(self,obsTstamp):
        fn = self.makeFileName(obsTstamp,'cosmicMask')
        if fn is None:
            return ''
        else:
            return fn.cosmicMask()

    def timeAdjustments(self,obsTstamp):
        fn = self.makeFileName(obsTstamp,'timeAdjustments')
        if fn is None:
            return ''
        else:
            return fn.timeAdjustments()

    def beammap(self,obsTstamp):
        fn = self.makeFileName(obsTstamp,'beammap')
        if fn is None:
            return ''
        else:
            return fn.beammap()

    def centroidList(self,obsTstamp):
        fn = self.makeFileName(obsTstamp,'centroidList')
        if fn is None:
            return ''
        else:
            return fn.centroidList()

    def populateTimeMasks(self,runPath):
        for root,dirs,files in os.walk(runPath):
            obsFilenames = fnmatch.filter(files,'obs*.h5')
            for obsFilename in obsFilenames:
                obsPath = os.path.join(root,obsFilename)
                
                try:
                    obs = ObsFile(obsPath)
                    fn = FileName(obsFile=obsPath)
                    params = {}
                    params['obs_run'],params['obs_date'],params['obs_tstamp'] = fn.getComponents()
                    params['timeMask_run'],params['timeMask_date'],params['timeMask_tstamp'] = fn.getComponents()
                    print params
                    self.updateObs(newParams=params)
                except:
                    pass

    def populateWaveCals(self,runPath):
        for root,dirs,files in os.walk(runPath):
            obsFilenames = fnmatch.filter(files,'obs*.h5')
            for obsFilename in obsFilenames:
                obsPath = os.path.join(root,obsFilename)
                
                try:
                    obs = ObsFile(obsPath)
                    fn = FileName(obsFile=obsPath)
                    params = {}
                    params['obs_run'],params['obs_date'],params['obs_tstamp'] = fn.getComponents()
                    params['waveSoln_run'] = params['obs_run']

                    obs.loadBestWvlCalFile()
                    waveCalPath = obs.wvlCalFileName
                    waveCalFilename = os.path.basename(waveCalPath)
                    if waveCalFilename.startswith('master'):
                        params['waveSoln_isMasterCal'] = True
                        tstamp = waveCalFilename.split('_')[1].split('.')[0]
                        params['waveSoln_tstamp'] = tstamp
                        params['waveSoln_date'] = ''
                    else:
                        params['waveSoln_isMasterCal'] = False
                        tstamp = waveCalFilename.split('_')[1].split('.')[0]
                        params['waveSoln_tstamp'] = tstamp
                        dirs = os.path.dirname(os.path.normpath(waveCalPath)).split(os.sep)
                        params['waveSoln_date'] = dirs[-1]
                    print params
                    self.updateObs(newParams=params)
                except:
                    pass
    

def writeNewCalLookupFile(path=None):
    if path is None:
        path = os.environ['MKID_CAL_LOOKUP']
    if os.path.exists(path):
        confirmResp = utils.confirm('File {} exists. Are you sure you want to DELETE it and make a new one?'.format(path),defaultResponse=False)
        if not confirmResp:
            return
            
    try:
        file = tables.openFile(path,mode='w')
    except:
        print 'Error: Couldn\'t create file, ',path
        return
    print 'writing to',path

    zlibFilter = tables.Filters(complevel=1, complib='zlib', fletcher32=False)
    table = file.createTable(file.root, 'lookup', CalLookup_Description,title='Cal Lookup Table',filters=zlibFilter)
    file.flush()
    file.close()

if __name__=='__main__':
    path = os.environ['MKID_CAL_LOOKUP']
    if not os.path.exists(path):
        writeNewCalLookupFile(path)
    lookup = CalLookupFile(path,mode='a')
    #lookup.populateWaveCals('/ScienceData/PAL2014')
    lookup.populateTimeMasks('/ScienceData/PAL2014')

    print lookup.obs('20141021-091336')
    
