
#Photometry Super class
import os
from util.readDict import readDict


class Photometry(object):

    def __init__(self,path='/Scratch/DisplayStack/RUN_TEMPLATE/TARGET_TEMPLATE'):
        '''
        Constructs a list of obs FileName objects from the dictionary in the path.
        'run' is assumed to be the second to last directory in the path
        
        Inputs:
            path - path to the display stack target info
        '''
        self.path = path
        run = os.path.basename(os.path.dirname(path))
        
        for f in os.listdir(path):
            if f.endswith(".dict"):
                #if self.verbose: print 'Loading params from ',path+os.sep+f
                self.params = readDict(path+os.sep+f)
                self.params.readFromFile(path+os.sep+f)
                
        self.obsFNs = []
        for i in range(len(self.params['sunsetDates'])):
            for j in range(len(self.params['obsTimes'][i])):
                obsFN = FileName(run=run, date=self.params['sunsetDates'][i], tstamp=self.params['utcDates'][i]+'-'+self.params['obsTimes'][i][j])
                self.obsFNs.append(obsFN)
        


















