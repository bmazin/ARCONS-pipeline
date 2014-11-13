
#Photometry Super class
import os
from util.readDict import readDict


class Photometry(object):

    def __init__(self,path='/Scratch/DisplayStack/RUN_TEMPLATE/TARGET_TEMPLATE',scienceDataPath = '/ScienceData',run = 'PAL2014'):
        '''
        Inputs:
            path - path to the display stack target info
        '''
        self.path = path
        
        for f in os.listdir(path):
            if f.endswith(".dict"):
                #if self.verbose: print 'Loading params from ',path+os.sep+f
                self.params = readDict(path+os.sep+f)
                self.params.readFromFile(path+os.sep+f)
                
        self.obs_name_list = []
        for i in range(len(self.params['sunsetDates'])):
            for j in range(len(self.params['obsTimes'][i])):
                self.obs_name_list.append(scienceDataPath+os.sep+run+os.sep+self.params['sunsetDates'][i]+os.sep+'obs_'+self.params['utcDates'][i]+'-'+self.params['obsTimes'][i][j]+'.h5')
        


















