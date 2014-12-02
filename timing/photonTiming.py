'''
Author: Matt Strader                Date: Nov 19, 2014
For barycentering photon timestamps and folding them by a (pulsar) timing model to calculate phase
Uses the programs TEMPO or TEMPO2.
The main function is processTimestamps.  Other functions provide interfaces allowing barycentering on
text files of timestamps (barycenterTimestampFile) or photonLists with timestamps (timePhotonList).
'''

import matplotlib.pyplot as plt
import numpy as np
import os
import time
import subprocess
import numexpr
import shutil
import threading
import Queue
import multiprocessing
import warnings
from residuals import read_residuals
from util.popup import PopUp
from photonlist.photlist import PhotList


UNIX_EPOCH_JD = 2440587.5 #used to convert from unix timestamps to Julian Days
EPOCH_MJD = 2400000.5 #used to convert from Julian days to Modified Julian Days
SECS_PER_DAY = 86400 

def barycenterTimestampFile(timestampPath,outputPath,parFile,**kwargs):
    '''
    Barycenters timestamps in a text file
    INPUTS:
        timestampPath - a textfile containing timestamps to be barycentered. Must be MJD format (presumably TT(UTC)), one per line
        outputPath - a path to put the output text file containing the barycentered timestamps, this may be TDB or TCB depending on whether TEMPO or TEMPO2 are used
        parFile - a parameter file with values needed for timing, specific to the observed target.  See an example at ARCONS-pipeline/timing/test/B0531+21.par (Crab pulsar)
        kwargs - other arguments passed to processTimestamps.  See processTimestamps().
    '''
    mjdTimestamps = np.loadtxt(timestampPath)
    pDict = processTimestamps(mjdTimestamps=mjdTimestamps,parFile=parFile,**kwargs)
    np.savetxt(outputPath,pDict['baryMjdTimestamps'])

def timePhotonList(photonListPath,parFile,bPulsarTiming=False,nPhotonsPerProcess=1e6,verbose=False,nProcesses=(multiprocessing.cpu_count()-1),**kwargs):
    '''
    This function takes previously created photon lists (either standard photon lists or pulsar photon lists with empty columns for timing information) and feeds them through tempo/tempo2 to fill the mjd,baryMjd,(pulseNumber,phase) columns.  It then creates indices in the photon lists so that these columns may be searched efficiently. 

    INPUTS:
        photonListPath - path to an existing PhotonList HDF5 file, created with photlist.writePhotonList
        parFile - a parameter file with values needed for timing, specific to the observed target.  See an example at ARCONS-pipeline/timing/test/B0531+21.par (Crab pulsar)
        bPulsarTiming - boolean to determine whether pulsar specific columns (pulsePhase,pulseNumber) need to be filled.  If set to False, only the mjd,baryMjd columns are filled
        nPhotonsPerProcess - the number of photons that should be processed at a time.
        verbose - if True, debugging information is output to the terminal
        nProcesses - The number of processors the calculation should be split amongst.
        kwargs - other arguments passed to processTimestamps or saveResultsWorker.  See processTimestamps().
    '''
    #open the file
    photonList = PhotList(photonListPath,mode='r+')
    #assume the JD in the photon list header is correct
    jd = photonList.header.cols.jd[0]
    mjdStart = jd - EPOCH_MJD

    nPhotons = len(photonList.photTable)
    nChunks = int(np.ceil(1.*nPhotons/nPhotonsPerProcess))
    if verbose:
        print nPhotons,' photons to be timed'

    if nProcesses > 1: #Split the work over multiple processes

        #We will split the work into chunks of nPhotonsPerProcess.
        #For each chunk, we'll start a process which will call processTimestamps to barycenter it's chunk
        #we'll also start a thread which will take the barycentered timestamps and other values and stuff them
        #back into the photonList.

        #The first part is done by separate processes so each can have its own environment with working directory.
        #The second part is done by threads so they can easily share access to the photonList object

        inputQueue = multiprocessing.Queue() #we'll put unbarycentered timestamps, and metadat in this one
        resultsQueue = multiprocessing.Queue() #barycentered timestamps and metadata will be placed in this queue
        #by the various processes, to be stuffed into the photonList
        photonListLock = threading.Lock() #this will be used to ensure that multiple threads don't try to write to the photonList at the same time

        timingWorkers = [] #this will be a list of active processes working on timestamps
        resultsWorkers = [] #this will be a list of active threads stuffing values into the photonList
        for iChunk in xrange(nChunks):
            if verbose:
                print 'processing chunk #',iChunk,'of',nChunks

            #start a process that will invoke processTimestamps and start a tempo instance once we give it an input
            #via the inputQueue.  It will put the results in a resultsQueue to be handled by a resultsWorker thread.
            timingWorkerKwargs = kwargs.copy()
            timingWorkerKwargs['parFile'] = parFile
            timingWorkerKwargs['verbose'] = verbose
            timingWorker = multiprocessing.Process(target=__processTimestampsWorker,args=(inputQueue,resultsQueue),kwargs=kwargs)
            timingWorker.start()#start so it will be waiting for something in the inputQueue
            timingWorkers.append(timingWorker)#add to list of active processes

            #start a thread that will receive results from the above process via the resultsQueue, and put it in the photonList
            resultsWorkerKwargs = {'bPulsarTiming':bPulsarTiming,'verbose':verbose}
            resultsWorker = threading.Thread(target=__saveResultsWorker,args=(resultsQueue,photonListLock,photonList),kwargs=resultsWorkerKwargs)
            resultsWorker.start()#start so it will be waiting for results to process
            resultsWorkers.append(resultsWorker)#add to list of active threads

            #calculate indices that will be used to read chunks of the photonList, and to put processed chunks in
            #to the proper part of each column
            iPhotonStart = iChunk*nPhotonsPerProcess
            iPhotonEnd = (iChunk+1)*nPhotonsPerProcess
            data = photonList.photTable.read(iPhotonStart,iPhotonEnd)
            nPhotonsRead = len(data)

            if verbose:
                print nPhotonsRead,' photons to process in chunk'
            #read timestamps from photonList in secs since the beginning of the file.
            #convert to Modified Julian Days
            timestamps = np.array(data['arrivalTime'],dtype=np.float64)
            mjdTimestamps = mjdStart+timestamps/SECS_PER_DAY

            #send the timingWorker some timestamps to process.  The iPhotonStart,iPhotonEnd will be passed along with the results to the resultsWorker thread, so it knows where to save results
            inputDict = {'mjdTimestamps':mjdTimestamps,'iPhotonStart':iPhotonStart,'iPhotonEnd':iPhotonEnd}
            inputQueue.put(inputDict)

            #once we've filled our available processors
            #wait for the existing workers to finish
            if iChunk % nProcesses == nProcesses-1:
                if verbose:
                    print 'waiting for existing processes to finish'
                for worker in timingWorkers:
                    worker.join()
                for worker in resultsWorkers:
                    worker.join()
                timingWorkers = []
                resultsWorkers = []

        #wait for anything leftover to finish
        if verbose:
            print 'waiting for existing processes to finish'
        for worker in timingWorkers:
            worker.join()
        for worker in resultsWorkers:
            worker.join()
        timingWorkers = []
        resultsWorkers = []
        #by now the workers have processes all the timestamp chunks and saved them to the photonList file
        
    else:
        for iChunk in xrange(nChunks):
            if verbose:
                print 'processing chunk #',iChunk,'of',nChunks
            #read the chunk from the photonList
            iPhotonStart = iChunk*nPhotonsPerProcess
            iPhotonEnd = (iChunk+1)*nPhotonsPerProcess
            data = photonList.photTable.read(iPhotonStart,iPhotonEnd)
            nPhotonsRead = len(data)

            #convert the timestamps to Modified Julian Date
            if verbose:
                print nPhotonsRead,' photons to process in chunk'
            timestamps = np.array(data['arrivalTime'],dtype=np.float64)#in secs since beginning of file
            mjdTimestamps = mjdStart+timestamps/SECS_PER_DAY

            #Process the timestamps with tempo/tempo2
            if verbose:
                print 'running timing program'
            timingDict = processTimestamps(mjdTimestamps,parFile,verbose=verbose,**kwargs)
            if verbose:
                print 'saving results for chunk'
            
            #stuff the results into the appropriate column
            photonList.photTable.cols.mjd[iPhotonStart:iPhotonEnd] = mjdTimestamps
            photonList.photTable.cols.baryMjd[iPhotonStart:iPhotonEnd] = timingDict['baryMjdTimestamps']
            if bPulsarTiming:
                photonList.photTable.cols.pulseNumber[iPhotonStart:iPhotonEnd] = timingDict['pulseNumbers']
                photonList.photTable.cols.pulsePhase[iPhotonStart:iPhotonEnd] = timingDict['pulsePhases']
                
    photonList.file.flush()
    #index the newly filled columns
    try:
        if bPulsarTiming:
            colsToIndex = ['mjd','baryMjd','pulseNumber','pulsePhase']
        else:
            colsToIndex = ['mjd','baryMjd']
            
        for colName in colsToIndex:
            if photonList.colinstances[colName].is_indexed:
                warnings.warn('Column is already indexed: '+colName+' - skipping.')
            else:
                if verbose:
                    print 'Indexing column: '+colName
                photonList.colinstances[colName].createCSIndex()
    finally:
        photonList.file.flush()
    del photonList


def __processTimestampsWorker(inputQueue,resultsQueue,**kwargs):
    '''
    This function will be invoked as a separate process.  It takes a dict with mjdTimestamp values from inputQueue,
    sends them through processTimestamps, then packages the results into a dict to be stuffed in resultsQueue

    INPUTS:
        inputQueue - multiprocessing.Queue containing a dictionary with keys 'mjdTimestamps','iPhotonStart','iPhotonEnd'
        resultsQueue - multiprocessing.Queue where a dictionary containing results and inputs is put
        kwargs - other arguments passed to processTimestamps.  See processTimestamps().
    '''
    inputDict = inputQueue.get()
    if kwargs.get('verbose',False):
        print 'new data in queue, start timing program'
    resultsDict = processTimestamps(mjdTimestamps=inputDict['mjdTimestamps'],**kwargs)
    resultsDict.update(inputDict)#stuff mjdTimestamps,iPhotonStart,iPhotonEnd into results 
    #to be passed along
    resultsQueue.put(resultsDict)
    if kwargs.get('verbose',False):
        print 'done processing data'

def __saveResultsWorker(resultsQueue,photonListLock,photonList,bPulsarTiming=False,verbose=False):
    '''
    This function will be invoked as a thread.  It takes a dict from resultsQueue,
    and stuffs them into the appropriate columns of a photonList

    INPUTS:
        resultsQueue - multiprocessing.Queue containing a dictionary with results to be saved
        photonListLock - a lock shared with other threads, so that threads do not write to the photonList at the same time
        photonList - the photonList to be saved to
        bPulsarTiming - indicates whether the pulsar specific columns need to be saved
        verbose - if True, will print debugging information
    '''
    resultsDict = resultsQueue.get()
    iPhotonStart = resultsDict['iPhotonStart']
    iPhotonEnd = resultsDict['iPhotonEnd']
    photonListLock.acquire() #wait for other threads to finish writing to photonList, and then lock it for our own writing
    if verbose:
        print 'saving photons',iPhotonStart,'through',iPhotonEnd
    photonList.photTable.cols.mjd[iPhotonStart:iPhotonEnd] = resultsDict['mjdTimestamps']
    photonList.photTable.cols.baryMjd[iPhotonStart:iPhotonEnd] = resultsDict['baryMjdTimestamps']
    if bPulsarTiming:
        photonList.photTable.cols.pulseNumber[iPhotonStart:iPhotonEnd] = resultsDict['pulseNumbers']
        photonList.photTable.cols.pulsePhase[iPhotonStart:iPhotonEnd] = resultsDict['pulsePhases']
    photonListLock.release() #let other threads have access to photonList again
    if verbose:
        print 'done saving photons',iPhotonStart,'through',iPhotonEnd
        
        
def processTimestamps(mjdTimestamps,parFile,workingDir=os.getcwd(),timingProgram='tempo2',bCleanTempDir=True,verbose=False,timingProgramLog='/dev/null'):
    '''
    Takes in timestamps in unbarycentered MJD format, and runs them through a timing program (tempo/tempo2) to barycenter them and if analyzing a pulsar, to compute phase and pulseNumber for every (photon) timestamp

    INPUTS:
        mjdTimestamps - timestamps to be barycentered and/or folded.  Presumably in TT(UTC) format
        parFile - a parameter file with values needed for timing, specific to the observed target.  See an example at ARCONS-pipeline/timing/test/B0531+21.par (Crab pulsar)
        workingDir - a directory where a temporary directory can be created.  You must have permission to write to this directory
        timingProgram - 'tempo2' to run TEMPO2 and 'tempo' to run TEMPO.  This determines whether the output is in TCB (tempo2) or TDB (tempo) standard format as well as what pulsar binary models are allowed
        bCleanTempDir - when finished, the temporary directory created for inputs and outputs is deleted
        verbose - if True, will print debugging information

    OUPUTS:
        Return value is a dictionary with tags:
            'baryMjdTimestamps' - the barycentered timestamps in Modified Julian Days.
                TDB if tempo used, TCB if tempo2 used
            'pulsePhases' - the pulsar spin phase at a timestamp as a fraction of a period
            'pulseNumbers' - which pulses each timestamp corresponds to relative to the TZR reference time from the
                par file
    '''
    print 'workingDir',workingDir,os.path.exists(parFile)
    print 'parFile',parFile,os.path.exists(parFile)
    print 'abs',os.path.abspath(workingDir)
    nTimestamps = len(mjdTimestamps)
    workingTime = '{:.4f}'.format(time.time())
    tempDir = os.path.join(workingDir,'dir'+workingTime)
    print 'tempDir',tempDir
    parPath = os.path.abspath(parFile)
    parFilename = os.path.basename(parPath)
    print 'parPath',parPath,parFilename

    
    #make a temporary directory and switch to it, so the various output files will all be here
    priorCwd = os.getcwd()
    os.mkdir(tempDir)
    os.chdir(tempDir)
    #make a copy of the par file in the new temp directory, so tempo won't change the original file

    print 'newDir',os.getcwd(),os.path.exists(os.getcwd())
    shutil.copy(parPath,parFilename)
    print 'parPath',parPath,os.path.exists(parPath)
    if timingProgram == 'tempo2':
        tempJDPath = 'temp_toa.txt' #input file
        tempBJDPath = 'temp_bary_toa.txt' #output file

        #save the timestamps to the tempo2 input file
        np.savetxt(tempJDPath,mjdTimestamps)
    
        #form the tempo2 command using the arcons plugin
        strCommand = 'tempo2 -gr arcons -f {} -a {} -o {} > {}'.format(parFilename,tempJDPath,tempBJDPath,timingProgramLog)

        #run tempo2
        if verbose:
            print 'running timing program',strCommand
        proc = subprocess.Popen(strCommand,shell=True)
        proc.wait()

        if verbose:
            print 'parsing results'
        
        #parse the results
        tempoOutput = np.loadtxt(tempBJDPath,usecols=[1,2])
        nBaryTimestamps = len(tempoOutput)
        if nBaryTimestamps != nTimestamps:
            raise Exception('number of photons sent to tempo2 does not match number of photons returned')
        baryMjdTimestamps = tempoOutput[:,0]
        phases = tempoOutput[:,1]
        pulseNumbers = np.array(phases,dtype=np.uint32)
        pulsePhases = numexpr.evaluate('phases%1.')
        
    elif timingProgram == 'tempo':
        tempJDPath = 'temp_toa.tim' #input file
        tempPulseNumbersPath = 'pulseNumbers.txt' #optional output file

        #.tim format
        ##first line indicates using tempo2 style TOA's with keyword FORMAT
        #following lines follow the format:
        #photonLabel observing_frequency toa_mjd(utc) error observatory_site_code
        #the first toa line will be the reference time TZR, which defines phase 0

        #Example:
        #FORMAT 1 
        #TZR 0.000000000000 56273.164314974114 0.0 GB
        #arcons 0.0 56273.163858501241 0.0 PL

        formatLine = b'FORMAT 1\n'
        #read reference point TZR from par file
        parDict = readParFile(parFilename)
        tzrMjd = parDict['TZRMJD'][0]
        tzrFreq  = parDict['TZRFRQ'][0]
        tzrSite = parDict['TZRSITE'][0]
        tzrError = 0.
        tzrLine = b'TZR {} {} {:.1f} {}\n'.format(tzrFreq,tzrMjd,tzrError,tzrSite)

        #create the tempo .tim input file
        with open(tempJDPath,'wb') as tempJDFile:
            tempJDFile.write(formatLine)
            tempJDFile.write(tzrLine)
            toaLine = 'arcons 0.0 %.12f 0.0 PL' #PL for Palomar, frequency=0, label=arcons, error=0
            np.savetxt(tempJDFile,mjdTimestamps,fmt=toaLine)

        paramsToFit = 1
        #form the tempo command
        nToasToAllocate = nTimestamps+1
        strCommand = 'tempo -f {} -no {} -m{} -l{} {} > {}'.format(parFilename,tempPulseNumbersPath,nToasToAllocate,paramsToFit,tempJDPath,timingProgramLog)

        #Run tempo
        if verbose:
            print 'running timing program',strCommand
        proc = subprocess.Popen(strCommand,shell=True)
        proc.wait()

        if verbose:
            print 'parsing results'
    
        #parse barycentered timestamps, pulseNumbers, phases from output files
        residuals = read_residuals()
        #the reference TZR is the first photon, remove it
        baryMjdTimestamps = residuals.bary_TOA[1:] 
        pulsePhases = residuals.prefit_phs[1:]
        pulseNumbers = np.loadtxt(tempPulseNumbersPath)[1:]
        #some phases may be reported as negative, positive value would be 1+(-phase)
        pulsePhases[pulsePhases<0.] = pulsePhases[pulsePhases<0.]+1.
        if len(baryMjdTimestamps) != nTimestamps or len(pulseNumbers) != nTimestamps or len(pulsePhases) != nTimestamps:
            print 'number of pulseNumbers doesn\'t match input',nTimestamps,len(pulseNumbers)
            print tempPulseNumbersPath,os.getcwd()
            print pulseNumbers
            raise Exception('number of photons sent to tempo ({}) does not match number of photons returned ({})'.format(nTimestamps,len(pulseNumbers)))
        print 'results parsed'

        
    if bCleanTempDir:
        print 'removing tempDir'
        shutil.rmtree(tempDir)
        print 'tempDir',tempDir,os.path.exists(tempDir)
    os.chdir(priorCwd)
    return {'baryMjdTimestamps':baryMjdTimestamps,'pulsePhases':pulsePhases,'pulseNumbers':pulseNumbers}

def readParFile(parPath):
    '''
    parses a par file.  See the documentation for Tempo/Tempo2 for the format and avaliable keywords

    INPUTS:
        parFile - a parameter file with values needed for timing, specific to the observed target.  See an example at ARCONS-pipeline/timing/test/B0531+21.par (Crab pulsar)

    OUPUTS:
        a dictionary with a list of values for every key/line in the par file
    '''
    with open(parPath,'r') as parFile:
        parDict = {}
        for line in parFile:
            lineArr = line.split()
            if len(lineArr) >= 2:
                parDict[lineArr[0]] = lineArr[1:]
    return parDict

    

if __name__=='__main__':
    parFile = 'test/B0531+21.par'
    mjd0 = 56273.163858501241
    mjdTimestamps = np.arange(20)+mjd0
    pDict = processTimestamps(mjdTimestamps,parFile,timingProgram='tempo2')
    for key in pDict.keys():
        print key,pDict[key]
    

#def processTimestamps(mjdTimestamps,parFile,workingDir=os.getcwd(),timingProgram='tempo2',bCleanTempDir=True,referenceTime=None,verbose=False):
        

