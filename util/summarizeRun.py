import os,sys,glob,re
from util import FileName
from util import ObsFile
class SummarizeRun():
    def __init__(self,run, dateLimkidDataDir=None,dateList="all"):
        self.mkidDataDir = os.getenv('MKID_RAW_PATH', '/ScienceData')
        self.lines = []
        if dateList == "all":
            dirList = os.listdir(os.path.join(self.mkidDataDir,run))
        else:
            dirList = dateList
        dirList.sort()
        print "dirList=",dirList
        self.info = []
        for thisDir in dirList:
            fullDir = os.path.join(self.mkidDataDir,run,thisDir)
            print "Now search fillDir=",fullDir
            if len(thisDir) == 8:
                fileList = filter(os.path.isfile, glob.glob(fullDir + "/*.h5"))
            #fileList.sort(key=lambda x: os.path.getmtime(x))
                fileList.sort(key=lambda x: re.split(r'[/._]+', x)[-2])
                for fileName in fileList:
                    tokens = re.split(r'[/._]+', fileName)
                    date = tokens[-4]
                    flavor = tokens[-3]
                    ts = tokens[-2]
                    fileName = FileName.FileName(run,date,ts)
                    target = ""
                    description = ""
                    if flavor == "obs":
                        h5 = fileName.obs()
                        try:
                            of = ObsFile.ObsFile(fileName.obs())
                            target = of.getFromHeader('target')
                            description = of.getFromHeader('description')
                        except:
                            pass
                        del h5
                    line = '%s,%s,%s,%s,"%s","%s"'%(run,date,ts,flavor,target,description)
                    self.lines.append(line)
if __name__ == '__main__':
    try:
        run = sys.argv[1]
        #print "len is ",len(sys.argv)
        if len(sys.argv) > 2:
            dateList = sys.argv[2:]
            if dateList[0] == 'all':
                dateList = 'all'
        else:
            dateList = "all"
    except:
        print "usage:  summarizeRun.py run dateList"
        print " where "
        print "  run is the run (PAL2014, for example)"
        print "  dateList is the list of dates; default all get all of them"
        exit(0)
    print 'run=',run, 'dateList=',dateList
    sr = SummarizeRun(run,dateList)
    print "number of lines is ",len(sr.lines)
    csvFile = open(run+".csv","w")
    for line in sr.lines:
        print line
        csvFile.write(line+"\n")
    csvFile.close()
