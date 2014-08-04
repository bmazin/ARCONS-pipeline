import sys, os
import tables
from tables import *
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import pickle
#from RO import StringUtil
import time
from matplotlib import dates
import datetime
import time
import matplotlib.dates as mdates

class Location(IsDescription):
    fileName  = StringCol(30)
    dateObs   = StringCol(30)
    ra        = StringCol(12)
    dec       = StringCol(10)
    lst       = StringCol(8)
    cassang   = FloatCol()

def listFiles():

    h5file = openFile("guideCamera.h5", mode = "w", title = "Guide Camera Locations")
    group = h5file.createGroup("/", "group", "group")
    table = h5file.createTable(group, "location", Location, "Locations")
    row = table.row

    dataDir = os.getenv('MKID_RAW_PATH','.')
    print "listFiles in dataDir=",dataDir
    files = os.listdir(dataDir)
    print "number of files found:",len(files)
    files.sort()
    print "done sorting"
    iFile = 0
    for file in files:
        fullFileName=os.path.join(dataDir,file)
        hdulist = pyfits.open(fullFileName)
        ra  = hdulist[0].header.get('RA',default="unknown")
        dec = hdulist[0].header.get('DEC',default="unknown")
        lst = hdulist[0].header.get('LST',default="unknown")
        cassang = hdulist[0].header.get('CASSANG',default="unknown")
        dateObs = hdulist[0].header.get("DATE-OBS",default="unknown")
        if (iFile%100 == 0) :
            print "iFile=",iFile,"/",len(files)," file=",file,"  dateObs=",dateObs," ra=",ra, " dec=",dec," lst=",lst, "cassang=",cassang
        iFile += 1
        row["fileName"] = file
        row["dateObs"] = dateObs
        row["ra"] = ra
        row["dec"] = dec
        row["lst"] = lst
        row["cassang"] = cassang
        row.append()
    hdulist.close()

def plotRa():
    fmt = "%Y-%m-%dT%H:%M:%S"
    ral = []
    decl = []
    datel = []
    fid = tables.openFile("guideCamera.h5")
    for location in fid.root.group.location:
        ras = location['ra']
        decs = location['dec']
        dateObss = location['dateObs']
        ra = 15*myDegFromDMSStr(ras)
        dec = myDegFromDMSStr(decs)
        date = time.strptime(dateObss,fmt)
        ral.append(ra)
        decl.append(dec)
        # time.mktime(date) is seconds since the Epoch
        # dates.date2num(datetime.datetime.fromtimestamp(time.mktime(d))) is the matplotlib float format
        dateMpl = dates.date2num(datetime.datetime.fromtimestamp(time.mktime(date)))
        datel.append(dateMpl)

    # now make the actual plot
    plt.close('all')
    fig, ax = plt.subplots(1)
    ax.plot(datel, ral)

    # set the format of the dates on the x axis
    fig.autofmt_xdate()
    hfmt = dates.DateFormatter('%m/%d %H:%M')
    ax.xaxis.set_major_locator(dates.DayLocator())
    ax.xaxis.set_major_formatter(hfmt)
    # set the format of dates in the toolbar
    ax.fmt_xdata = mdates.DateFormatter('%y-%m-%dT%H:%M:%S')
    plt.show()
    return ral,decl,datel

def myDegFromDMSStr(DMSStr):
    a = DMSStr.split(":")
    d = float(a[0])
    m = float(a[1])
    s = float(a[2])
    retval = abs(d)+m/60.0+s/3600.0
    if d<0 : retval *= -1
    return retval

def myDegToHMS(deg):
    hours = deg/15.0
    absHours = abs(hours)



    
