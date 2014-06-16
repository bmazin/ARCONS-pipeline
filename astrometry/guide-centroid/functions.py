import pyfits
import numpy as np
from astropy import wcs
from scipy.optimize import leastsq,fsolve
import warnings
import os

def convert(coordinates,returnFormat=False):
    """
    Input a list of RA/DEC coordinate in degrees or hour format and converts it to the other format. To do, use ( ) library to eliminate these functions

    Keyword arguments:
    coordinate -- RA/DEC in degrees or hr (example,cooridnate=[[100.3,39],[54,21]] or coordinate=[[-1:2:3,22:0:52.2]])
    returnFormat -- if False, the function returns the list of converted coordinate, else, the function returns a dictionary with keywords   
                    coordinates and format.
    """    

    def hr2deg(string):
        nlist = string.split(':')
        hr = int(nlist[0])
        minute = int(nlist[1])
        sec = np.float(nlist[2])
        deg = (hr * (360.0/24.0)) + (minute * (360.0/(24.0*60))) + (sec * (360.0/(24.0*3600.0)))
        return deg

    def arc2deg(string):    
        nlist = string.split(':')
        negative = False
        if nlist[0][0] == '-':
            negative = True                 
        deg = int(nlist[0])
        arcm = int(nlist[1])
        arcs = np.float(nlist[2])
        deg = deg + (arcm * (1/60.0)) + (arcs * (1/3600.0))
        if negative:
            deg = -deg
        return deg

    def deg2hr(deg):
        deg = np.float(deg)
        sec = deg *(24*3600.0/360)
        minutes = int(sec/60)
        sec = sec - minutes * 60
        hr = int(minutes/60.0)
        minutes = minutes - hr *60
        ra = '%i:%i:%.2f' %(hr,minutes,sec)
        return ra
    
    def deg2arc(deg):
        if deg < 0:
            negative = True
        else:
            negative = False
        deg = abs(np.float(deg))
        arcs = deg * 3600.0
        arcm = int(arcs/60.0)
        arcs = arcs - arcm*60
        deg = int(arcm/60.0)
        arcm = arcm - deg*60
        if negative:
            dec = '-%i:%i:%.2f' %(deg,arcm,arcs)
        else:
            dec = '%i:%i:%.2f' %(deg,arcm,arcs)
        return dec

    def deg2std(sList):         
        stdList = []            
        for cr in range(len(sList)):
            ra = deg2hr(sList[cr][0])
            dec = deg2arc(sList[cr][1])
            stdList.append([ra,dec])
        return stdList  

    def std2deg(sList):
        degList = []
        for cr in range(len(sList)):
            deg1 = hr2deg(sList[cr][0])
            deg2 = arc2deg(sList[cr][1])
            degList.append([deg1,deg2])
        return degList
    
    if type(coordinates[0][0]) == str or type(coordinates[0][0]) == np.string_: 
        coordinates = std2deg(coordinates)
        format = 'standard'
    else:   
        format = 'standard'
        coordinates = deg2std(coordinates)
        format = 'degrees'
    if returnFormat:
        return {'coorinates':coordinates,'format':format}
    else:
        return coordinates

def poly2(p,x1,x2):
    """
    Polynomial distortion model following SIP convention
        
    Keywords argument:
    p --- ordered coeffient list
    x1 --- x coordinate
    x2 --- y coordinate
    """
    a10,a01,a20,a02,a11,a21,a12,a30,a03 = p
    y = a10*x1 + a01*x2 + a20*x1**2 + a02*x2**2 + a11*x1*x2 + a21*(x1**2)*x2 + a12*x1*x2**2 + a30*x1**3 + a03*x2**3
    return y
    
def pix2world(fitsImageName,pix,paramFile=None,crval1=None,crval2=None):
    """
    Input a LIST of pixel coordinates and output a list of world coordinates in degrees (example, pix=[[55,23],[34,43],[142,888]]). 
    If param file is provided, apply distortion to the coordinate conversion.
    """
    
    print pix
    #apply distortion params if it exist prior to conversion
    if paramFile != None:
        pix = distApp(pix,paramFile,crval1,crval2)
    
    imageList = pyfits.open(fitsImageName)
    header = imageList[0].header    

    #this fixes the error in header in which the length of CTYPE1 doesn't match the len of CTYPE2
    if header['CTYPE1'] == 'RA--TAN':
        header['CTYPE1'] = 'RA---TAN'
    
    #coordinates conversion
    w = wcs.WCS(imageList[0].header)
    world = w.wcs_pix2world(pix,1)
    world = world.tolist()
    return world

def world2pix(fitsImageName,world,paramFile=None,crval1=None,crval2=None):
    """
    Input a LIST of world coordinates in degrees and output a list of pixel coordinates (example, world=[[103.2,67],[28,33],[69,11]])
    If param file and crvals are provided, apply distortion to the coordinate conversion. Crvals are needed because SIP is based on the coordinates
    that are relative to the reference coordinate. If crvals provided here are differ from the ones in param file, error will be raised.
    """   
    #load file
    imageList = pyfits.open(fitsImageName)
    header = imageList[0].header
    #fix error
    if header['CTYPE1'] == 'RA--TAN':
        header['CTYPE1'] = 'RA---TAN'
    #coordinates conversion
    w = wcs.WCS(header)
    pix = w.wcs_world2pix(world,1)
    pix = pix.tolist()
    
    #Then if distortion params provided, find the inverse of that and apply it to find the final pixel coordinates      
    if paramFile != None:
        paramFile = np.load(paramFile)
        dt1 = paramFile['dt1']
        dt2 = paramFile['dt2']
        x1ref = paramFile['xref'][0]
        x2ref = paramFile['xref'][1]
        #Prevent python errors of too many files opened
        paramFile.close()
                    
        if crval2 == None or crval1 == None:
            raise ValueError('You must provide reference pixels to apply distortion!')    
        if int(crval1) != int(x1ref) or int(crval2) != int(x2ref):
            raise ValueError('Mismatch in reference coordinates between image and paramters file!')
        
        def _equations(p,dt1,dt2,x,y,x1ref,x2ref):
            """ 
            Return a system of equations roots will give the reverse transformation of distortion. Remember everything is in relative coordinates.
            """
            x1,x2 = p
            eqn1 = x - poly2(dt1,x1,x2) - x1 
            eqn2 = y - poly2(dt2,x1,x2) - x2 
            return (eqn1,eqn2)
        
        
        for index in range(len(pix)):
            x = pix[index][0] - x1ref
            y = pix[index][1] - x2ref
        
            x1,x2 = fsolve(_equations,(x,y),args=(dt1,dt2,x,y,x1ref,x2ref))
            
            pix[index][0] = x1 + x1ref
            pix[index][1] = x2 + x2ref
            
    return pix
    
def sortList(fileList,sec=True):
    '''
    Sort the file list in time ascending order. The first six letters in the file name has to be the time stamp (eg.061234)
    It assumes all files have same suffix in this case.
    '''
    returnList = [] 
    tempList1 = []
    
    for nfile in fileList:
        
        timeStamp = nfile[0:6]
        sec = timeConvert(timeStamp)
        tempList1.append(sec)
        #actually this only has to be done once since it assume all files have same suffix
        suffix = nfile[6:]
    
    #time-ascending order
    tempList1.sort() 
    
    #return in seconds format if sec is True
    if Sec:
        return tempList1
    else:
        for ntime in tempList1:
            timeStr = timeConvert(ntime)
            outputFileName = '%s%s' %(timeStr,suffix)
            returnList.append(outputFileName)
            return returnList
 
def distApp(pix,paramFile,crval1,crval2):
    """
    Apply distortion polynomial to the coordinates in objectList. (crval1,crval2) represents the reference coordinate in which the paraFile is calculated from.
    """
    paramFile = np.load(paramFile)
    dt1 = paramFile['dt1']
    dt2 = paramFile['dt2']
    x1ref = paramFile['xref'][0]
    x2ref = paramFile['xref'][1]
    #Prevent python errors of too many files opened
    paramFile.close()
    if crval2 == None or crval1 == None:
        raise ValueError('You must provide reference pixels to apply distortion!')    
    if int(crval1) != int(x1ref) or int(crval2) != int(x2ref):
        raise ValueError('Mismatch in reference coordinates between image and paramters file!')      
    #Apply distortion        
    for index in range(len(pix)):
        #Change to relative coordinates then apply distortion parameters
        x = pix[index][0] - x1ref
        y = pix[index][1] - x2ref
        pix[index][0] += poly2(dt1,x,y)
        pix[index][1] += poly2(dt2,x,y)
    return pix

def updateHeader(fitsImage,keyword,arg):
    warnings.filterwarnings("ignore")
    imageList = pyfits.open(fitsImage,mode='update')
    header = imageList[0].header
    header.update(keyword,arg)
    #notice for pyfits in turk, header[keyword]=arg will raise error if keyword is not already present. However, it won't be a problem in a newer version of pyfits. But header.update(..) always works.
    #header[keyword] = arg
    imageList.flush()

def distHeaderUpdate(fitsImage,saveName,paramFile):
    paramFile = np.load(paramFile)
    dt1 = paramFile['dt1']
    dt2 = paramFile['dt2']
    coeffList = paramFile['coeffList']
    coeffheaderUpdate(saveName,fitsImage,dt1,dt2,coeffList)
 
def coeffUpdate(saveName,fitsImageName,coeffx,coeffy,coeffList):           
    #Open the image
    imageList = pyfits.open(fitsImageName)
    header = imageList[0].header    
    #Updating header
    header['A_ORDER'] = 3
    header['B_ORDER'] = 3
    header['CTYPE1'] = 'RA---TAN-SIP'
    header['CTYPE2'] = 'DEC--TAN-SIP'
    count = 0
    for string in coeffList:
        p = int(string[0])
        q = int(string[1])
        keyword1 = 'A_%s_%s' %(p,q)
        keyword2 = 'B_%s_%s' %(p,q)
        header[keyword1] = coeffx[count]
        header[keyword2] = coeffy[count]

        count += 1
    #Try to remove any existing fits files  
    try:
        os.remove(saveName)
    except:
        pass           
    imageList.writeto(saveName,output_verify='ignore')    
    
def readCat(catFile,flagList=[0,1,4]):
    #read the cat file generated by sextractor. Only return stars that have flags 0, 1, or 4
    nfile = open(catFile)
    starList = []

    while True:
        line = nfile.readline()
        #end the reading of the file at the last line
        if line == '':
            break
        
        #ignore the comment sentence
        if str(line.split()[0])[0] == '#':
            pass
        else:
        #flag 0 1 and 4 usually indicate good source whereas other flags represent bad or possibily misidentified source
            try:
                if int(line.split()[3]) in flagList:
                    x = float(line.split()[1])
                    y = float(line.split()[2])
                    starList.append([x,y])
        #this argument is intended for manual catalog where flag and index are not requied to present. They are splitted by blacnk spaces
            except:
                x = float(line.split()[0])
                y = float(line.split()[1])
                starList.append([x,y])
    return starList
      

def timeConvert(timeStamp):
    '''
    eg convert '021456' expression to seconds in integer or convert seconds in integer into string expression
    '''
    if isinstance(timeStamp,str):
        hr = int(timeStamp[0:2])
        minu = int(timeStamp[2:4])
        sec = int(timeStamp[4:6])
        sec = hr*3600 + minu*60 + sec
        return sec
    elif isinstance(timeStamp,int):
        hr = timeStamp/3600
        minu = (timeStamp-(hr*3600))/60
        sec = timeStamp-(hr*3600)-(minu*60)       
        hr = str(hr)
        minu = str(minu)
        sec = str(sec)        
        if len(hr) == 1:
            hr = '0' + hr
        if len(minu) == 1:   
            minu = '0' + minu
        if len(sec) == 1:
            sec = '0' + sec
        return hr+minu+sec
    else:
        raise ValueError, 'Error in Timestamp Converting Type!'
