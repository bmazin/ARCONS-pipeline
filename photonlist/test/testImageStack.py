'''
Author: Julian van Eyken            Date: May 31 2013
A bit of photon-list image stacking testing....
'''


import photonlist.photlist as pl
import photonlist.RADecImage as rdi
import os.path
import glob

def makeImageStack(fileNames='photon_*.h5',dir='/Users/vaneyken/Data/UCSB/ARCONS/turkDataCopy/Intermediate/photonLists/20121211'):
    '''
    Create an image stack
    INPUTS:
        filenames - string, list of photon-list .h5 files. Can either
                    use wildcards (e.g. 'mydirectory/*.h5') or if string
                    starts with an @, supply a text file which contains
                    a list of file names to stack. (e.g.,
                    'mydirectory/@myfilelist.txt', where myfilelist.txt 
                    is a simple text file with one file name per line.)
    '''
    
    #Get the list of filenames
    if fileNames[0]=='@':
        #(Note, actually untested, but should be more or less right...)'
        files=[]
        with open(fileNames[1:]) as f:
            for line in f:
                files.append(os.path.join(dir,line))
    else:
        files = glob.glob(os.path.join(dir, '*.h5'))

    virtualImage = rdi.RADecImage()
    for eachFile in files:
        if os.path.exists(eachFile):
            print 'Loading: ',os.path.basename(eachFile)
            #fullFileName=os.path.join(dir,eachFile)
            phList = pl.PhotList(eachFile)
            virtualImage.loadImage(phList,stack=True)
            virtualImage.display()
        else:
            print 'File doesn''t exist: ',eachFile
    
    return virtualImage

