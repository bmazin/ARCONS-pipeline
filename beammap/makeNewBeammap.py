'''
Author: Matt Strader                    Date: May 21, 2014
Tool to remap pixels from an existing beammap to create a new existing beammap
'''
from util.ObsFile import ObsFile
from util.FileName import FileName
import tables
from beammap import remapPixels

def main():
    #load the current PAL2012 beammap from an obs file from that run,
    #a sky file for hr9087
    run = 'PAL2012'
    date = '20121210'
    obsTimestamp = '20121211-051650'
    obsFN = FileName(run=run,date=date,tstamp=obsTimestamp)
    obsFileName = obsFN.obs()
    obs = ObsFile(obsFileName)

    #Load an existing definition of which pixels to remap where
    pixMap = remapPixels.PixelMap(obsFN.pixRemap())
    pixMapSourceList,pixMapDestList = pixMap.getRemappedPix()

    #create a new h5 file and copy the existing beammap node hierarchy into it
    newBeammapPath = 'beamimage_PAL2012_corrected.h5'
    newBeammapFile = tables.openFile(newBeammapPath, mode='w')
    newBeammapFile.copyNode(obs.file.root.beammap, newparent=newBeammapFile.root, newname='beammap', recursive=True)
    newBeammap = newBeammapFile.root.beammap.beamimage

    #for each array in the h5, move the pixels according to pixMap
    newBeammap[:,:] = pixMap.remapArray(obs.beamImage)
    newBeammapFile.root.beammap.resfreq[:,:] = pixMap.remapArray(newBeammapFile.root.beammap.resfreq.read())
    newBeammapFile.root.beammap.atten[:,:] = pixMap.remapArray(newBeammapFile.root.beammap.atten.read())

    newBeammapFile.flush()
    newBeammapFile.close()
    obs.file.close()

if __name__ == '__main__':
    main()
