#!/bin/python
#This program prints header information and description for every obs*.h5 file in a given folder.

import tables
import sys
import os
import glob

if len(sys.argv) != 2:
    print 'Usage: python ',sys.argv[0],' folderPath'
    exit(1)
obsPath = sys.argv[1]
for obs in sorted(glob.glob(os.path.join(obsPath,'obs*.h5'))+glob.glob(os.path.join(obsPath,'cal*.h5'))):

    try:
        f=tables.openFile(obs,'r')
        try:
            hdrNode = f.getNode('/header/header')
            hdr=hdrNode.read()
            print 'Filename: ',os.path.basename(obs)
            print 'Target: ',hdr['target'][0]
            print 'Local time: ',hdr['localtime'][0]
            print 'UTC: ',hdr['utc'][0]
            print 'Description: ',hdr['description'][0]
            print 'Filter: ',hdr['filt'][0]
            print 'Exposure Time: ',hdr['exptime'][0]
            print ''
            f.close()
        except ValueError:
            print 'Missing Header Entries'
        except tables.exceptions.HDF5ExtError:
            print os.path.basename(obs)
            print 'Can\'t Read Header Data'
        except IndexError:
            print 'Missing Header'
    except tables.exceptions.HDF5ExtError:
        print 'Can\'t open file'



