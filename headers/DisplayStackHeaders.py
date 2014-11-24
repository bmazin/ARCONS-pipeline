
#Headers for DisplayStack
import tables
import numpy as np

# Create header and data group and table names
headerGroupName = 'header'
headerTableName = 'header'
imagesGroupName = 'stack'
imagesTableName = 'stack'
pixIntTimeTableName = 'pixel_integration_time'
timeTableName = 'time'
intTimeTableName = 'integration_time'

# Create lookup names for header information
runColName = 'run'
targetColName = 'targetName'
obsFileColName = 'obsFileName'
wvlCalFileColName = 'wvlCalFileName'
flatCalFileColName = 'flatCalFileName'
nRowColName = 'nRow'
nColColName = 'nCol'
RAColName = 'RA'
DecColName = 'Dec'
deadPixColName = 'deadPixFileName'
hotPixColName = 'hotPixFileName'
lowWvlColName = 'lowWvlCutoff'
highWvlColName = 'highWvlCutoff'
expTimeColName = 'exptime'
lstColName = 'lst'
integrationTimeColName = 'integrationTime'
HA_offsetColName = 'HA_offset'


class ImageStackHeaderDescription(tables.IsDescription):
    targetName = tables.StringCol(100, dflt='')
    run = tables.StringCol(100, dflt='')
    obsFileName = tables.StringCol(100, dflt='')
    wvlCalFileName = tables.StringCol(100, dflt=np.nan)
    flatCalFileName = tables.StringCol(100, dflt='')
    deadPixFileName = tables.StringCol(100, dflt='')
    hotPixFileName = tables.StringCol(100, dflt='')
    nCol = tables.UInt32Col(dflt=-1)
    nRow = tables.UInt32Col(dflt=-1)
    lowWvlCutoff = tables.Float64Col(dflt=np.nan)
    highWvlCutoff = tables.Float64Col(dflt=np.nan)
    exptime = tables.Float64Col(dflt=np.nan)
    lst = tables.StringCol(100, dflt='')
    integrationTime = tables.Float64Col(dflt=np.nan)
    RA = tables.StringCol(100, dflt='')
    Dec = tables.StringCol(100, dflt='')
    HA_offset = tables.Float64Col(dflt=0.0)

    
