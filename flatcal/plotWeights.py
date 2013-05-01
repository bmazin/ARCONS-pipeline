import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from util.FileName import FileName

#        npzFileName = os.path.splitext(fullFlatCalFileName)[0]+'.npz'
#
#        #calculate total spectra and medians for programs that expect old format flat cal
#        spectra = np.array(np.sum(self.spectralCubes,axis=0))
#
#        wvlMedians = np.zeros(self.nWvlBins)
#        spectra2d = np.reshape(spectra,[self.nRow*self.nCol,self.nWvlBins ])
#        for iWvl in xrange(self.nWvlBins):
#            spectrum = spectra2d[:,iWvl]
#            goodSpectrum = spectrum[spectrum != 0]#dead pixels need to be taken out before calculating medians
#            wvlMedians[iWvl] = np.median(goodSpectrum)
#        np.savez(npzFileName,median=wvlMedians,medianSpectra=np.array(self.medianSpectra),binEdges=self.wvlBinEdges,spectra=spectra,weights=np.array(self.flatWeights.data))

def main():
    run='PAL2012'
    #flatSunsetDates = ['20121211','20121207','20121210','20121211']
    #flatTimestamps = ['20121212-084725','','','']
    #flatLabels = ['geminga 1211','twi 1207','twi 1210','twi 1211']
    flatSunsetDates = ['20121207','20121210','20121211']
    flatTimestamps = ['','','']
    flatLabels = ['twi 1207','twi 1210','twi 1211']
    oldFlatSunsetDates = ['20121207','20121210','20121211']
    oldFlatTimestamps = ['','','']
    oldFlatLabels = ['old twi 1207','old twi 1210','old twi 1211']

    flatFileNames = [FileName(run=run,date=date,tstamp=tstamp).flatInfo() for date,tstamp in zip(flatSunsetDates,flatTimestamps)]
    oldFlatFileNames = [FileName(run=run,date=date,tstamp=tstamp).oldFlatInfo() for date,tstamp in zip(oldFlatSunsetDates,oldFlatTimestamps)]
    flatInfos = [np.load(filename) for filename in flatFileNames]
    oldFlatInfos = [np.load(filename) for filename in oldFlatFileNames]

    pdfFileName = '/Scratch/flatCalSolnFiles2/compareFlatWeights.pdf'
    pp = PdfPages(pdfFileName)
    nPlotsPerRow = 2 
    nPlotsPerCol = 4 
    nPlotsPerPage = nPlotsPerRow*nPlotsPerCol
    iPlot = 0 

    matplotlib.rcParams['font.size'] = 4 

    nRow=46
    nCol=44
    for iRow in xrange(nRow):
        print iRow
        for iCol in xrange(nCol):

            if iPlot % nPlotsPerPage == 0:
                fig = plt.figure(figsize=(10,10),dpi=100)

            ax = fig.add_subplot(nPlotsPerCol,nPlotsPerRow,iPlot%nPlotsPerPage+1)
            ax.set_ylim(.5,2.)
            ax.set_xlabel(r'$\lambda$ ($\AA$)')
            ax.set_ylabel(r'Weights')
            ax.set_title('Flat Weights')

            ax.set_title('p %d,%d'%(iRow,iCol))

            
            allWeights = np.array([flatInfo['weights'][iRow,iCol] for flatInfo in flatInfos])
            
            if np.any(allWeights):
                for iFlatInfo,flatInfo in enumerate(oldFlatInfos):
                    #weights = flatInfo['weights'][iRow,iCol]
                    weights = flatInfo['weights'][iRow,iCol]
                    flatSpectra = flatInfo['spectra'][iRow,iCol]
                    flatMedians = flatInfo['median']
                    deltaFlatSpectra = np.sqrt(flatSpectra)
                    deltaWeights = weights*deltaFlatSpectra/flatSpectra
                    color=matplotlib.cm.jet((iFlatInfo+1.)/(len(flatInfos)*2))
                    #wvlBinCenters = self.wvlBinEdges[:-1]+np.diff(self.wvlBinEdges)/2.
                    #ax.step(self.wvlBinEdges[:-1],weights,linestyle='-',label=self.params['flatInfoFiles'][iFlat],color=color,)

                    ax.errorbar(flatInfo['binEdges'][0:-1],weights,linestyle='-',yerr=deltaWeights,color=color,label=oldFlatLabels[iFlatInfo],alpha=.7)
                for iFlatInfo,flatInfo in enumerate(flatInfos):
                    #weights = flatInfo['weights'][iRow,iCol]
                    weights = allWeights[iFlatInfo]
                    flatSpectra = flatInfo['spectra'][iRow,iCol]
                    flatMedians = flatInfo['median']
                    deltaFlatSpectra = np.sqrt(flatSpectra)
                    deltaWeights = weights*deltaFlatSpectra/flatSpectra
                    color=matplotlib.cm.jet((iFlatInfo+1.+len(flatInfos))/(len(flatInfos)*2))
                    #wvlBinCenters = self.wvlBinEdges[:-1]+np.diff(self.wvlBinEdges)/2.
                    #ax.step(self.wvlBinEdges[:-1],weights,linestyle='-',label=self.params['flatInfoFiles'][iFlat],color=color,)

                    ax.errorbar(flatInfo['binEdges'][0:-1],weights,linestyle='-',yerr=deltaWeights,color=color,label=flatLabels[iFlatInfo])

                ax.legend(loc='lower right')
                if iPlot%nPlotsPerPage == nPlotsPerPage-1 or (iRow == nRow-1 and iCol == nCol-1):
                    pp.savefig(fig)
                iPlot += 1

    pp.close()


if __name__=='__main__':
    main()
