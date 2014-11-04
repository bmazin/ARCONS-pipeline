from subprocess import Popen, PIPE
import shlex
import httplib as hlib
import urllib as ulib
import xml.etree.ElementTree as ET	
import astropy.io.votable as VT
import warnings


def queryVizier(fitsTableName,format='asu-fits',source='USNO-B1.0',pos='104.9566125,14.2341555',radius=2):
	
	'''
	documentation: http://cdsarc.u-strasbg.fr/doc/asu-summary.htx
	'''
	
	print 'loading Vizier fits table'
	
	conn2 = hlib.HTTPConnection('vizier.u-strasbg.fr')
	params2 = ulib.urlencode({'-source':source,'-c':pos,'-c.rm':radius,'-out':'_RA _DE'})
	address = '/viz-bin/'+format
	
	conn2.request('POST',address,params2)
	resp2 = conn2.getresponse()

	fitsTable = open(fitsTableName,'w')
	fitsTable.write(resp2.read())

	conn2.close()
	fitsTable.close()
	print 'done'

def queryFitsImage(fitsImage,voTableName,pos='104.9566125,14.2341555'):

	'''
	retrieve fits image from 2MASS http://irsa.ipac.caltech.edu/applications/2MASS/IM/docs/siahelp.html
	note that fitsImage has to be in gz file in order to be extracted and opened correctly
	'''
	
	print 'loading 2MASS fits images'
	conn1 = hlib.HTTPConnection('irsa.ipac.caltech.edu')
	params1 = ulib.urlencode({'POS':pos,'INTERSECT':'ENCLOSED'})
	conn1.request('POST','/cgi-bin/2MASS/IM/nph-im_sia',params1)
	resp1 = conn1.getresponse()

	xml = open(voTableName,'w')
	xml.write(resp1.read())

	conn1.close()
	xml.close()

	#parsing the votable and download the selected image
	tree = ET.parse(voTableName)
	root = tree.getroot()

	#map all null values to 0 which fixes the error when open the votable with astropy
	for name in root.iter('TD'):
		if name.text == 'null':
			name.text = 0
	tree.write(voTableName)

	#warning suppression when open votable
	print 'downloading fits image...'
	warnings.resetwarnings()
	warnings.filterwarnings('ignore', category=Warning, append=True)
	voTable = VT.parse(voTableName)
	warnings.resetwarnings()
	warnings.filterwarnings('always', category=Warning, append=True)
	
	#raise warning if the image cannot be found
	try:
		table = voTable.get_first_table()
	except:
		raise ValueError('Cannot Locate the fits image!')
	 
	imageUrls = table.array['download']
	print imageUrls
	download = imageUrls[0]
	print download	
	ulib.urlretrieve(download,fitsImage)
	print 'done'


if __name__ == '__main__':
	
	
	pos = 'arp147'
	queryVizier('test.fits',pos=pos)
	queryFitsImage('test_image.fits','test_vo.xml',pos=pos)


	
	

