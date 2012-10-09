#!/usr/bin/env python
# resolve an object name to a position using the Sesame service.
# Author: Noel Winstanley noel.winstanley@manchester.ac.uk 2006 - 2008
# run with --help for usage information.
import xmlrpclib
import sys
import os
import optparse
#verify we're running a suitable version of python
if not (sys.version_info[0] > 2 or sys.version_info[1] >= 5):
    print """This script runs best on python 2.5 or above. 
    You're running """ + sys.version
    
#connect to acr
fname = os.path.expanduser("~/.astrogrid-desktop")
assert os.path.exists(fname),  'No AR running: Please start your AR and rerun this script'
prefix = file(fname).next().rstrip()
acr = xmlrpclib.Server(prefix + "xmlrpc")

parser = optparse.OptionParser(usage="%prog [options] <object-name>",
                       description="resolve an object name to position using the CDS Sesame service")
parser.add_option('-a','--all',action='store_true',default=False,
                  help='Display full information for this object')
parser.add_option('-d','--decimal',action='store_true',default=True
                  , help="Return position in decimal degrees (default)")
parser.add_option('-s','--sexa',action='store_false',dest='decimal'
                  ,help="Return position in sexagesimal")
parser.add_option('-e','--examples', action='store_true', default=False
                  , help='display some examples of use and exit.')

(opts,args) = parser.parse_args()
if opts.examples:
    parser.exit(0, """
examples:
sesame.py m32
    : resolve position of m32 in decimal degrees
sesame.py --sexa m54
    : resolve position of m54 in sexagesimal
sesame.py --all crab
    : display full information about object 'crab'
""")
#default object name                  
object = "m1"
#find object name
if len(args) == 0:
     parser.print_help()
     parser.error("no objectname provided")
else:
    object = args[0]

#call sesame service
result = acr.cds.sesame.resolve(object)

#format the results
if opts.all:
    print """
Object Name: %(OName)s, Type: %(OType)s
RA,DEC: %(ra)f,%(dec)f
Sexagesimal: %(posStr)s
Aliases: %(aliases)s""" %  result
elif opts.decimal :
    print "%(ra)f %(dec)f" % result
else:
    print "%(posStr)s" % result
