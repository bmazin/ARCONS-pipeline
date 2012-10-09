from RO import StringUtil
import xmlrpclib
import sys
import os
import optparse

for line in open("mkstds-input.txt"):
    object = "hd"+line[:6]

    #connect to acr
    fname = os.path.expanduser("~/.astrogrid-desktop")
    assert os.path.exists(fname),  'No AR running: Please start your AR and rerun this script'
    prefix = file(fname).next().rstrip()
    acr = xmlrpclib.Server(prefix + "xmlrpc")

    result = acr.cds.sesame.resolve(object)
    newLine=result['posStr']+" "+line
    print newLine.strip()
