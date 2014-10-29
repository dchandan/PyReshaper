#! /usr/bin/env python 

import numpy
import sys,os,shutil

import netCDF4 as nc4
 
import glob
import math
import json
import time

import parUtils as par
from collections import OrderedDict

# The blocksize of the GLADE filesystem
bsize = float(4*1024*1024)

#Debug   = True
Debug   = False
WriteTI = True
WriteTV = True

# Pull in environment variables
case     = os.environ['CASE']
dataroot = os.environ['DATAROOT']
subdir   = os.environ['SUBDIR']
outdir   = os.environ['OUTDIR']
model    = os.environ['COMPID']
startyr  = int(os.environ['STARTYR'])
endyr    = int(os.environ['ENDYR'])
metafile = os.environ['METAFILE']

#
# Read in the JSON file which describes the 
# Metadata present in the input file
#
fd = open(metafile)
metainfo = json.load(fd)
mkeysTimeInvarient = metainfo['TimeInvarient']
mkeysTimeVarient   = metainfo['TimeVarient']

# Strings for the input years
years = ['{0:=04d}'.format(i) for i in range(startyr,endyr+1)]

#print 'YEARS: ',years
#print 'CASE: ',case 

# Disable the default PreFill option on NetCDF
#opt = Nio.options()
#opt.PreFill = False
#opt.HeaderReserveSpace = 1024*1024

# Get the rank 
rank = par.GetRank()

timerkey = ['openInput','sortVar','openOutput','metaTIcreate','metaTI','metaTVcreate','metaTV','timeseries']
elapseTime = {}
startTime={}
endTime={}
for k in timerkey:
  elapseTime[k]=0.0


# open all the files for reading
startTime['openInput']=time.time()
ifiles = []
ifilenames = []
for y in years: 
   pattern = y + '-*.nc'
   fileglob = dataroot + case + subdir + case + '.' + model + '.' + pattern
   # print 'fileglob: ',fileglob
   ifilenames = glob.glob(fileglob)
   for iname in ifilenames:
       ifiles.append(nc4.Dataset(iname,"r"))
elapseTime['openInput'] = elapseTime['openInput'] + (time.time() - startTime['openInput'])
#
# total number of input files
#
nslices = len(ifiles) 

idimensions = ifiles[0].dimensions
iattributes = ifiles[0].__dict__

basename = 'foo.nc'


# Break up the dictionary of variables into three different types
#
#   TimeInvarient:  Those variables which do not change overtime 
#
#   TimeVarient:    Varibles that change over tiem but should not have their own timeseries file	
#
#   TimeSeries:	    Variables that should have their own timeseries file
#  
startTime['sortVar']=time.time()
mvarsTimeInvarient={}
mvarsTimeVarient={}
varsTimeSeries={}
temp = ifiles[0].variables
area = temp["area"]
for k,v in temp.iteritems():
    if k in mkeysTimeInvarient:
        # print v.rank
        #if v.rank > 0:
           # print 'rank: ',v.rank
           # bar = v.get_value()
	   # print bar
	   # mvarsTimeInvarient[k] = v
        #else:
	mvarsTimeInvarient[k] = v
    elif k in mkeysTimeVarient:
        mvarsTimeVarient[k] = temp[k]
    else:
        varsTimeSeries[k] = temp[k] 

varsTimeSeries = OrderedDict(sorted(varsTimeSeries.items(), key=lambda t: t[0]))

elapseTime['sortVar'] = elapseTime['sortVar'] + (time.time() - startTime['sortVar'])

if rank == 0:
  print varsTimeSeries.keys()
  print 'Total number of timeseries files : ',len(varsTimeSeries)

lvars = par.Partition(varsTimeSeries)
lvarnames = []
for k,v in lvars.iteritems():
   lvarnames.append(k)

# limit the number of variables created to just U10 for debugging purposes
if Debug:
   lvarnames = ['V','U']
   lvars = {}
   for k in lvarnames:
      lvars[k] = varsTimeSeries[k]


print 'IAM: ',rank,' processing variables: ',lvarnames
#sys.exit()

# create a subdirectory to hold the output files 
if par.IsMaster(): 
   if not os.path.exists(outdir) :
       os.mkdir(outdir)
par.Sync()


vars={}
mvars={}

request={}
requestTI={}
requestTV={}

actual={}
actualTI={}
actualTV={}


ofiles = {}
# iterate over all timeseries output variables
for k,v in lvars.iteritems():
   oname = outdir + '/' + k + '.' + basename 
   # print sys.getsizeof(v[:])
   # print v[:].nbytes
   # print float(v[:].nbytes)/float(4*1024*1024)
   request[k] = 0.0
   actual[k]  = 0.0
	
  
   if 'time' in v.dimensions and not k in mkeysTimeVarient:
      

      print 'IAM: ',rank, ' output filename: ', oname
      if os.path.exists(oname):
         os.unlink(oname)
      startTime['openOutput']=time.time()
      ofile = nc4.Dataset(oname,"w",format='NETCDF3_CLASSIC')
      for ka,va in iattributes.iteritems():
        setattr(ofile,ka,va)
      for k2,v2 in idimensions.iteritems():
         if k2 == 'time':
            # properly setup the unlimited dimension
            ofile.createDimension(k2,None)
         else:
            ofile.createDimension(k2,len(v2))
         # print v
      elapseTime['openOutput'] = elapseTime['openOutput'] + (time.time() - startTime['openOutput'])
      #
      # Create and Write the TimeInvarient meta variables
      #
      startTime['metaTIcreate']=time.time()
      if WriteTI:
          for mk,mv in mvarsTimeInvarient.iteritems():
	      ofile.createVariable(mk,mv.dtype,mv.dimensions)
          # print 'variable: ', mk
          # print 'rank: ', mv
      elapseTime['metaTIcreate'] = elapseTime['metaTIcreate'] + (time.time() - startTime['metaTIcreate'])

      startTime['metaTI']=time.time()
      if WriteTI:
          for mk,mv in mvarsTimeInvarient.iteritems():
             for ka,va in mv.__dict__.iteritems():
                setattr(ofile.variables[mk],ka,va)
	     if mv.ndim > 0:
	        ofile.variables[mk][:] = mv[:]
                # request[mk] = request[mk] + mv[:].nbytes
	        # actual[mk]  = actual[mk] + bsize * math.ceil(mv[:].nbytes/bsize)
             else:
	        ofile.variables[mk] = mv.getValue()
	        # request[mk]         = mv.nbytes
      elapseTime['metaTI'] = elapseTime['metaTI'] + (time.time() - startTime['metaTI'])

      startTime['metaTVcreate']=time.time()
      #
      # Create the Time Varient meta variables 
      #
      for mk,mv in mvarsTimeVarient.iteritems():
          mvars[mk]   = ofile.createVariable(mk,mv.dtype,mv.dimensions)
          
	  # ofile.variables[mk].varAttName = mv.varAttName
      elapseTime['metaTVcreate'] = elapseTime['metaTVcreate'] + (time.time() - startTime['metaTVcreate'])
      #
      # Create the timeseries variable
      #
      vars[k]   = ofile.createVariable(k,v.dtype,v.dimensions)
      ofiles[k] = ofile

# print 'Finished opening files'
for k,v in ofiles.iteritems():
   print 'IAM: ', rank, ' Processing variable: ',k
   temp = v.variables[k]
#   for ka,va in ifiles[0].variables[k].__dict__.iteritems():
#     setattr(temp,ka,va)
   for ka in ifiles[0].variables[k].ncattrs():
     if ka != '_FillValue':
       va = ifiles[0].variables[k].__getattribute__(ka)
       setattr(temp,ka,va)
   vrank = ifiles[0].variables[k].ndim
   if vrank > 0: 
       idx=0
       for j in ifiles: 
	   #
           # Write out the TimeVarient metadata
	   #
           startTime['metaTV'] = time.time()
           if WriteTV:
	       for mk,mv in mvars.iteritems():
	           temp2 = v.variables[mk]
#                   if idx == 0:
#                     for ka,va in j.variables[mk].__dict__.iteritems():
#                       setattr(temp2,ka,va)
#                     for ka in j.variables[mk].ncattrs():
#                       if ka != '_FillValue':
#                         va = j.variables[mk].__getattribute__(ka)
#                         setattr(temp2,ka,va)
		   amount = j.variables[mk][:].nbytes
                   temp2[idx] = j.variables[mk][:]  
		   request[k] = request[k] + amount
                   actual[k]  = actual[k] + bsize * math.ceil(amount/bsize) 
           elapseTime['metaTV'] = elapseTime['metaTV'] + (time.time() - startTime['metaTV'])
	   #
	   # Write out the Timeseries variables
	   #
           startTime['timeseries'] = time.time()
           temp[idx]  = j.variables[k][:]
           amount = j.variables[k][:].nbytes
           request[k] = request[k] + amount
           actual[k]  = actual[k] + bsize * math.ceil(amount/bsize) 
           elapseTime['timeseries'] = elapseTime['timeseries'] + (time.time() - startTime['timeseries'])
	   idx=idx+1
   #
   # close the file 
   #
   v.close()

#
# Perform some performance analysis
#
#actualTotal=0.0
#for k,v in actual.iteritems():
#   actualTotal = actualTotal + actual[k]   
actualTotal  = par.Sum(actual)

#requestTotal=0.0
#for k,v in request.iteritems():
#   requestTotal = requestTotal + request[k]   
requestTotal = par.Sum(request)

openi        = par.Max(elapseTime['openInput'])
openo        = par.Max(elapseTime['openOutput'])
sortvar      = par.Max(elapseTime['sortVar'])
metaTIcreate = par.Max(elapseTime['metaTIcreate'])
metaTI       = par.Max(elapseTime['metaTI'])
metaTVcreate = par.Max(elapseTime['metaTVcreate'])
metaTV       = par.Max(elapseTime['metaTV'])
timeseries   = par.Max(elapseTime['timeseries'])
# Sync up before printing results 
par.Sync()
#

#  Print performance analysis
#
if par.IsMaster():
   print 'actual: ',actual
   print 'request: ',request
   print '==============================================='
   print '       Performance data '
   print '==============================================='
   print ' Opening Inputfiles: ',openi
   print ' Sorting out variables: ',sortvar
   print ' Opening Outputfiles: ',openo
   print ' Create TimeInvarient metadata: ',metaTIcreate
   print ' Write  TimeInvarient metadata: ',metaTI
   print ' Create TimeVarient metadata: ', metaTVcreate
   print ' Write TimeVarient metadata: ', metaTV
   print ' Writing Timeseries data: ',timeseries
   print '   Actual MB: ',actualTotal/float(1024*1024)
   print 'Requested MB: ',requestTotal/float(1024*1024)
   print 'Min useful fraction: ',requestTotal/actualTotal

