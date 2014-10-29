#! /usr/bin/env python 
import sys, getopt,os 
from scipy import stats as stats
import numpy as np 

import Nio 
import myfpzip


def usage():
   print '\n Implements quality metrics tests from the HPDC paper \n'
   print ' fastjuice.py '
   print '   -h              : Prints out this usage message'
   print '   -f, --orig      : Original data file'
   print '   --label         : The label to use for the printout'
   print '   --recon         : Reconstrucetd data file'
   print '   -m $method      : Applies a compression algorithm and creates the reconstructed dataset'
   print '                     Currently supported values of $method are:'
   print '                         grib2, fpzip32, fpzip24, fpzip16, netcdf4'
   print '   --genfile       : If set to be true then will generate an output file with the reconstructed data'
   print '   --max           : Calculates the maximum value in field'
   print '   --min           : Calculates the minimum value in field'
   print '   --std           : Calculates the standard deviation in field'
   print '   --mean,--avg    : Calculates the mean in field'
   print '   --corr          : Calculates the Pearson Correlation'
   print '   --rmse          : Calculates the average error'
   print '   --maxnorm       : Calculates the maximum error'
   print '   --maxnormens    : Calculates the RMSE ensemble'
   print '   --rmsz          : Calculates the RMS-Z score'
   print '   --rmszens       : Calcualtes the RMS-Z ensemble'
   print '   --bias          : Calculates the Bias test'
   print '   --prec          : Supply a precision bit for fpzip'
   print '   '
   print 'Version 0.1'

def notSupportedQuality():
   print ' Sorry the requested quality metric is not yet supported'
def notSupportedGenerate():
   print ' Sorry the generation of reconstruted data files is not yet supported'

def main(argv):
  
  delim = ' : '
  space = ' '
  reqmeth  = []
  setmeth  = False
  label = 'orig'
  SetReconFile = False
  com_size={}
  recon_data={}
  label_temp={}

  # Logical to determine if a file with datafile
  # with reconstructed data is generated
  GenFile   = False

  # Logicals for all the different quality metrics 
  CalcMin        = False
  CalcMax        = False
  CalcStd        = False
  CalcMean       = False
  CalcCorr       = False
  CalcRmse       = False
  CalcMaxnorm    = False
  CalcMaxnormens = False
  CalcRmsz       = False
  CalcRmszens    = False
  CalcBias       = False 
  
  # The start and end time index on which to perform the calculation  
  tstart    = 1
  tend      = 1

  # Pull in and analyze the command line arguements
  try: 
     opts, args = getopt.getopt(argv,"hf:m:",["orig=","recon=","genfile", "label=",\
          "max","min","std","mean","avg","corr","rmse","maxnorm","maxnormens","rmsz","rmszens","bias","prec="])
  except getopt.GetoptError:
     usage()
     sys.exit(2)
  for opt, arg in opts:
     if opt == "-h":
        usage()
        sys.exit()
     elif opt in ("-f","--orig"):
        origfname = arg
     elif opt in ("--recon"):
        reconfname = arg
        SetReconFile = True
     elif opt in ("-m"):
        reqmeth.append(arg)
        setmeth = True
     elif opt in ("--label"):
        label = arg
     elif opt in ("--genfile"):
        GenFile = True
     elif opt in ("--max"):
	CalcMax = True
     elif opt in ("--min"):
        CalcMin = True
     elif opt in ("--std"):
        CalcStd = True
     elif opt in ("--mean","--avg"):
        CalcMean = True
     elif opt in ("--corr"):
        CalcCorr  = True
     elif opt in ("--rmse"):
        CalcRmse  = True
     elif opt in ("--maxnorm"):
        CalcMaxnorm  = True
     elif opt in ("--maxnormens"):
        CalcMaxnormens  = True
     elif opt in ("--rmsz"):
        CalcRmsz  = True
     elif opt in ("--rmszens"):
        CalcRmszens  = True
     elif opt in ("--bias"):
        CalcBias  = True
     elif opt in ("--prec"):
        precision = arg
        prec_int=int(precision)


  #
  # Open up the original data file
  #
  if(os.path.isfile(origfname)):
     fid = Nio.open_file(origfname,"r")
  else:
     print 'file: ',origfname,' Not found'
     sys.exit(2)
  #
  # Open up the reconstructed data file
  #
  if SetReconFile:
     if(os.path.isfile(reconfname)):
           rfid = Nio.open_file(reconfname,"r")
     else:
	if GenFile: 
           rfid = Nio.open_file(reconfname,"w")
	else:
           print 'file: ',reconfname,' Not found'
           sys.exit(2)
  
   
  # Pull out a list of all the timeseries fields on the file
  otimeSeries  = fid.variables
  odimensions  = fid.dimensions
  if SetReconFile and not GenFile:
    rtimeSeries = rfid.variables
    rdimensions = rfid.dimensions

  # TODO: Need to check to see if the varible in the original and reconstructed file match 
  #
  #  MatchDictionary(otimeSeries,rTimeSeries)
  #
  CalcQualityMetrics = CalcMin or CalcMax or CalcMean or CalcStd or CalcCorr or CalcRmse or \
		       CalcMaxnorm or CalcMaxnormens or CalcRmsz or CalcRmszens or CalcBias

  if CalcQualityMetrics:

     for k,v in otimeSeries.iteritems():
        #
        # KLUDGE: The following if test is a rather fragile way to 
        #         limit the calculations to the time-series variables 
        #         only.
        if 'time' in v.dimensions and v.rank > 1 and v.typecode() == 'f':
           dp=0
           orig = v[:]
           #
           # Pull out the reconstructed field
	   #
	   if SetReconFile:
              recon = rtimeSeries[k][:]
           for it in range(tstart,tend+1): 
             if CalcMin:
               min = np.min(orig[it])
               strtmp = k + delim + label + space + 'min' + delim + '{0:9e}'.format(min)
               print strtmp
	     if CalcMax:
	       max = np.max(orig[it])
               strtmp = k + delim + label + space +  'max' + delim + '{0:9e}'.format(max)
               print strtmp
	     if CalcMean:
	       mean = np.mean(orig[it],dtype=np.float64)
               strtmp = k + delim + label + space + 'avg' + delim + '{0:9e}'.format(mean)
               print strtmp
	     if CalcStd:
               std = np.std(orig[it],dtype=np.float64)
               strtmp = k + delim + label + space + 'stddev' + delim + '{0:9e}'.format(std)
               print strtmp
             if CalcMaxnorm:
               if SetReconFile:
		 rng = np.max(orig[it]) - np.min(orig[it])
		 if rng < threshold:
  		   maxnormm = 0.0 
                 else:
		   maxnorm = np.max(orig[it]-recon[it])/rng
		 strtmp = k + delim + label + space + 'maxnorm' + delim + '{0:9e}'.format(maxnorm)
		 print strtmp
             if CalcRmse:
                 orig_size=1
                 for j in orig[it].shape:
                    orig_size=orig_size*j
                 sumsqr=np.sum(np.square(orig[it]-recon[i]))
                 orig_range=np.max(orig[it]-recon[i])
                 if abs(orig_range) < 1e-12:
                   rmse=0.0
                 else:
                   rmse=np.sqrt(sumsqr)/orig_size/orig_range
		 strtmp = k + delim + label_temp[i] + space + 'rmse' + delim + '{0:9e}'.format(rmse)
                 print strtmp
	     if CalcCorr:
                 orig_1d=orig[it].ravel()
                 recon_1d=recon[i].ravel()
                 if abs(np.std(orig[it])) < 1e-12 or abs(np.std(recon[i])) < 1e-12:
                   corr=[0.0, 0.0]
                 else:
                   corr=stats.pearsonr(orig_1d,recon_1d)
		 strtmp = k + delim + label_temp[i] + space + 'corr' + delim + '{0:9e}'.format(corr[0])
                 print strtmp
	     if CalcMaxnormens or CalcRmsz or CalcRmszens or CalcBias:
	       # Perform calculation on the corresponding reconstructed field 
	       notSupportedQuality()


	       if CalcMaxnormens or CalcRmsz or CalcRmszens or CalcBias:
		 # Perform calculation on the corresponding reconstructed field 
		 notSupportedQuality()
  elif GenFile:
     for i in reqmeth:
         if i == 'fpzip' :
            com_size[i],recon[i]=myfpzip.my_python_fpzip(orig[it],prec_int,dp)
            label_temp[i]=i+'{0}'.format(prec_int)
         else:
            notSupportedGenerate()
   
  fid.close()

if __name__ == "__main__":
    main(sys.argv[1:])
