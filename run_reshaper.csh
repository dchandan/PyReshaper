#!/bin/csh
#BSUB -n 1 
#BSUB -q geyser
#BSUB -N 
#BSUB -W 11:00
#BSUB -P STDD0002


setenv CASE b.e12.B1850C5CN.ne30_g16.init.ch.027
setenv DATAROOT /glade/u/tdd/asap/data/
setenv SUBDIR /atm/hist/
setenv OUTDIR /glade/scratch/mickelso/Data/CAMSE-1deg
setenv COMPID "cam.h0"
setenv STARTYR 1
setenv ENDYR 10
setenv METAFILE camse-metainfo.json
setenv FORMAT 'NetCDF4Classic'
setenv COMPRESSLEVEL 1

setenv LID "`date +%y%m%d-%H%M%S`"
setenv HOST "`hostname`"
setenv LOGFILE slice2series.pynio.$HOST.CAMSE-1deg.$LID.log

rm -rf $OUTDIR

module load python
module load all-python-libs

./reShaper.py > & $LOGFILE
#strace -tt -v -T ./reShaper.py > & $LOGFILE


