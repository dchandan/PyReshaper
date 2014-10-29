#!/usr/bin/env python

import sys
import os
import string

import Nio
import glob
import math
import json
import time

import parUtils as par


#===============================================================================
# Check to see if we need to append for this date stamp
#===============================================================================
def time_series_append_check(existing_stamps, n_start):

    # Parse the new time stamp
    new_start = string.split(n_start)
    new_start_len = len(new_start)
    new_start_year = (new_start[0:3])
    if (new_start_len > 4):
        new_start_mon = (new_start[4:5])
    if(new_start_len > 6):
        new_start_day = (new_start[6:7])
    if(new_start_len > 8):
        new_start_time = (new_start[8:11])

    new_start_stamp = n_start

    for stamp in existing_stamps:
        # Parse the existing time stamp and pull out year, mon, etc for start and end
        start_end = string.split(stamp, '-')
        start = string.split(start_end[0])
        start_len = len(start)
        end = string.split(start_end[1])
        end_len = len(end)

        start_year = (start[0:3])
        if (start_len > 4):
            start_mon = (start[4:5])
        if(start_len > 6):
            start_day = (start[6:7])
        if(start_len > 8):
            start_time = (start[8:11])

        end_year = (end[0:3])
        if (end_len > 4):
            end_mon = (end[4:5])
        if(end_len > 6):
            end_day = (end[6:7])
        if(end_len > 8):
            end_time = (end[8:11])

        # Compare to see if the start year for the new stamp exists in a
        # preveously computed timeseries file
        if (new_start_year <= end_year) and (new_start_year >= start_year):
            new_start_stamp = string.split(stamp, '-')[0]

    return new_start_stamp


#===============================================================================
#
# Make a list of time series date stamps from the time series that have already
# been computed
#
#===============================================================================
def existing_series_list(output_directory, root_name, start_time, end_time, output_time_option):
    # Figure out the time series set steps
    # A time series set includes all of the variable files within a set time stamp range

    # Parse and set start times and end times that the user wants to convert to time series
    start = string.split(start_time, '-')
    end = string.split(end_time, '-')

    # Parse the start stamp of the history files the user wants to convert to time series
    start_len = len(start)
    start_year = int(start[0])
    if (start_len >= 2):
        start_month = int(start[1])
    else:
        start_month = -99
    if (start_len >= 3):
        start_day = int(start[2])
    else:
        start_day = -99
    if (start_len >= 4):
        start_hour = int(start[3])
    else:
        start_hour = -99

    # Parse the end stamp of the history files the user wants to convert to time series
    end_len = len(end)
    end_year = int(end[0])
    if (end_len >= 2):
        end_month = int(end[1])
    else:
        end_month = -99
    if (end_len >= 3):
        end_day = int(start[2])
    else:
        end_day = -99
    if (end_len >= 4):
        end_hour = int(start[3])
    else:
        end_hour = -99

    # Make a time stamp list to list the existing timeseries date ranges that have
    # already been computed.    This will help us determine if we need to append to
    # an existing timeseries, over write time slices, or create a new set of files.
    # Ex timeseries file name b40.20th.track1.1deg.006.cam2.h0.T.185001-185912.nc
    pattern = output_directory + root_name + '*.nc'
    dir_list = glob.glob(pattern)
    existing_time_stamps_temp = []
    for fn in dir_list:
        fn_p = string.split(fn, '.')  # Split the string up by '.'
        existing_time_stamps_temp.append(fn_p[len(fn_p) - 2])

    #Remove the duplicates by adding to a set, then delete old list
    existing_time_stamps = set(existing_time_stamps_temp)
    del existing_time_stamps_temp

    #Make a list of all time series date stamp ranges that will be converted
    new_time_stamps = []
    n_ts = 0
    if (output_time_option == 'nyears') or (output_time_option == 'nyear'):
        i = start_year
        while i <= end_year:
            # Try to figure out the time stamp based on what *logically* would come next
            if ((i % output_time_n) != 0):
                increment = output_time_n - (i % output_time_n)
            else:
                increment = output_time_n

            guess_stamp = str(i).zfill(4)
            if (start_month != -99):
                if (i == start_year):
                    guess_stamp += str(start_month).zfill(2)
                else:
                    guess_stamp += '01'
            if(start_day != -99):
                if (i == start_year):
                    guess_stamp += str(start_day).zfill(2)
                else:
                    guess_stamp += '01'
            if(start_hour != -99):
                if (i == start_year):
                    guess_stamp += str(start_hour).zfill(4)
                else:
                    guess_stamp += '0001'

            # Check the time stamp guess aginst previously created time series to see
            # if we'll need to append or create a new file.    If appending, reset the
            # time stamp to reflect the start date for the previously started
            # timeseires that already exists.
            if (len(existing_time_stamps) > 0):
                result = time_series_append_check(existing_time_stamps, guess_stamp)
            else:
                result = guess_stamp

            new_end_year = (i + increment) - 1

            #Append the new time stamp together
            if (new_end_year >= end_year):
                full_end_stamp = str(end_year).zfill(4)
                if (end_month != -99):
                    full_end_stamp += str(end_month).zfill(2)
                if(end_day != -99):
                    full_end_stamp += str(end_day).zfill(2)
                if (end_hour != -99):
                    full_end_stamp += str(end_hour).zfill(4)
                full_stamp = string.join((result, "-", full_end_stamp), "")
            else:
                full_stamp = string.join((result, "-", str(new_end_year).zfill(4), "12"), "")

            new_time_stamps.append(full_stamp)
            i += increment

    return new_time_stamps


#===============================================================================
# Check that the time period between time slice files are even
#===============================================================================


def check_time_spacing(time_fields, start_month):
    day_to_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    previous = time_fields[0]
    incr = 0
    i = 1
    m = start_month
    if (m == 12):
        m = 0
    # Loop through the times and get the differnece between the previous and
    # next times and see if the values are correct
    while (i < len(time_fields)):
        next = time_fields[i]
        incr = next - previous
        okay = 1
        # Compare agaist the correct increment based on the input_time_period
        # that was passed to the reshaper code
        if (input_time_period == 'year'):
            if (incr != 365):
                okay = 0
        elif (input_time_period == 'mon'):
            if (incr != day_to_month[m]):
                okay = 0
        elif(input_time_period == 'day'):
            if (incr != 1):
                okay = 0
        #### NEED TO ADD HOURLY INCREMENT !!!!! ######
        # Check to see if the icrement between the steps were expected
        if (okay == 0):
            print "Incorrect time period between the files.    Exiting.",
            print str(next) + str(previous) + str(incr) + str(day_to_month[m])
            sys.exit(1)
        # increase indexes and set a new previous value
        i += 1
        # Check to see if the end of the year has been reached.    If so, reset the
        # day_to_month index to point to January.
        if ((i % 12) == 0):
            m = 0
        else:
            m += 1
        previous = next


#===============================================================================
#
# John's original reshaper code
#
#===============================================================================
def convert_slice_to_series(case, dataroot, outdir, ts, metafile):
    # The blocksize of the GLADE filesystem
    bsize = float(4 * 1024 * 1024)

    Debug = True

    WriteTI = True
    WriteTV = True

    #
    # Read in the JSON file which describes the
    # Metadata present in the input file
    #
    fd = open(metafile)
    metainfo = json.load(fd)
    mkeysTimeInvarient = metainfo['TimeInvarient']
    mkeysTimeVarient = metainfo['TimeVarient']

    # Get the start and stop years for this set of time series
    start_end = string.split(ts, '-')
    start = string.split(start_end[0])
    end = string.split(start_end[1])
    startyr = int(start[0][:4])
#    if (len(start_end[0]) > 4):
    start_mon = int(start[0][4:6])
    endyr = int(end[0][:4])
#    if (len(start_end[1]) > 4):
    end_mon = int(end[0][4:6])

    # Strings for the input years
    years = ['{0:=04d}'.format(i) for i in range(startyr, endyr + 1)]

    # Disable the default PreFill option on NetCDF
    opt = Nio.options()
    opt.PreFill = False
    opt.Format = ncformat
    opt.CompressionLevel = compression_level

    # Get the rank
    rank = par.GetRank()

    timerkey = ['openInput', 'sortVar', 'openOutput', 'metaTIcreate', 'metaTI',
                            'metaTVcreate', 'metaTV', 'timeseries']
    elapseTime = {}
    startTime = {}
    endTime = {}
    for k in timerkey:
        elapseTime[k] = 0.0

    # open all the files for reading
    startTime['openInput'] = time.time()
    ifiles = []
    for y in years:
        ifilenames = []
        # Check to see if the month boundaries for the timeseries end on year
        # bounds.    If not, make sure only part of the year is added to the
        # file list
        if (y == str(startyr).zfill(4)) and (start_mon != 1):
            m = start_mon
            n_files = (12 - start_mon) + 1
            while (m <= 12):
                filename = dataroot + case + str(y).zfill(4) + "-" + str(m).zfill(2) + '.nc'
                # double check to see if the file exists on disk, exit if it's not there
                if (os.path.isfile(filename)):
                    ifilenames.append(filename)
                else:
                    print "COULD NOT LOCATE FILE    " + filename + "    EXITING"
                    sys.exit()
                m += 1
        elif (y == str(endyr).zfill(4)) and (end_mon != 12):
            m = 1
            n_files = end_mon
            while (m <= end_mon):
                filename = dataroot + case + str(y).zfill(4) + "-" + str(m).zfill(2) + '.nc'
                # double check to see if the file exists on disk, exit if it's not there
                if (os.path.isfile(filename)):
                    ifilenames.append(filename)
                else:
                    print "COULD NOT LOCATE FILE    " + filename + "    EXITING"
                    sys.exit()
                m += 1
        else:
            n_files = 12  # Should get a full years worth of time slice files
            pattern = y + '-*.nc'
            fileglob = dataroot + case + pattern
            ifilenames = glob.glob(fileglob)

        # Print out input file list for this year
        for fn in ifilenames:
            print "Adding: " + fn
        print '\n'

        #Check to see if the correct number of files were found
        if (len(ifilenames) != n_files):
            print "COULD NOT LOCATE ALL OF THE NEEDED FILES FOR YEAR " + y + "    EXITING"
            sys.exit(1)
        # Open the input files
        for iname in ifilenames:
            ifiles.append(Nio.open_file(iname, "r"))

        # Delete ifilenames to create a new one for the next year in the loop
        del ifilenames

    elapseTime['openInput'] = elapseTime['openInput'] + (time.time() - startTime['openInput'])
    #
    # total number of input files
    #
    nslices = len(ifiles)

    idimensions = ifiles[0].dimensions
    iattributes = ifiles[0].attributes

    # Get the time value for each file and make sure the time periods between the slices are identical
    time_fields = []
    for f in ifiles:
        time_var = f.variables['time']
        for value in time_var.get_value():
            time_fields.append(value)
    check_time_spacing(time_fields, start_mon)

    # Break up the dictionary of variables into three different types
    #
    #     TimeInvarient:    Those variables which do not change overtime
    #
    #     TimeVarient:        Varibles that change over tiem but should not have their own timeseries file
    #
    #     TimeSeries:         Variables that should have their own timeseries file
    #
    startTime['sortVar'] = time.time()
    mvarsTimeInvarient = {}
    mvarsTimeVarient = {}
    varsTimeSeries = {}
    temp = ifiles[0].variables
    for k, v in temp.iteritems():
        if k in mkeysTimeInvarient:
            mvarsTimeInvarient[k] = v
        elif k in mkeysTimeVarient:
            mvarsTimeVarient[k] = temp[k]
        else:
            varsTimeSeries[k] = temp[k]
    elapseTime['sortVar'] = elapseTime['sortVar'] + (time.time() - startTime['sortVar'])

    print 'Total number of timeseries files : ', len(varsTimeSeries)

    lvars = par.Partition(varsTimeSeries)
    lvarnames = []
    for k, v in lvars.iteritems():
        lvarnames.append(k)

    # limit the number of variables created to just U10 for debugging purposes
    if Debug:
        lvarnames = ['V', 'U']
        lvars = {}
        for k in lvarnames:
            lvars[k] = varsTimeSeries[k]

    print 'IAM: ', rank, ' processing variables: ', lvarnames

    # create a subdirectory to hold the output files
    if par.IsMaster():
        if not os.path.exists(outdir):
            os.mkdir(outdir)
    par.Sync()

    vars = {}
    mvars = {}

    request = {}
    requestTI = {}
    requestTV = {}

    actual = {}
    actualTI = {}
    actualTV = {}

    ofiles = {}
    # iterate over all timeseries output variables
    for k, v in lvars.iteritems():
        oname = outdir + '/' + root_name + k + '.' + ts + '.nc'
        request[k] = 0.0
        actual[k] = 0.0

        if 'time' in v.dimensions and not k in mkeysTimeVarient:
            print 'IAM: ', rank, ' output filename: ', oname
            if os.path.exists(oname):
                os.unlink(oname)
            startTime['openOutput'] = time.time()
            ofile = Nio.open_file(oname, "w", options=opt)
            for ka, va in iattributes.iteritems():
                setattr(ofile, ka, va)
            for k2, v2 in idimensions.iteritems():
                if k2 == 'time':
                    # properly setup the unlimited dimension
                    ofile.create_dimension(k2, None)
                else:
                    ofile.create_dimension(k2, v2)
            elapseTime['openOutput'] = elapseTime['openOutput'] + (time.time() - startTime['openOutput'])

            #
            # Create and Write the TimeInvarient meta variables
            #
            startTime['metaTIcreate'] = time.time()
            if WriteTI:
                for mk, mv in mvarsTimeInvarient.iteritems():
                    temp3 = ofile.create_variable(mk, mv.typecode(), mv.dimensions)
                    for ka, va in mv.attributes.iteritems():
                        setattr(temp3, ka, va)
            elapseTime['metaTIcreate'] = elapseTime['metaTIcreate'] + (time.time() - startTime['metaTIcreate'])

            startTime['metaTVcreate'] = time.time()
            #
            # Create the Time Varient meta variables
            #
            for mk, mv in mvarsTimeVarient.iteritems():
                mvars[mk] = ofile.create_variable(mk, mv.typecode(), mv.dimensions)
                for ka, va in mv.attributes.iteritems():
                    setattr(mvars[mk], ka, va)
            elapseTime['metaTVcreate'] = elapseTime['metaTVcreate'] + (time.time() - startTime['metaTVcreate'])
            #
            # Create the timeseries variable
            #
            vars[k] = ofile.create_variable(k, v.typecode(), v.dimensions)
            ofiles[k] = ofile

    # print 'Finished opening files'
    for k, v in ofiles.iteritems():
        print 'IAM: ', rank, ' Processing variable: ', k
        temp = v.variables[k]
        for ka, va in ifiles[0].variables[k].attributes.iteritems():
            setattr(temp, ka, va)
        vrank = ifiles[0].variables[k].rank
        startTime['metaTI'] = time.time()
        if WriteTI:
            for mk, mv in mvarsTimeInvarient.iteritems():
                # print 'Writing TimeInvarient variable: ',mk,' with rank: ',mv.rank
                if mv.rank > 0:
                    v.variables[mk][:] = mv[:]
                else:
                    v.variables[mk].assign_value(mv.get_value())
        elapseTime['metaTI'] = elapseTime['metaTI'] + (time.time() - startTime['metaTI'])
        if vrank > 0:
            idx = 0
            for j in ifiles:
                #
                # Write out the TimeVarient metadata
                #
                startTime['metaTV'] = time.time()
                if WriteTV:
                    for mk, mv in mvars.iteritems():
                        temp2 = v.variables[mk]
                        temp2[idx] = j.variables[mk][:]
                        amount = j.variables[mk][:].nbytes
                        request[k] = request[k] + amount
                        actual[k] = actual[k] + bsize * math.ceil(amount / bsize)
                elapseTime['metaTV'] = elapseTime['metaTV'] + (time.time() - startTime['metaTV'])

                #
                # Write out the Timeseries variables
                #
                startTime['timeseries'] = time.time()
                temp[idx] = j.variables[k][:]
                amount = j.variables[k][:].nbytes
                request[k] = request[k] + amount
                actual[k] = actual[k] + bsize * math.ceil(amount / bsize)
                elapseTime['timeseries'] = elapseTime['timeseries'] + (time.time() - startTime['timeseries'])
                idx = idx + 1
        #
        # close the file
        #
        v.close()

    return actual, request, elapseTime

#########################################
#########################################
#
#     START MAIN PROGRAM
#
#########################################
#########################################

#########################################
# Set all of the control variables
########################################
# Get control variable values from outside source (hard coded for now until the
# namelist to python parser is done)

file_type = 'netcdf'  # NOT NEEDED FOR NOW?
input_directory = \
    '/glade/u/tdd/asap/data/b.e12.B1850C5CN.ne30_g16.init.ch.027/atm/hist/'
output_directory = '/glade/scratch/mickelso/Data/CAMSE-1deg/'
in_time_format = 'yyyy-mm'
out_time_format = 'yyyymm'
start_time = '0004-03'
end_time = '0012-12'
root_name = 'b.e12.B1850C5CN.ne30_g16.init.ch.027.cam.h0.'
input_time_period = 'mon'
output_time_option = 'nyears'
output_time_n = 10
compset = 'B1850C5CN'
grid = 'ne30_g16'
meta_vars = ' '
metafile = 'camse-metainfo.json'
ncformat = 'NetCDF4Classic'
compression_level = 1
##### NEED A VARIABLE FOR LEAP YEARS PRESENT?#####

# Set constant time lengths
month_to_year = 12
### ADD ABILITY FOR LEAP YEARS FOR NEXT 4 VARIABLES?
day_to_year = 365
sixH__to_year = 1460
threeH_to_year = 2920
day_to_month = {'01': 31, '02': 28, '03': 31, '04': 30, '05': 31, '06': 30,
                '07': 31, '08': 31, '09': 30, '10': 31, '11': 30, '12': 31}

#########################################
#
# Open the first file, pull out the variable list and determine if the variable
# is a time series variable, a dimension, or belongs in the once file
#
#########################################

# Open the first time slice file
first_file_name = input_directory + root_name + start_time + ".nc"
first_file = Nio.open_file(first_file_name)

# GET NEEDED FILE INFORMATION HERE


# Close the first time slice file
first_file.close()

#########################################
#
# Figure out the time series set steps
# A time series set includes all of the variable files
# within a set time stamp range
#
#########################################

# Make a time stamp list to list the existing timeseries date ranges that have
# already been computed
new_time_stamps = existing_series_list(output_directory, root_name, start_time,
                                       end_time, output_time_option)

#########################################
#
# Start the Main Loop that will create one timeseries set at a time
# (A time series set is a set a variable time series files that all contain the
# same date range)
#
#########################################

for ts in new_time_stamps:
    # Call the function that creates the time series files and return
    # timing/accounting information
    actual, request, elapseTime = convert_slice_to_series(
        root_name, input_directory, output_directory, ts, metafile)



#########################################
#
# Perform some performance analysis
#
#########################################
actualTotal = par.Sum(actual)

requestTotal = par.Sum(request)

openi = par.Max(elapseTime['openInput'])
openo = par.Max(elapseTime['openOutput'])
sortvar = par.Max(elapseTime['sortVar'])
metaTIcreate = par.Max(elapseTime['metaTIcreate'])
metaTI = par.Max(elapseTime['metaTI'])
metaTVcreate = par.Max(elapseTime['metaTVcreate'])
metaTV = par.Max(elapseTime['metaTV'])
timeseries = par.Max(elapseTime['timeseries'])
# Sync up before printing results
par.Sync()


#    Print performance analysis

if par.IsMaster():
    print 'actual: ', actual
    print 'request: ', request
    print '==============================================='
    print '             Performance data '
    print '==============================================='
    print ' Opening Inputfiles: ', openi
    print ' Sorting out variables: ', sortvar
    print ' Opening Outputfiles: ', openo
    print ' Create TimeInvarient metadata: ', metaTIcreate
    print ' Write    TimeInvarient metadata: ', metaTI
    print ' Create TimeVarient metadata: ', metaTVcreate
    print ' Write TimeVarient metadata: ', metaTV
    print ' Writing Timeseries data: ', timeseries
    print '     Actual MB: ', actualTotal / float(1024 * 1024)
    print 'Requested MB: ', requestTotal / float(1024 * 1024)
    print 'Min useful fraction: ', requestTotal / actualTotal
