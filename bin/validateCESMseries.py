#!/usr/bin/env python

"""
This script aids in the validation of data for the timeseries files
created using "slice2series". It uses multiple parallel processes to
perform the validation using the multiprocessing module. This module 
is used instead of mpi4py because it provides convenient communincation
vectors (such as Queues) between processes that are suitable to the
task at hand. 
"""

import os, sys
import os.path as ospath
import numpy as np
import multiprocessing
from multiprocessing import Process, Queue
import logging
import argparse
import glob
import gc

try:
    import Nio
    nio_available = True
except ImportError:
    nio_available = False

try:
    import netCDF4
    netcdf_available = True
except ImportError:
    netcdf_available = False

from PyCESM.case.CESMCase import CESMCase



def get_unlimited_dimension(series_file):
    """
    DESCRIPTION
        Get the name of the unlimited dimension in the timeseries files.
        Asusmes that:
            (i) there is only one unlimited dimension
            (ii) that the unlimited dimension is the same in all timeseries files
    ARGUMENTS
        series_file - One timeseries file
    RETURNS 
        Name of the unlimited dimension
    """
    ncid = Nio.open_file(series_file, "r")

    unlimitd_dim = None
    for dim in ncid.dimensions:
        if ncid.unlimited(dim): 
            unlimited_dim = dim
            break
    ncid.close()
    if (unlimited_dim == None):
        raise ValueError("No unlimited dimension found in timeseries file {0}".format(series_file))
    return unlimited_dim



def worker(work_queue, return_queue, unlim_dim, slice_files, backend):
    """
    DESCRIPTION
        The main code that will be executed by all workers processes. 
    ARGUMENTS
        work_queue  - A multiprocessing.Queue object from which the function will
                      extract a timeseries filename to test. 
        return_queue- A queue in which to return if the validation succeeded
        unlim_dim   - Name of the unlimited dimension in the timeseries files
        slice_files - List of time-slice files that were used in the construction of
                      the timeseries files. 
        backend     - Backend netcdf interface library
    RETURNS

    """
    logger = multiprocessing.get_logger()  # get the logger for the multiprocessing module

    opened_files = []
    for fname in slice_files:
        try:
            if (backend == "nio"):
                opened_files.append(Nio.open_file(fname, "r"))
            else:
                opened_files.append(netCDF4.Dataset(fname, "r"))
        except:
            raise IOError("Error while opening {0}".format(fname))
    num_vars_processed = 0
    while True:
        # Infinite loop. The worker stops when the poison pill is encountered
        series_file = work_queue.get()  # Get a piece of item
        
        # Poison pill
        if series_file == "STOP": break
        
        logger.info(series_file)

        all_okay   = True
        num_errors = 0
        result     = {}
        error_vars = set()

        if (backend == "nio"):
            ncfile = Nio.open_file(series_file, "r")
        else:
            ncfile = netCDF4.Dataset(series_file, "r")

        time_indep_vars = []   # Time dependent variables
        time_dep_vars   = []   # Time independent variables
        variables = ncfile.variables
        for var_name in variables.keys():
            var = variables[var_name]
            if unlim_dim not in var.dimensions:
                time_indep_vars.append(var_name)
            else:
                time_dep_vars.append(var_name)

        # Loop through all the time dependent variables
        for var_name in time_dep_vars:
            var = variables[var_name]

            # Determining the type of the variable
            if (backend == "nio"):
                var_dtype = np.dtype(var.typecode())
            else:
                var_dtype = var.dtype

            # Checking 'non-string' variables only
            if (var_dtype != np.dtype('S1')):
                if len(var.dimensions) != 1: # not scalar
                    for i, _slice in enumerate(opened_files):
                        try:
                            assert(np.allclose(var[i,:], _slice.variables[var_name][0,:]))
                        except AssertionError:
                            all_okay = False
                            num_errors += 1
                            error_vars.add(var_name)
                else: # scalar
                    for i, _slice in enumerate(opened_files):
                        try:
                            assert(np.allclose(var[i], _slice.variables[var_name][:]))
                        except AssertionError:
                            all_okay = False
                            num_errors += 1
                            error_vars.add(var_name)


        ncfile.close()
        del ncfile
        num_vars_processed += 1

        result["pass"] = all_okay
        result["num_errors"] = num_errors
        result["error_vars"] = list(error_vars)
        result["filename"]   = series_file
        return_queue.put(result)
    return


def subcommand_range(args):
    case        = CESMCase(args.case)

    total_years = args.end - args.start + 1
    years       = range(args.start, args.end + 1)
    
    # A mapping between model names and file corresponding file names
    mtypes    = {"atm":"cam2", "lnd":"clm2", "ocn":"pop", "ice":"cice"}
    freqtypes = {"atm":"h0", "lnd":"h0", "ocn":"h", "ice":"h"}


    #comp_direc = osp.join("/scratch/p/peltier/dchandan/ctest", model, "hist")
    comp_direc = ospath.join(case.DOUT_S_ROOT, args.model, "hist")
    # os.chdir(comp_direc)

    # Generating the list of files that need to be worked on
    slice_files = []
    for year in years:
        pattern = "{0}/{1}.{2}.{3}.{4:04d}*.nc".format(comp_direc,
                                                      args.case, 
                                                      mtypes[args.model], 
                                                      freqtypes[args.model], 
                                                      year)
        slice_files.extend(glob.glob(pattern))

    slice_files.sort()

    print("+--------------------------------------+")
    print("| Number of files to operate upon: {0:3d} |".format(len(slice_files)))
    print("+--------------------------------------+")

    main(slice_files, args.prefix, args.n, args.backend)


def subcommand_files(args):
    main(args.slice_files, args.prefix, args.n, args.backend)


def main(slice_files, prefix, nprocs, backend):
    num_slice_files = len(slice_files)

    print("Configurations for validator script:")
    print("Filename prefix : {0}".format(prefix))
    print("Netcdf backend  : {0}".format(backend))
    print("Number of time-slice files found : {0}".format(num_slice_files))


    print("Generating collection of timeseries files ... ")
    series_files = glob.glob(prefix + "*")

    # Query the unlimited dimension name on one of the files
    unlim_dim = get_unlimited_dimension(series_files[0])

    # These are the Queues we'll use to send and receive data from the worker
    # processes
    work_queue   = Queue()
    return_queue = Queue()

    # Logger that takes care of pretty printing to the screen from parallel processes
    logger = multiprocessing.log_to_stderr()
    logger.setLevel(logging.INFO)

    # Putting 'jobs', i.e. timeseries files into the work queue
    for fname in series_files: work_queue.put(fname)

    processes = []  # List of parallel processes
    for i in range(nprocs):
        p = Process(target=worker, args=(work_queue, return_queue,
                                         unlim_dim, slice_files, backend))
        p.start()
        processes.append(p)
        work_queue.put('STOP')  # Putting a poison pill for the workers

    # Wait for jobs to complete
    for p in processes: p.join()

    # Checking results >>>>>>>>>
    files_pass = 0  # Counter for number of timeseries files that passed validation
    files_fail = 0  # Counter for number of timeseries files that failed validation
    
    for i in xrange(len(series_files)):
        result = return_queue.get()
        if result['pass'] == False:
            print("FAIL: Variables {0} in file {1}".format(result['error_vars'], result['filename']))
            files_fail += 1
        else:
            files_pass += 1

    print("-"*43)
    print("Number of timeseries files validated: {0}".format(len(series_files)))
    print("Number of files that passed test    : {0}".format(files_pass))
    print("Number of files that failed test    : {0}".format(files_fail))
    print("-"*43)
    # <<<<<<<<< Checking results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="validateCESMseries.py")
    parser.add_argument('--prefix',  type=str, help='filename prefix', required=True)

    backend_choices = []
    if netcdf_available: backend_choices.append("netcdf")
    if nio_available: backend_choices.append("nio")
    if not (nio_available or netcdf_available):
        print("ERROR: Neither the PyNIO nor the netCDF4 backend is available.")
        sys.exit(-1)

    parser.add_argument('--backend', type=str, choices=backend_choices, 
                        default=backend_choices[0], help='netcdf backend to use')

    subparsers = parser.add_subparsers(help='help for subcommand')

    parser_a = subparsers.add_parser('files', help='files help')
    parser_a.add_argument('slice_files', nargs='+', type=str, help="All time slice files used "
                                       "in the creation of the timeseries files")
    parser_a.set_defaults(func=subcommand_files)
    parser_b = subparsers.add_parser('range', help='help for range')
    parser_b.add_argument('case',  type=str, help='CESM case name')
    parser_b.add_argument('model', type=str, choices=["atm", "lnd", "ocn", "ice"], help='models component to process')
    parser_b.add_argument('start', type=int, help='start year')
    parser_b.add_argument('end',   type=int, help='end year')
    parser_b.set_defaults(func=subcommand_range)

    parser.add_argument('-n',    type=int, default=8, help='number of parallel processes')

    args = parser.parse_args()
    args.func(args)

