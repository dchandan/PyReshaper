#!/usr/bin/env python
"""
This script provides the command-line interface (CLI) to the PyReshaper

This script uses the optparse module for Python 2.6 and the newer argparse
module for Python 2.7.

Copyright 2015, University Corporation for Atmospheric Research
See the LICENSE.txt file for details
"""

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


import os, glob, sys
import os.path as ospath
from pyreshaper.specification import create_specifier
from pyreshaper.reshaper import create_reshaper

import argparse

from asaptools.simplecomm import create_comm, SimpleComm
from PyCESM.case.CESMCase import CESMCase


def cli():
    # Choices for the netCDF backend library
    backend_choices = []
    if netcdf_available: backend_choices.append("netcdf")
    if nio_available: backend_choices.append("nio")
    if not (nio_available or netcdf_available):
        print("ERROR: Neither the PyNIO nor the netCDF4 backend library is installed.")
        sys.exit(-1)

    parser = argparse.ArgumentParser(prog="{0}: Convert CESM slice files to series files".format(__file__))
    parser.add_argument('case',  type=str, help='CESM case name')
    parser.add_argument('model', type=str, choices=["atm", "lnd", "ocn", "ice"], help='models component to process')
    parser.add_argument('start', type=int, help='start year')
    parser.add_argument('end',   type=int, help='end year')
    parser.add_argument("odir",  type=str, action='store', help="Directory to write output files into")

    parser.add_argument('--backend', type=str, choices=backend_choices, 
                    default=backend_choices[0], help='netcdf backend to use')
    parser.add_argument('-m', '--metadata', type=str, nargs="+", default=[],
                      help='Names of a variable to be included in all output '
                           'files. '
                           '[Default: none]')
    parser.add_argument('-d', '--deflate', default=3, type=int, action='store',
                      help='Compression level for the output files. Only used '
                           'when netcdf_format = netcdf4c.'
                           '[Default: 3]')
    parser.add_argument('--serial', default=False,
                      action='store_true', dest='serial',
                      help='Whether to run in serial (True) or parallel (False). '
                           '[Default: False]')
    parser.add_argument('--once', default=False,
                      action='store_true', dest='once',
                      help='Whether to write a "once" file with all metadata. '
                           '[Default: False]')
    parser.add_argument('--skip_existing', default=False,
                      action='store_true', dest='skip_existing',
                      help='Whether to skip time-series generation for '
                           'variables with existing output files. '
                           '[Default: False]')
    parser.add_argument('--overwrite', default=False,
                      action='store_true', dest='overwrite',
                      help='Whether to overwrite existing output files. '
                           '[Default: False]')
    parser.add_argument('--timecode', default=False,
                      action='store_true', dest='timecode',
                      help='Whether to time the code internally. '
                           '[Default: False]')
    parser.add_argument('--preprocess', default=False,
                      action='store_true', dest='preprocess',
                      help='Whether to preprocess the input files for validation purposes. '
                           '[Default: False]')
    parser.add_argument('-v', '--verbosity', default=1, type=int, action='store',
                      help='Verbosity level for level of output.  A value of 0 '
                           'means no output, and a value greater than 0 means '
                           'more output detail. [Default: 1]')
    parser.add_argument('-l', '--limit', default=0, type=int, action='store',
                      help='The limit on the number of time-series files per '
                           'processor to write.  Useful when debugging.  A '
                           'limit of 0 means write all output files.'
                           '[Default: 0]')



    # Parse the CLI options and assemble the Reshaper inputs
    # print parser.parse_args()
    return parser.parse_args()


#==============================================================================
# Main script
#==============================================================================
def main(args):
    simplecomm = create_comm(serial=args.serial)

    # Only the manager process needs to write the configuration to the screen
    # and generate the list of input files.
    if simplecomm.is_manager():
      print("Configuration  >>>>>>>")
      print("  Case name      : {0}".format(args.case))
      print("  Model          : {0}".format(args.model))
      print("  Start year     : {0}".format(args.start))
      print("  End year       : {0}".format(args.end))
      print("")
      print("  Serial         : {0}".format(args.serial))
      print("  Once           : {0}".format(args.once))
      print("  Skip Existing  : {0}".format(args.skip_existing))
      print("  Overwrite      : {0}".format(args.overwrite))
      print("  Time Code      : {0}".format(args.timecode))
      print("  Preprocess     : {0}".format(args.preprocess))
      print("  Backend        : {0}".format(args.backend))
      print("  Compression    : {0}".format(args.deflate))
      print("  Verbosity      : {0}".format(args.verbosity))
      print("<<<<<<<<<<<<<<<<<<<<<<")

      # The CESM case object
      case        = CESMCase(args.case)

      total_years = args.end - args.start + 1
      years       = range(args.start, args.end + 1)
      
      # A mapping between model names and file corresponding file names
      mtypes    = {"atm":"cam2", "lnd":"clm2", "ocn":"pop", "ice":"cice"}
      freqtypes = {"atm":"h0",   "lnd":"h0",   "ocn":"h",   "ice":"h"}


      #comp_direc = osp.join("/scratch/p/peltier/dchandan/ctest", model, "hist")
      comp_direc = ospath.join(case.DOUT_S_ROOT, args.model, "hist")
      # os.chdir(comp_direc)

      # Generating the list of files that need to be worked on
      list_of_files = []
      for year in years:
          pattern = "{0}/{1}.{2}.{3}.{4:04d}*.nc".format(comp_direc,
                                                        args.case, 
                                                        mtypes[args.model], 
                                                        freqtypes[args.model], 
                                                        year)
          list_of_files.extend(glob.glob(pattern))

      list_of_files.sort()

      print("+--------------------------------------+")
      print("| Number of files to operate upon: {0:3d} |".format(len(list_of_files)))
      print("+--------------------------------------+")

    else:
      list_of_files = None

    # We broadcast the list of files to all processes
    list_of_files = simplecomm._comm.bcast(list_of_files, root=0)


    output_prefix = "tseries_{0}_{1}.".format(args.start, args.end)

    # Main run function for python 2.7
    # Create the input object for the Reshaper
    spec = create_specifier(infiles=list_of_files,
                            ncfmt="netcdf4c",
                            deflate=args.deflate,
                            prefix=output_prefix,
                            suffix=".nc",
                            outdir=args.odir,
                            metadata=args.metadata)

    # Create the PyReshaper object
    reshpr = create_reshaper(spec,
                             serial=args.serial,
                             verbosity=args.verbosity,
                             skip_existing=args.skip_existing,
                             overwrite=args.overwrite,
                             once=args.once,
                             simplecomm=simplecomm,
                             backend=args.backend,
                             timecode=args.timecode,
                             sort_files=False,
                             preprocess=args.preprocess)

    # Run the conversion (slice-to-series) process
    reshpr.convert(output_limit=args.limit)

    # Print timing diagnostics
    reshpr.print_diagnostics()


#==============================================================================
# Command-line Opeartion
#==============================================================================
if __name__ == '__main__':
    args = cli()

    # No need to do checks here for positional arguments as this is 
    # automatically taken care of by argparse

    # Checking for output directories
    if not ospath.isdir(args.odir):
      raise ValueError("Invalid output directory {0}. ".format(args.odir))

    main(args)

