#!/usr/bin/env python
"""
This script provides the command-line interface (CLI) to the PyReshaper

This script uses the optparse module for Python 2.6 and the newer argparse
module for Python 2.7.

Copyright 2015, University Corporation for Atmospheric Research
See the LICENSE.txt file for details
"""

import os, glob
import os.path as ospath
from pyreshaper.specification import create_specifier
from pyreshaper.reshaper import create_reshaper

import argparse

from PyCESM.case.CESMCase import CESMCase


def cli():
    parser = argparse.ArgumentParser(prog="{0}: Convert CESM slice files to series files".format(__file__))
    parser.add_argument('case',  type=str, help='CESM case name')
    parser.add_argument('start', type=int, help='start year')
    parser.add_argument('end',   type=int, help='end year')
    parser.add_argument('model', type=str, choices=["atm", "lnd", "ocn", "ice"], help='models component to process')
    parser.add_argument("odir",  type=str, action='store', help="Directory to write output files into")

    # ncflags = parser.add_argument_group('nccopy flags')
    # ncflags.add_argument('-d',   type=int, default=6, help="compression level")


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
    print("Configuration  >>>>>>>")
    print("  Case name  : {0}".format(args.case))
    print("  Model      : {0}".format(args.model))
    print("  Start year : {0}".format(args.start))
    print("  End year   : {0}".format(args.end))
    print("  Compression: {0}".format(args.deflate))
    print("  Verbosity  : {0}".format(args.verbosity))
    print("<<<<<<<<<<<<<<<<<<<<<<")

    # The CESM case object
    case        = CESMCase(casename)

    total_years = args.end - args.start + 1
    years       = range(args.start, args.end + 1)
    
    # A mapping between model names and file corresponding file names
    mtypes = {"atm":"cam2", "lnd":"clm2", "ocn":"pop", "ice":"cice"}


    #comp_direc = osp.join("/scratch/p/peltier/dchandan/ctest", model, "hist")
    comp_direc = ospath.join(case.DOUT_S_ROOT, args.model, "hist")
    os.chdir(comp_direc)

    # Generating the list of files that need to be worked on
    list_of_files = []
    for year in years:
        pattern = "{0}.{1}.h*.{2:04d}*.nc".format(args.case, mtypes[args.model], year)
        list_of_files.extend(glob.glob(pattern))

    print("+--------------------------------------+")
    print("| Number of files to operate upon: {0:3d} |".format(len(list_of_files)))
    print("+--------------------------------------+")

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
                             once=args.once)

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

