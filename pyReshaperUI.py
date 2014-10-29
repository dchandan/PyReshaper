"""
pyReshaperUI

This is the user-interface code for the PyNIO-MPI Reshaper.  It is a class
that reads and defines the user input for the code.
"""

import os
import xml.etree.ElementTree as ET


#===============================================================================
# User-Interface Class
#===============================================================================
class pyReshaperUI(object):
    """
    Class for reading/writing user-input XML files for the pyReshaper.  It
    reads an XML definition file, against which the specified input XML input
    is validated.
    """

    def __init__(self, def_file='pyreshaper_definition.xml'):
        """
        Constructor

        INPUT:
        ------

        def_file     The pyReshaper XML UI definition file, listing the user
                     options by name, type, default values, and validation
                     conditions.  By default, it looks for a file in the
                     current working directory.
        """

        # Check if the definition file exists
        abs_def_file = os.path.abspath(def_file)
        if (not os.path.exists(abs_def_file)):
            errmsg = 'Cannot find pyReshaper configuration definition file:\n'\
                   + str(abs_def_file)
            raise ValueError(errmsg)

        # Read the definition file and store in an ElementTree
        file_etree = ET.parse(abs_def_file)
        root = file_etree.getroot()
        if (root.tag != 'pyreshaper'):
            errmsg = 'Cannot find root of pyReshaper XML definition file'
            raise ValueError(errmsg)

        # Convert the ElementTree into a dictionary for internal use
        self.definition = root
        print ET.tostring(self.definition)

    def read_config(self, config_file):
        """
        Read a user-input XML configuration file.  This looks for a single
        "pyreshaper" node in the XML file, which is used as the root of the
        XML element-tree.  Multiple child "config" nodes below the root node
        can be specified.  If multiple "pyreshaper" nodes exist, it takes the
        first.

        INPUT:
        ------

        config_file  The pyReshaper XML input file, specifying how the
                     reshaper is to be run.  Must be specified.
        """
        # Check if the configuration file exists
        abs_config_file = os.path.abspath(config_file)
        if (not os.path.exists(abs_config_file)):
            errmsg = 'Cannot find pyReshaper XML configuration file:\n'\
                   + str(abs_config_file)
            raise ValueError(errmsg)

        # Read the configuration file and store in an ElementTree
        file_etree = ET.parse(abs_config_file)
        root = file_etree.getroot()
        if (root.tag != 'pyreshaper'):
            errmsg = 'Cannot find root of pyReshaper XML configuration file'
            raise ValueError(errmsg)

        # Set the configuration ElementTree to the pyReshaper root node
        self.config = root

        # Validate the configuration
        print ET.tostring(self.config, method=0)


#import argparse
#import fileinput
#import re
#import sys
#
#
#def parseInput():
#    """
#    The standard input parser.
#
#    This looks for command-line input, and if it is not specified,
#    then it assumes that namelist-formatted input will be supplied
#    on stdin.
#    """
#
#    # Get the command-line input
#    cliArgs = parseCommandLineInput()
#
#    # Get the arguments from a namelist (stdin)
#    nmlArgs = parseNamelistInput()
#
#    # Combine the two inputs (CLI overrides NML)
#    args = nmlArgs.copy()
#    for key in cliArgs:
#        args[key] = cliArgs[key]
#
#    # Validate the input
#    validateInput(args)
#
#    # If and ONLY if dictionary is incomplete
#    return args
#
#
#def parseCommandLineInput():
#    """
#    Use the ArgumentParser to parse command-line interface input.
#    """
#
#    parser = argparse.ArgumentParser(
#        description='''
#        Convert a set of time-slice (history) NetCDF files into time-series
#        NetCDF files, essentially a transpose operation on the data.
#        ''')
#
#    # Definethe CLI arguments (and create the help menu)
#    parser.add_argument("--file_type", dest="file_type", type=str,
#                            help="Type of netCDF to use (pNetCDF, netcdf, netcdf4)",
#                            metavar="FILE_TYPE")
#    parser.add_argument("--input_directory", dest="input_directory", type=str,
#                            help="Location of time-slice files to be read as input",
#                            metavar="INPUT_DIR")
#    parser.add_argument("--output_directory", dest="output_directory", type=str,
#                            help="Location of time-series files to be written as output",
#                            metavar="OUTPUT_DIR")
#    parser.add_argument("--start_time", dest="start_time", type=str,
#                            help="Beginning date stamp of time-slice file to be read",
#                            metavar="START_TIME")
#    parser.add_argument("--end_time", dest="end_time", type=str,
#                            help="Ending date stamp of time-slice file to be read",
#                            metavar="END_TIME")
#    parser.add_argument("--root_name", dest="root_name", type=str,
#                            help="Root name string of the time-slice files",
#                            metavar="ROOT_NAME")
#    parser.add_argument("--input_time_peried", dest="input_time_period", type=str,
#                            help="Time-slice periodicity in the input files",
#                            metavar="INPUT_TIME_PERIOD")
#    parser.add_argument("--output_time_option", dest="output_time_option", type=str,
#                            help="Time span units for a series file length",
#                            metavar="OUTPUT_TIME_OPTION")
#    parser.add_argument("--output_time_n", dest="output_time_n", type=int,
#                            help="Number of time span units in a series file",
#                            metavar="OUTPUT_TIME_N")
#    parser.add_argument("--meta_vars", dest="meta_vars", type=str,
#                            help="Non-time-series variables to be included in all time-series files",
#                            metavar="META_VARS")
#
#    parser.add_argument("--bounds_name", dest="bounds_names", type=str,
#                            help="Names for a list of bounds regions",
#                            metavar="BOUNDS_NAME", nargs='+')
#    parser.add_argument("--bounds_definition", dest="bounds_definition", type=str,
#                            help="Bounds definition for each bounds region",
#                            metavar="BOUNDS_DEFINITION", nargs='+')
#
#    # Parse the user arguments as a dictionary
#    args = vars(parser.parse_args())
#
#    # Remove unspecified values from the argument dictionary
#    for key in list(args.keys()):
#        if (args[key] == None):
#            del args[key]
#
#    # Return the argument dictionary
#    return args
#
#
#def parseNamelistInput():
#    """
#    Parse the namelist input from stdin
#    """
#
#    # Store standard input into a string
#    reading = False
#    args = {}
#    for line in fileinput.input('-'):
#
#        # Strip out comments
#        loc = line.find('!')
#        if (loc >= 0):
#            line = line[:loc]
#        line = line.strip()
#
#        # Search for start/end of namelist block
#        if (line[:10] == '&input_nml'):
#            reading = True
#        elif (line[0] == '/'):
#            reading = False
#
#        # Continue only if 'reading' inside the namelist
#        elif (reading):
#            loc = line.find('=')
#
#            # If there is an '=' character, then split into (object, value)
#            if (loc >= 0):
#                (obj, val) = line.split('=', 1)
#                obj = obj.strip()
#                val = val.strip()
#
#                # Remove quotes from string literals
#                val = val.replace("'", '').replace('"', '')
#
#                # Exchange ',' separators for whitespace
#                val = val.replace(',', ' ')
#
#                # Split on whitespace
#                valarray = val.split()
#
#                # De-listify if only 1 element present
#                if (len(valarray) == 1):
#                    val = valarray[0]
#                else:
#                    val = valarray
#
#                # Add the object and value to the dictionary
#                args[obj] = val
#
#    # Convert integer/float values, if possible
#    def _convertVal(obj, t):
#        if (obj in args):
#            try:
#                args[obj] = t(args[obj])
#            except:
#                pass  # Keep as string if fails
#
#    def _convertArr(obj, t):
#        if (obj in args):
#            try:
#                for i in range(len(args[obj])):
#                    args[obj][i] = t(args[obj][i])
#            except:
#                pass  # Keep as string if fails
#    _convertVal('output_time_n', int)
#    _convertArr('bounds_min', float)
#    _convertArr('bounds_max', float)
#
#    # Return the argument dictionary
#    return args
#
#
#def validateInput(args):
#    """
#    Check to make sure the input is complete and correctly formatted.
#    """
#
#    # String error message storage (do not print errors until end)
#    errorMessages = []
#
#    # Simple checking of key/value types
#    def _validateReqKey(key):
#        if (key not in args):
#            errorMessages.append('ERROR: Required input parameter "' + key + '" not specified.')
#            return False
#        return True
#
#    def _validateStrVal(key):
#        if (type(args[key]) is not str):
#            errorMessages.append('ERROR: Input parameter "' + key + '" must be a string.')
#            errorMessages.append('         ' + key + ' = ' + str(args[key]))
#            return False
#        return True
#
#    def _validateIntVal(key):
#        if (type(args[key]) is not int):
#            errorMessages.append('ERROR: Input parameter "' + key + '" must be an integer.')
#            errorMessages.append('         ' + key + ' = ' + str(args[key]))
#            return False
#        return True
#
#    def _validateFloatVal(key):
#        if (type(args[key]) is not float):
#            errorMessages.append('ERROR: Input parameter "' + key + '" must be a float.')
#            errorMessages.append('         ' + key + ' = ' + str(args[key]))
#            return False
#        return True
#
#    # Required value checking
#    def _validateReqStrVal(key):
#        if (_validateReqKey(key)):
#            return _validateStrVal(key)
#        return False
#
#    def _validateReqIntVal(key):
#        if (_validateReqKey(key)):
#            return _validateIntVal(key)
#        return False
#
#    def _validateReqFloatVal(key):
#        if (_validateReqKey(key)):
#            return _validateFloatVal(key)
#        return False
#
#    # Special value checking
#    def _validateReqStrValOpts(key, opts):
#        if (_validateReqStrVal(key)):
#            if (args[key] not in opts):
#                errorMessages.append('ERROR: Input parameter "' + key + '" must be one of:')
#                errorMessages.append('         ' + str(opts))
#                errorMessages.append('         ' + key + ' = ' + str(args[key]))
#                return False
#            return True
#        return False
#
#    def _validateReqDateVal(key):
#        if (_validateReqStrVal(key)):
#            lenArg = len(args[key])
#            if (lenArg == 4):
#                pattern = re.compile(r'\d{4}')
#            if (lenArg == 7):
#                pattern = re.compile(r'\d{4}-\d{2}')
#            if (lenArg == 10):
#                pattern = re.compile(r'\d{4}-\d{2}-\d{2}')
#            if (lenArg == 16):
#                pattern = re.compile(r'\d{4}-\d{2}-\d{2}-\d{5}')
#            if (pattern.match(args[key]) is None):
#                errorMessages.append('ERROR: Date-Time string "' + key + '" must have the format YYYY-MM-DD-SSSSS')
#                errorMessages.append('         ' + key + ' = ' + str(args[key]))
#                return False
#            return True
#        return False
#
#    # Optional argument checks
#    def _validateOptStrArr(key):
#        if (key in args):
#            valid = True
#            i = 0
#            while (valid == True and i < len(args[key])):
#                if (type(args[key][i]) is not str):
#                    valid = False
#                i = i + 1
#            if (not valid):
#                errorMessages.append('ERROR: Input parameter "' + key + '" must be an array of type str')
#                errorMessages.append('         ' + key + ' = ' + str(args[key]))
#                return False
#        return True
#
#    def _validateOptIntArr(key):
#        if (key in args):
#            valid = True
#            i = 0
#            while (valid == True and i < len(args[key])):
#                if (type(args[key][i]) is not int):
#                    valid = False
#                i = i + 1
#            if (not valid):
#                errorMessages.append('ERROR: Input parameter "' + key + '" must be an array of type int')
#                errorMessages.append('         ' + key + ' = ' + str(args[key]))
#                return False
#        return True
#
#    def _validateOptFloatArr(key):
#        if (key in args):
#            valid = True
#            i = 0
#            while (valid == True and i < len(args[key])):
#                if (type(args[key][i]) is not float):
#                    valid = False
#
#                i = i + 1
#            if (not valid):
#                errorMessages.append('ERROR: Input parameter "' + key + '" must be an array of type float')
#                errorMessages.append('         ' + key + ' = ' + str(args[key]))
#                return False
#        return True
#
#    # Validate required arguments
#    _validateReqStrValOpts('file_type', ['netcdf', 'netcdf4'])
#    _validateReqStrVal('input_directory')
#    _validateReqStrVal('output_directory')
#    _validateReqDateVal('start_time')
#    _validateReqDateVal('end_time')
#    _validateReqStrVal('root_name')
#    _validateReqStrValOpts('input_time_period', ['mon', 'daily', 'hourly'])
#    _validateReqStrValOpts('output_time_option', ['nyears', 'nmonths', 'ndays', 'nhours'])
#    _validateReqIntVal('output_time_n')
#
#    # Validate optional arguments
#    _validateOptStrArr('meta_vars')
#    _validateOptStrArr('bounds_name')
#    _validateOptFloatArr('bounds_min')
#    _validateOptFloatArr('bounds_max')
#
#    # Check the concurrence (all or none) of the optional bounds arrays
#    namePresent = 'bounds_name' in args
#    minPresent = 'bounds_min' in args
#    maxPresent = 'bounds_max' in args
#    allPresent = namePresent and minPresent and maxPresent
#    nonePresent = (not namePresent) and (not minPresent) and (not maxPresent)
#    if (not (allPresent or nonePresent)):
#        errorMessages.append('ERROR: Input parameters "bounds_name", "bounds_min", and "bounds_max" must all be present if used.')
#
#    # Check length of optional arrays
#    if (allPresent):
#        lenName = len(args['bounds_name'])
#        lenMin = len(args['bounds_min'])
#        lenMax = len(args['bounds_max'])
#        allSame = (2 * lenName == lenMin) and (lenMin == lenMax)
#        if (not allSame):
#            errorMessages.append('ERROR: Input parameters "bounds_min" and "bounds_max" must be twice the length of "bounds_name"')
#
#    # Print error messages and exit, if any
#    if (len(errorMessages) > 0):
#        print "There were errors in the input supplied.  See below:"
#        print
#        for msg in errorMessages:
#            print msg
#        sys.exit(1)
