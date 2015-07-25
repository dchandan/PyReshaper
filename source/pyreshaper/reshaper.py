"""
The module containing the Reshaper class

This is the main reshaper module.  This is where the specific operations
are defined.  Currently, only one operation has been implemented (i.e.,
the time-slice to time-series operation).

Copyright 2015, University Corporation for Atmospheric Research
See the LICENSE.txt file for details
"""

# Built-in imports
import abc
import os
import itertools

# Third-party imports
import Nio
import netCDF4
import numpy
from asaptools.simplecomm import create_comm, SimpleComm, SimpleCommMPI
from asaptools.timekeeper import TimeKeeper
from asaptools.partition import WeightBalanced
from asaptools.vprinter import VPrinter

# PyReshaper imports
from specification import Specifier


#==============================================================================
# create_reshaper factory function
#==============================================================================
def create_reshaper(specifier, serial=False, verbosity=1,
                    skip_existing=False, overwrite=False,
                    once=False, simplecomm=None, backend="netcdf",
                    timecode=True, preprocess=True, sort_files=True):
    """
    Factory function for Reshaper class instantiations.

    Parameters:
        specifier (Specifier): An instantiation of a Specifier class that 
            defines the type of operation to be performed.  (That is, a
            Slice2SeriesSpecifier will invoke the creation of a
            matching Slice2SeriesReshaper object.)
            Alternatively, this can be a list of Specifier objects,
            or dictionary of named Specifier objects.
            In this case, a reshaper will be created for each
            specifier in the list, and each reshaper will be
            created and run in sequence.

    Keyword Arguments:
        serial (bool): True or False, indicating whether the Reshaper object
            should perform its operation in serial (True) or
            parallel (False).
        verbosity (int): Level of printed output (stdout).  A value of 0 means
            no output, and a higher value means more output.  The
            default value is 1.
        skip_existing (bool): Flag specifying whether to skip the generation
            of time-series for variables with time-series files that already
            exist.  Default is False.
        overwrite (bool): Flag specifying whether to forcefully overwrite
            output files if they already exist.  Default is False.
        once (bool): True or False, indicating whether the Reshaper should
            write all metadata to a 'once' file (separately).
        simplecomm (SimpleComm): A SimpleComm object to handle the parallel
            communication, if necessary

    Returns:
        Reshaper: An instance of the Reshaper object requested
    """
    # Determine the type of Reshaper object to instantiate
    if type(specifier) is Specifier:
        return Slice2SeriesReshaper(specifier,
                                    serial=serial,
                                    verbosity=verbosity,
                                    skip_existing=skip_existing,
                                    overwrite=overwrite,
                                    once=once,
                                    simplecomm=simplecomm,
                                    backend=backend,
                                    timecode=timecode,
                                    preprocess=preprocess,
                                    sort_files=sort_files)
    elif isinstance(specifier, list):
        spec_dict = dict([(str(i), s) for (i, s) in enumerate(specifier)])
        return create_reshaper(spec_dict,
                               serial=serial,
                               verbosity=verbosity,
                               skip_existing=skip_existing,
                               overwrite=overwrite,
                               once=once,
                               simplecomm=simplecomm,
                               backend=backend,
                               timecode=timecode,
                               preprocess=preprocess,
                               sort_files=sort_files)
    elif isinstance(specifier, dict):
        spec_types = set([type(s) for s in specifier.values()])
        if len(spec_types) > 1:
            err_msg = 'Multiple specifiers must all have the same type'
            raise TypeError(err_msg)
        spec_type = spec_types.pop()
        if spec_type is Specifier:
            return MultiSpecReshaper(specifier,
                                     serial=serial,
                                     verbosity=verbosity,
                                     skip_existing=skip_existing,
                                     overwrite=overwrite,
                                     once=once,
                                     simplecomm=simplecomm)
        else:
            err_msg = 'Multiple specifiers of type {0} are not valid'.format(spec_type)
            raise TypeError(err_msg)
    else:
        err_msg = 'Specifier type {0} not a valid Specifier object'.format(type(specifier))
        raise TypeError(err_msg)


#==============================================================================
# _pprint_dictionary - Helper method for printing diagnostic data
#==============================================================================
def _pprint_dictionary(title, dictionary, order=None):
    """
    Hidden method for pretty-printing a dictionary of numeric values,
    with a given title.

    Parameters:
        title (str): The title to give to the printed table
        dictionary (dict): A dictionary of numeric values

    Keyword Arguments:
        order (list): The print order for the keys in the dictionary (only 
            items that are in both the order list and the dictionary will be
            printed)

    Return:
        str: A string with the pretty-printed dictionary data
    """
    # Type checking
    if (type(title) is not str):
        err_msg = 'Title must be a str type'
        raise TypeError(err_msg)
    if (not isinstance(dictionary, dict)):
        err_msg = 'Input dictionary needs to be a dictionary type'
        raise TypeError(err_msg)
    if (order != None and not isinstance(order, list)):
        err_msg = 'Order list needs to be a list type'
        raise TypeError(err_msg)

    # Determine the print order, if present
    print_order = dictionary.keys()
    if (order != None):
        print_order = []
        for item in order:
            if (item in dictionary):
                print_order.append(item)

    # Header line with Title
    hline = '-' * 50 + os.linesep
    ostr = hline + ' ' + title.upper() + ':' + os.linesep + hline

    # Determine the longest timer name
    # and thus computer the column to line up values
    valcol = 0
    for name in print_order:
        if (len(str(name)) > valcol):
            valcol = len(str(name))
    valcol += 2

    # Print out the timer names and the accumulated times
    for name in print_order:
        spacer = ' ' * (valcol - len(str(name)))
        ostr += str(name) + ':' + spacer
        ostr += str(dictionary[name]) + os.linesep
    ostr += hline

    return ostr


#==============================================================================
# Reshaper Abstract Base Class
#==============================================================================
class Reshaper(object):

    """
    Abstract base class for Reshaper objects
    """

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def convert(self):
        """
        Method to perform the Reshaper's designated operation.
        """
        return

    @abc.abstractmethod
    def print_diagnostics(self):
        """
        Print out timing and I/O information collected up to this point
        """
        return


#==============================================================================
# Reshaper Class
#==============================================================================
class Slice2SeriesReshaper(Reshaper):

    """
    The time-slice to time-series Reshaper class

    This is the class that defines how the time-slice to time-series 
    reshaping operation is to be performed.
    """

    def __init__(self, specifier, serial=False, verbosity=1,
                 skip_existing=False, overwrite=False,
                 once=False, simplecomm=None, backend="netcdf",
                 timecode=True, preprocess=True, sort_files=True):
        """
        Constructor

        Parameters:
            specifier (Specifier): An instance of the Specifier class, 
                defining the input specification for this reshaper operation.

        Keyword Arguments:
            serial (bool): True or False, indicating whether the operation
                should be performed in serial (True) or parallel
                (False).  The default is to assume parallel operation
                (but serial will be chosen if the mpi4py cannot be
                found when trying to initialize decomposition.
            verbosity(int): Level of printed output (stdout).  A value of 0 
                means no output, and a higher value means more output.  The
                default value is 1.
            skip_existing (bool): Flag specifying whether to skip the generation
                of time-series for variables with time-series files that already
                exist.  Default is False.
            overwrite (bool): Flag specifying whether to forcefully overwrite
                output files if they already exist.  Default is False.
            once (bool): True or False, indicating whether the Reshaper should
                write all metadata to a 'once' file (separately).
            simplecomm (SimpleComm): A SimpleComm object to handle the parallel 
                communication, if necessary
        """

        # Type checking (or double-checking)
        if not isinstance(specifier, Specifier):
            err_msg = "Input must be given in the form of a Specifier object"
            raise TypeError(err_msg)
        if type(serial) is not bool:
            err_msg = "Serial indicator must be True or False."
            raise TypeError(err_msg)
        if type(verbosity) is not int:
            err_msg = "Verbosity level must be an integer."
            raise TypeError(err_msg)
        if type(skip_existing) is not bool:
            err_msg = "Skip_existing flag must be True or False."
            raise TypeError(err_msg)
        if type(once) is not bool:
            err_msg = "Once-file indicator must be True or False."
            raise TypeError(err_msg)
        if simplecomm is not None:
            if not (isinstance(simplecomm, SimpleComm) or \
                    isinstance(simplecomm, SimpleCommMPI)):
                err_msg = ("Simple communicator object is not a SimpleComm or ",
                           "SimpleCommMPI")
                raise TypeError(err_msg)

        # Whether to write a once file
        self._use_once_file = once
        self._preprocess = preprocess
        self._timecode   = timecode
        self._need_to_sort_files = sort_files
        self._backend = backend
        self._validate_backend()


        # We need the self.preprocess to be True if we need to sort files, as
        # sorting files is a part of preprocessing.
        if self._need_to_sort_files: self._preprocess = True

        # Internal timer data
        if self._timecode: self._timer = TimeKeeper()

        # Dictionary storing read/write data amounts
        self.assumed_block_size = float(4 * 1024 * 1024)
        self._byte_counts = {}

        if self._timecode: self._timer.start('Initializing Simple Communicator')
        if simplecomm is None:
            simplecomm = create_comm(serial=serial)
        # Reference to the simple communicator
        self._simplecomm = simplecomm
        if self._timecode: self._timer.stop('Initializing Simple Communicator')

        # Contruct the print header
        header = ''.join(['[', str(self._simplecomm.get_rank()),
                          '/', str(self._simplecomm.get_size()), '] '])

        # Reference to the verbose printer tool
        self._vprint = VPrinter(header=header, verbosity=verbosity)

        # Debug output starting
        if self._simplecomm.is_manager():
            self._vprint('Initializing Reshaper', verbosity=1)

        # Validate the user input data
        if self._timecode: self._timer.start('Specifier Validation')
        specifier.validate()
        if self._timecode: self._timer.stop('Specifier Validation')
        if self._simplecomm.is_manager():
            self._vprint('Specifier validated', verbosity=1)


        if (self._backend == "nio"):
            # Setup PyNIO options (including disabling the default PreFill option)
            opt = Nio.options()
            opt.PreFill = False

            # Determine the Format and CompressionLevel options
            # from the NetCDF format string in the Specifier
            if specifier.netcdf_format == 'netcdf':
                opt.Format = 'Classic'
            elif specifier.netcdf_format == 'netcdf4':
                opt.Format = 'NetCDF4Classic'
                opt.CompressionLevel = 0
            elif specifier.netcdf_format == 'netcdf4c':
                opt.Format = 'NetCDF4Classic'
                opt.CompressionLevel = specifier.netcdf_deflate
                if self._simplecomm.is_manager():
                    self._vprint('PyNIO compression level: {0}'.format(\
                        specifier.netcdf_deflate), verbosity=2)

            self._nio_options = opt
            if self._simplecomm.is_manager():
                self._vprint('PyNIO options set', verbosity=2)
        else:
            self._netcdf_dataset_options = {}
            self._netcdf_dim_options = {}
            self._netcdf_var_options = {}
            if specifier.netcdf_format == 'netcdf':
                self._netcdf_dataset_options["format"]  = "NETCDF3_64BIT"
            elif specifier.netcdf_format == 'netcdf4':
                self._netcdf_dataset_options["format"]  = "NETCDF4_CLASSIC"
            elif specifier.netcdf_format == 'netcdf4c':
                self._netcdf_dataset_options["format"]  = "NETCDF4"
                self._netcdf_var_options["zlib"] = True
                self._netcdf_var_options["complevel"] = specifier.netcdf_deflate
                if self._simplecomm.is_manager():
                    self._vprint('netCDF4 compression level: {0}'.format(\
                        specifier.netcdf_deflate), verbosity=2)

            if self._simplecomm.is_manager():
                self._vprint('netCDF4 options set', verbosity=2)

       
        # This is an abstracted attribute to open input files transparently
        # regardless of the backend
        if (self._backend == "nio"):
            self.Open_Input_File = Nio.open_file
        else:
            self.Open_Input_File = netCDF4.Dataset


        # Copying the list of input file names into this class
        self.input_file_list = specifier.input_file_list
        self.num_input_files = len(self.input_file_list)


        # Analyze the input files for validation and getting information. We do
        # this only if the preprocess flag is True
        if self._timecode: self._timer.start('Input File Validation')
        self._analyze_input_files()
        if self._timecode: self._timer.stop('Input File Validation')
        if self._simplecomm.is_manager():
            self._vprint('Input files validated', verbosity=2)


        # Retrieve and sort the variables in each time-slice file
        # (To determine if it is time-invariant metadata, time-variant
        # metadata, or if it is a time-series variable)
        if self._timecode: self._timer.start('Sort Variables')
        self._sort_variables(specifier)
        if self._timecode: self._timer.stop('Sort Variables')
        if self._simplecomm.is_manager():
            self._vprint('Variables sorted', verbosity=2)

        # Validate the output files
        if self._timecode: self._timer.start('Output File Validation')
        self._validate_output_files(specifier, skip_existing, overwrite)
        if self._timecode: self._timer.stop('Output File Validation')
        if self._simplecomm.is_manager():
            self._vprint('Output files validated', verbosity=2)

        # Helpful debugging message
        if self._simplecomm.is_manager():
            self._vprint('Reshaper initialized.', verbosity=1)

        # Sync before continuing..
        self._simplecomm.sync()


    def _get_unlim_dim_len(self, fhandle, unlim_dim_name):
        """
        Backend abstraction method to get the length of the unlimited dimension.

        ARGUMENTS
            fhandle - opened netcdf file handle returned from either backend library
            unlim_dim_name - name of the unlimited dimension
        RETURNS
            Length og the unlimited dimension
        """
        return fhandle.dimensions[unlim_dim_name] if (self._backend == "nio") else \
               len(fhandle.dimensions[unlim_dim_name])


    def _get_var_ndims(self, var):
        """
        Backend abstraction method to get the number of dimensions for a variable

        ARGUMENTS
            var - variable object returned from either backend library
        RETURNS
            Number of dimensions
        """
        return len(var.dimensions) if (self._backend == "nio") else var.ndim


    def _create_nc_variable(self, out_file, out_name, in_var):
        """
        Backend abstraction method to create a netcdf variable in an output file.

        ARGUMENTS
            out_file - output file handle from either backend
            out_name - name of the variable to create
            in_var   - variable object return from either backend and which corresponds
                       to a variable in an input file.
        RETURNS
            output variable object
        """
        if (self._backend == "nio"):
            return out_file.create_variable(out_name,
                                            in_var.typecode(), 
                                            in_var.dimensions)
        else:
            return out_file.createVariable(out_name,
                                           in_var.dtype,
                                           in_var.dimensions,
                                           **self._netcdf_var_options)


    def _copy_nc_var_attributes(self, out_var, in_var):
        """
        Copy netcdf attributes from one file to amother
        ARGUMENTS
            out_var - variable object from either backend
            in_var  - variable object from either backend
        RETURNS
            variable object for the output variable
        """
        if (self._backend == "nio"):
            for att_name, att_val in in_var.attributes.iteritems():
                setattr(out_var, att_name, att_val)
        else:
            out_var.setncatts(in_var.__dict__)
        return out_var



    def _validate_backend(self):
        if (self._backend not in ["nio", "netcdf"]):
            msg = "Netcdf backend library must be either 'nio' or 'netcdf'"
            raise ValueError(msg)


    def _analyze_input_files(self):
        """
        Perform validation of input data files themselves and analyzes them for
        information that is helpful in the pre-process stage.  
        """


        # The analysis and extraction of relevant information only needs to be
        # performed by one process, since all processes operate on the same list
        # of input files. The information that needs to be shared among all processes
        # in then "broadcasted" to all using MPI's broadcast communication.
        if self._simplecomm.is_manager():
            # Helpful debugging message
            self._vprint('Validating input files', verbosity=1)

            # In the first file, look for the 'unlimited' dimension and get its variables
            ifile = self.Open_Input_File(self.input_file_list[0], "r")
            self._unlimited_dim = None
            for dim in ifile.dimensions:
                if (self._backend == "nio"): 
                    _isunlim = ifile.unlimited(dim)
                else:
                    _isunlim = ifile.dimensions[dim].isunlimited()
                if _isunlim:
                    self._unlimited_dim = dim
                    break  # There can only be 1!
            if self._unlimited_dim == None:
                err_msg = 'Unlimited dimension not identified.'
                raise LookupError(err_msg)

            variables = ifile.variables
            var_names = set(variables.keys())

            ifile.close()

            missing_vars   = set()
            time_values    = []
            broadcast_data = {}


            # Make a pass through each file and:
            # (1) Make sure it has the 'unlimited' dimension
            # (2) Make sure this dimension is truely 'unlimited'
            # (3) Check that this dimension has a corresponding variable
            # (4) Make sure that the list of variables in each file is the same
            # (5) Get the time values for each file so that we can sort the files
            if self._preprocess:  # We only do this if we've been asked to preprocess
                for i in range(self.num_input_files):
                    ifile = self.Open_Input_File(self.input_file_list[i], "r")

                    # Make sure it has the 'unlimited' dimension
                    if self._unlimited_dim not in ifile.dimensions:
                        err_msg = 'Unlimited dimension not found in file ({0})'.\
                                  format(self.input_file_list[i])
                        raise LookupError(err_msg)

                    # Make sure this dimension is truely 'unlimited'
                    if (self._backend == "nio"):
                        if not ifile.unlimited(self._unlimited_dim):
                            err_msg = 'Unlimited dimension not unlimited in file ({0})'.\
                                      format(self.input_file_list[i])
                            raise LookupError(err_msg)
                    else:
                        pass
                    
                    # Check that this dimension has a corresponding variable
                    if self._unlimited_dim not in ifile.variables:
                        err_msg = 'Unlimited dimension variable not found in file ({0})'.\
                                  format(self.input_file_list[i])
                        raise LookupError(err_msg)
                    
                    # Make sure that the list of variables in each file is the same            
                    var_names_next = set(ifile.variables.keys())
                    missing_vars.update(var_names - var_names_next)

                    # Get the time values for each file so that we can sort the files
                    # Only needed if the control flag need_to_sort_files is True
                    if self._need_to_sort_files:
                        if (self._backend == "nio"):
                            time_values.append(ifile.variables[self._unlimited_dim].get_value())
                        else:
                            time_values.append(ifile.variables[self._unlimited_dim][:])

                    ifile.close()

                if len(missing_vars) != 0:
                    warning = "WARNING: The first input file has variables that are " \
                        + "not in all input files:" + os.linesep + '   '
                    warning += " ".join(missing_vars)
                    self._vprint(warning, header=True, verbosity=1)


            broadcast_data['unlimdim'] = self._unlimited_dim
            
            # Now we sort the names of the input files
            # Only needed if the control flag need_to_sort_files is True
            if self._need_to_sort_files:
                new_order = sorted(range(self.num_input_files), 
                                   key=lambda i: time_values[i][0])
                broadcast_data['new_order'] = new_order
                broadcast_data['time_values'] = time_values

        else:
            broadcast_data = None


        broadcast_data = self._simplecomm._comm.bcast(broadcast_data, root=0)
        self._unlimited_dim = broadcast_data['unlimdim']
        
        # We sort the files only if specified by the control flag need_to_sort_files
        if self._need_to_sort_files:
            new_order       = broadcast_data['new_order']
            time_values     = broadcast_data['time_values']

            order = range(self.num_input_files)
            new_filenames = [None] * len(new_order)
            new_values    = [None] * len(new_order)
            for i in order:
                new_filenames[i] = self.input_file_list[new_order[i]]
                new_values[i]    = time_values[new_order[i]]

            # Save this data in the new orders
            self.input_file_list = new_filenames

            # Now, check that the largest time in each file is less than the
            # smallest time in the next file (so that the time spans of each file
            # do not overlap)
            for i in order[:-1]:
                if new_values[i][-1] >= new_values[i + 1][0]:
                    err_msg = 'Times in input files {0} and {1} appear to overlap'
                    err_msg = err_msg.format(new_filenames[i], new_filenames[i+1])
                    raise ValueError(err_msg)

        # Now that this is validated, let's string together the numpy array
        # of all times (using the new_values array)
        # self._all_time_values = \
        #     numpy.fromiter(itertools.chain.from_iterable(new_values),
        #                    dtype='float')




    def _sort_variables(self, specifier):
        """
        Internal method for sorting the variables in each time-slice file

        This method determines if each variable is to be treated as 
        time-invariant metadata, time-variant metadata (user defined), or 
        time-series variables.  All metadata is written to every time-series 
        file, and any time-series variable is written to its own file.  
        The time-variant metadata variables are determined by user input, 
        and are contained in the Specifier data member:

            Specifier.time_variant_metadata.

        Parameters:
            specifier (Specifier): The reshaper specifier object
        """

        # Helpful debugging message
        if self._simplecomm.is_manager():
            self._vprint('Sorting variables', verbosity=1)

        # Initialize the dictionary of variable names for each category
        # (Keys are variable names, Values are variable sizes)
        self._time_variant_metadata = {}
        self._time_invariant_metadata = {}
        self._time_series_variables = {}

        # Categorize each variable (only looking at first file)
        ifile = self.Open_Input_File(self.input_file_list[0], "r")
        variables = ifile.variables
        for var_name in variables.keys():
            var = variables[var_name]
            if (self._backend == "nio"):
                size = numpy.dtype(var.typecode()).itemsize
            else:
                size = var.dtype.itemsize
            size = size * numpy.prod(var.shape)
            if self._unlimited_dim not in var.dimensions:
                self._time_invariant_metadata[var_name] = size
            elif var_name in specifier.time_variant_metadata:
                self._time_variant_metadata[var_name] = size
            else:
                self._time_series_variables[var_name] = size

        # Debug output
        if self._simplecomm.is_manager():
            self._vprint('Time-Invariant Metadata: ' +
                         str(self._time_invariant_metadata.keys()), verbosity=2)
            self._vprint('Time-Variant Metadata: ' +
                         str(self._time_variant_metadata.keys()), verbosity=2)
            self._vprint('Time-Series Variables: ' +
                         str(self._time_series_variables.keys()), verbosity=2)

        # Add 'once' variable if writing to a once file
        # NOTE: This is a "cheat"!  There is no 'once' variable.  It's just
        #       a catch for all metadata IFF the 'once-file' is enabled.
        if self._use_once_file:
            self._time_series_variables['once'] = 1

        ifile.close()



    def _validate_output_files(self, specifier,
                               skip_existing=False, overwrite=False):
        """
        Perform validation of output data files themselves.  

        We compute the output file name from the prefix and suffix, and then
        we check whether the output files exist.  By default, if the output
        file

        Parameters:
            specifier (Specifier): The reshaper specifier object

        Keyword Arguments:
            skip_existing (bool): Flag specifying whether to skip the generation
                of time-series for variables with time-series files that already
                exist.  Default is False.
            overwrite (bool): Flag specifying whether to forcefully overwrite
                output files if they already exist.  Default is False.
        """

        # Helpful debugging message
        if self._simplecomm.is_manager():
            self._vprint('Validating output files', verbosity=1)

        # Loop through the time-series variables and generate output filenames
        prefix = specifier.output_file_prefix
        suffix = specifier.output_file_suffix
        self._time_series_filenames = \
            dict([(variable, prefix + variable + suffix)
                  for variable in self._time_series_variables])


        # This only needs to be performed by one process
        if self._simplecomm.is_manager():
            # Find which files already exist
            existing = []
            for variable, filename in self._time_series_filenames.items():
                if os.path.isfile(filename):
                    existing.append(variable)

            # If overwrite is enabled, delete all existing files first
            if overwrite:
                self._vprint('WARNING: Deleting existing output files for '
                             'time-series variables: {0}'.format(existing),
                             verbosity=1)
                for variable in existing:
                    os.remove(self._time_series_filenames[variable])

            # Or, if skip_existing is set, remove the existing time-series
            # variables from the list of time-series variables to convert
            elif skip_existing:
                self._vprint('WARNING: Skipping time-series variables with '
                             'existing output files: {0}'.format(existing),
                             verbosity=1)
                for variable in existing:
                    self._time_series_variables.pop(variable)

            # Otherwise, throw an exception if any existing output files are found
            elif len(existing) > 0:
                err_msg = ("Found existing output files for time-series "
                           "variables: {0}").format(existing)
                raise RuntimeError(err_msg)



    def convert(self, output_limit=0):
        """
        Method to perform the Reshaper's designated operation.

        In this case, convert a list of time-slice files to time-series files.

        Keyword Arguments:
            output_limit (int): Limit on the number of output (time-series) 
                files to write during the convert() operation.  If set
                to 0, no limit is placed.  This limits the number
                of output files produced by each processor in a
                parallel run.
        """
        # Type checking input
        if type(output_limit) is not int:
            err_msg = 'Output limit must be an integer'
            raise TypeError(err_msg)

        # Start the total convert process timer
        self._simplecomm.sync()
        if self._timecode: self._timer.start('Complete Conversion Process')

        # Debugging output
        if self._simplecomm.is_manager():
            self._vprint('Converting time-slices to time-series', verbosity=1)

        # For data common to all input files, we reference only the first
        ref_infile = self.Open_Input_File(self.input_file_list[0], "r")

        # Store the common dimensions and attributes for each file
        # (taken from the first input file in the list)
        common_dims = ref_infile.dimensions
        if (self._backend == "nio"):
            common_atts = ref_infile.attributes
        else:
            common_atts = ref_infile.__dict__

        # Partition the time-series variables across all processors
        tsv_names_loc = self._simplecomm.partition(self._time_series_variables.items(),
                                                   func=WeightBalanced(),
                                                   involved=True)
        if output_limit > 0:
            tsv_names_loc = tsv_names_loc[0:output_limit]

        # Print partitions for all ranks
        dbg_msg = 'Local time-series variables are {0}'.format(tsv_names_loc)
        self._vprint(dbg_msg, header=True, verbosity=2)

        # Reset all of the timer values (as it is possible that there are no
        # time-series variables in the local list procuded above)
        if self._timecode: 
            self._timer.reset('Open Output Files')
            self._timer.reset('Create Time-Invariant Metadata')
            self._timer.reset('Create Time-Variant Metadata')
            self._timer.reset('Create Time-Series Variables')
            self._timer.reset('Write Time-Invariant Metadata')
            self._timer.reset('Write Time-Variant Metadata')
            self._timer.reset('Write Time-Series Variables')
            self._timer.reset('Close Output Files')

        # Initialize the byte count dictionary
        self._byte_counts['Requested Data'] = 0
        self._byte_counts['Actual Data'] = 0

        # Defining a simple helper function to determine whether to
        # write time-series data and/or write metadata.  This is useful
        # for adding the ability to write a "once" file
        def _get_once_info(vname):
            is_once_file = (vname == 'once')
            write_meta = True
            write_tser = True
            if self._use_once_file:
                write_meta = is_once_file
                write_tser = not is_once_file
            return is_once_file, write_meta, write_tser

        # NOTE: In the prototype, we check for the existance of the output
        # directory at this point.  If it does not exist, we create it (but
        # only from the master rank).  This requires synchronization with
        # the decomp utility.  Instead, we assume the output directory
        # already exists (and is checked by the Specifier's validation).  No
        # synchronization is needed.

        # For each time-series variable, create the corresponding output file
        # (Also defines the header info for each output file)
        out_files = {}
        out_tvm_vars = {}
        for out_name in tsv_names_loc:
            is_once_file, write_meta, write_tser = _get_once_info(out_name)

            # Determine the output file name for this variable
            out_filename = self._time_series_filenames[out_name]
            dbg_msg = 'Creating output file for variable: {0}'.format(out_name)
            if is_once_file: dbg_msg = 'Creating "once" file.'
            self._vprint(dbg_msg, header=True, verbosity=1)

            # Open each output file and create the dimensions and attributes
            # NOTE: If the output file already exists, abort!
            if self._timecode: self._timer.start('Open Output Files')
            if os.path.exists(out_filename):
                err_msg = 'Found existing output file: {0}'.format(out_filename)
                raise OSError(err_msg)

            if (self._backend == "nio"):
                out_file = Nio.open_file(out_filename, 'w',
                                         options=self._nio_options)
                for att_name, att_val in common_atts.iteritems():
                    setattr(out_file, att_name, att_val)
                for dim_name, dim_val in common_dims.iteritems():
                    if dim_name == self._unlimited_dim:
                        out_file.create_dimension(dim_name, None)
                    else:
                        out_file.create_dimension(dim_name, dim_val)
            else:
                out_file = netCDF4.Dataset(out_filename, 'w', **self._netcdf_dataset_options)
                # Setting the common attributes
                out_file.setncatts(common_atts)
                # Creating dimensions
                for dim_name, dim_val in common_dims.iteritems():
                    if dim_name == self._unlimited_dim:
                        out_file.createDimension(dim_name, None)
                    else:
                        out_file.createDimension(dim_name, len(dim_val))

            
            if self._timecode: self._timer.stop('Open Output Files')

            # Create the time-invariant metadata variables
            if (write_meta):
                if self._timecode: self._timer.start('Create Time-Invariant Metadata')
                for name in self._time_invariant_metadata:
                    in_var = ref_infile.variables[name]
                    out_var = self._create_nc_variable(out_file, name, in_var)
                    out_var = self._copy_nc_var_attributes(out_var, in_var)
                if self._timecode: self._timer.stop('Create Time-Invariant Metadata')

            # Create the time-variant metadata variables
            if write_meta:
                if self._timecode: self._timer.start('Create Time-Variant Metadata')
                for name in self._time_variant_metadata:
                    in_var  = ref_infile.variables[name]
                    out_tvm_vars[name] = self._create_nc_variable(out_file, name, in_var)
                    out_tvm_vars[name] = self._copy_nc_var_attributes(out_tvm_vars[name], in_var)
                if self._timecode: self._timer.stop('Create Time-Variant Metadata')

            # Create the time-series variable itself
            if write_tser:
                if self._timecode: self._timer.start('Create Time-Series Variables')
                in_var  = ref_infile.variables[out_name]
                out_var = self._create_nc_variable(out_file, out_name, in_var)
                if self._timecode: self._timer.stop('Create Time-Series Variables')

            # Append the output file to list
            out_files[out_name] = out_file

        # Now that each output file has been created, start writing the data
        # (Looping over output file index, which is common in name lists)
        

        # First, we loop over all output files to create assign the variable
        # attributes and to write time-invariant metadata.
        for out_name, out_file in out_files.iteritems():
            is_once_file, write_meta, write_tser = _get_once_info(out_name)

            # Create the attributes of the time-series variable
            if write_tser:
                in_var = ref_infile.variables[out_name]
                out_var = out_file.variables[out_name]
                out_var = self._copy_nc_var_attributes(out_var, in_var)

            # Write the time-invariant metadata
            if write_meta:
                if self._timecode: self._timer.start('Write Time-Invariant Metadata')
                for name in self._time_invariant_metadata:
                    in_meta = ref_infile.variables[name]
                    out_meta = out_file.variables[name]
                    if (self._backend == "nio"):
                        if in_meta.rank > 0:
                            out_meta[:] = in_meta[:]
                        else:
                            out_meta.assign_value(in_meta.get_value())
                    else:
                        out_meta[:] = in_meta[:]

                if self._timecode: self._timer.stop('Write Time-Invariant Metadata')


        series_step_index = 0
        
        
        for j in xrange(self.num_input_files): # Loop over input files
            in_file = self.Open_Input_File(self.input_file_list[j], "r")
            # Get the number of time steps in this slice file
            num_steps = self._get_unlim_dim_len(in_file, self._unlimited_dim)

            # Loop over output variables
            for out_name, out_file in out_files.iteritems():
                is_once_file, write_meta, write_tser = _get_once_info(out_name)
                
                # Loop over the time steps in this slice file
                for slice_step_index in range(num_steps):

                    # Write the time-varient metadata
                    if write_meta:
                        for name in self._time_variant_metadata:
                            in_meta  = in_file.variables[name]
                            out_meta = out_file.variables[name]
                            ndims = self._get_var_ndims(in_meta)
                            udidx = in_meta.dimensions.index(self._unlimited_dim)
                            in_slice = [slice(None)] * ndims
                            in_slice[udidx] = slice_step_index
                            out_slice = [slice(None)] * ndims
                            out_slice[udidx] = series_step_index
                            out_meta[tuple(out_slice)] = in_meta[
                                tuple(in_slice)]

                            requested_nbytes = in_meta[:].nbytes
                            self._byte_counts[
                                'Requested Data'] += requested_nbytes
                            actual_nbytes = self.assumed_block_size \
                                * numpy.ceil(requested_nbytes / self.assumed_block_size)
                            self._byte_counts['Actual Data'] += actual_nbytes

                    # Write the time-series variables
                    if write_tser:
                        in_var  = in_file.variables[out_name]
                        out_var = out_file.variables[out_name]
                        ndims   = self._get_var_ndims(in_var)
                        udidx   = in_var.dimensions.index(self._unlimited_dim)
                        in_slice= [slice(None)] * ndims
                        in_slice[udidx] = slice_step_index
                        out_slice = [slice(None)] * ndims
                        out_slice[udidx] = series_step_index
                        out_var[tuple(out_slice)] = in_var[tuple(in_slice)]

                        requested_nbytes = in_file.variables[
                            out_name][:].nbytes
                        self._byte_counts['Requested Data'] += requested_nbytes
                        actual_nbytes = self.assumed_block_size \
                            * numpy.ceil(requested_nbytes / self.assumed_block_size)
                        self._byte_counts['Actual Data'] += actual_nbytes

            # Increment the time-series step index
            series_step_index += 1

            # Now we close the input file as it is no longer needed.
            in_file.close()

            # Print a progress message
            progress_message = "Completed {0:04d} of {1:04d} input files".format(
                                   series_step_index, self.num_input_files)
            self._vprint(progress_message, header=True, verbosity=1)


        self._vprint("Process {0} Completed Work!".format(self._simplecomm.get_rank()), 
                     header=False, verbosity=1)

        # Information
        self._simplecomm.sync()
        if self._simplecomm.is_manager():
            self._vprint(
                'Finished converting time-slices to time-series.', verbosity=1)

        # Finish clocking the entire convert procedure
        if self._timecode: self._timer.stop('Complete Conversion Process')

    def print_diagnostics(self):
        """
        Print out timing and I/O information collected up to this point
        """
        if self._timecode: 
            # Get all totals and maxima
            my_times    = self._timer.get_all_times()
            max_times   = self._simplecomm.allreduce(my_times, op='max')
            my_bytes    = self._byte_counts
            total_bytes = self._simplecomm.allreduce(my_bytes, op='sum')

            # Synchronize
            self._simplecomm.sync()

            # Print timing maxima
            o = self._timer.get_names()
            time_table_str = _pprint_dictionary('TIMING DATA', max_times, order=o)
            if self._simplecomm.is_manager():
                self._vprint(time_table_str, verbosity=0)

            # Convert byte count to MB
            for name in total_bytes:
                total_bytes[name] = total_bytes[name] / float(1024 * 1024)

            # Print byte count totals
            byte_count_str = _pprint_dictionary('BYTE COUNTS (MB)', total_bytes)
            if self._simplecomm.is_manager():
                self._vprint(byte_count_str, verbosity=0)


#==============================================================================
# MultiSpecReshaper Class
#==============================================================================
class MultiSpecReshaper(Reshaper):

    """
    Multiple Slice-to-Series Reshaper class

    This class is designed to deal with lists of multiple 
    Slice2SeriesSpecifiers at a time.  Instead of being instantiated 
    (or initialized) with a single Slice2SeriesSpecifier,
    it takes a dictionary of Slice2SeriesSpecifier objects.
    """

    def __init__(self, specifiers, serial=False, verbosity=1,
                 skip_existing=False, overwrite=False,
                 once=False, simplecomm=None):
        """
        Constructor

        Parameters:
            specifiers (dict): A dict of named Specifier instances, each
                defining an input specification for this reshaper operation.

        Keyword Arguments:
            serial (bool): True or False, indicating whether the operation
                should be performed in serial (True) or parallel
                (False).  The default is to assume parallel operation
                (but serial will be chosen if the mpi4py cannot be
                found when trying to initialize decomposition.
            verbosity(int): Level of printed output (stdout).  A value of 0 
                means no output, and a higher value means more output.  The
                default value is 1.
            skip_existing (bool): Flag specifying whether to skip the generation
                of time-series for variables with time-series files that already
                exist.  Default is False.
            overwrite (bool): Flag specifying whether to forcefully overwrite
                output files if they already exist.  Default is False.
            once (bool): True or False, indicating whether the Reshaper should
                write all metadata to a 'once' file (separately).
            simplecomm (SimpleComm): A SimpleComm object to handle the parallel 
                communication, if necessary
        """

        # Check types
        if not isinstance(specifiers, dict):
            err_msg = "Input must be given in a dictionary of Specifiers"
            raise TypeError(err_msg)
        if type(serial) is not bool:
            err_msg = "Serial indicator must be True or False."
            raise TypeError(err_msg)
        if type(verbosity) is not int:
            err_msg = "Verbosity level must be an integer."
            raise TypeError(err_msg)
        if type(skip_existing) is not bool:
            err_msg = "Skip_existing flag must be True or False."
            raise TypeError(err_msg)
        if type(once) is not bool:
            err_msg = "Once-file indicator must be True or False."
            raise TypeError(err_msg)
        if simplecomm is not None:
            if simplecomm is not isinstance(simplecomm, SimpleComm):
                err_msg = "Simple communicator object is not a SimpleComm"
                raise TypeError(err_msg)

        # Whether to write to a once file
        self._use_once_file = once

        # Whether to write to a once file
        self._skip_existing = skip_existing

        # Whether to write to overwrite output files
        self._overwrite = overwrite

        # Store the list of specifiers
        self._specifiers = specifiers

        # Store the serial specifier
        self._serial = serial

        # Check for a SimpleComm, and if none create it
        if simplecomm is None:
            simplecomm = create_comm(serial=serial)

        # Pointer to its own messenger
        self._simplecomm = simplecomm

        # Store the verbosity
        self._verbosity = verbosity

        # Set the verbose printer
        self._vprint = VPrinter(verbosity=verbosity)

        # Storage for timing data
        self._times = {}

        # Orders for printing timing data
        self._time_orders = {}

        # Storage for all byte counters
        self._byte_counts = {}

    def convert(self, output_limit=0):
        """
        Method to perform each Reshaper's designated operation.

        Loops through and creates each Reshaper, calls each Reshaper's 
        convert() method, and pulls the timing data out for each convert 
        operation.

        Keyword Arguments:
            output_limit (int): Limit on the number of output (time-series) 
                files to write during the convert() operation.  If set
                to 0, no limit is placed.  This limits the number
                of output files produced by each processor in a
                parallel run.
        """
        # Type checking input
        if type(output_limit) is not int:
            err_msg = 'Output limit must be an integer'
            raise TypeError(err_msg)

        # Loop over all specifiers
        for spec_name in self._specifiers:
            if self._simplecomm.is_manager():
                self._vprint('--- Converting Specifier: '
                             + str(spec_name), verbosity=0)

            rshpr = create_reshaper(self._specifiers[spec_name],
                                    serial=self._serial,
                                    verbosity=self._verbosity,
                                    skip_existing=self._skip_existing,
                                    overwrite=self._overwrite,
                                    once=self._use_once_file,
                                    simplecomm=self._simplecomm)
            rshpr.convert(output_limit=output_limit)

            this_times = rshpr._timer.get_all_times()
            self._times[spec_name] = rshpr._simplecomm.allreduce(
                this_times, op='max')
            self._time_orders[spec_name] = rshpr._timer.get_names()
            this_count = rshpr._byte_counts
            self._byte_counts[spec_name] = rshpr._simplecomm.allreduce(
                this_count, op='sum')

            if self._simplecomm.is_manager():
                self._vprint('--- Finished converting Specifier: '
                             + str(spec_name) + os.linesep, verbosity=0)
            self._simplecomm.sync()

    def print_diagnostics(self):
        """
        Print out timing and I/O information collected up to this point
        """
        # Loop through all timers
        for name in self._specifiers:
            if self._simplecomm.is_manager():
                self._vprint('Specifier: ' + str(name), verbosity=0)

            times = self._times[name]
            o = self._time_orders[name]
            times_str = _pprint_dictionary('TIMING DATA', times, order=o)
            if self._simplecomm.is_manager():
                self._vprint(times_str, verbosity=0)

            counts = self._byte_counts[name]
            for name in counts:
                counts[name] = counts[name] / float(1024 * 1024)
            counts_str = _pprint_dictionary('BYTE COUNTS (MB)', counts)
            if self._simplecomm.is_manager():
                self._vprint(counts_str, verbosity=0)
