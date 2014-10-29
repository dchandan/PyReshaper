'''
This is the main reshaper module.  This is the manager class for all of the
reshaper operations.  It is responsible for creating as many conversion
objects as needed for parallel operation and running each conversion
utility.

________________________
Created on Apr 30, 2014

@author Kevin Paul <kpaul@ucar.edu>
'''

from specification import Specifier, Slice2SeriesSpecifier
from messenger import create_messenger
from timekeeper import TimeKeeper

import Nio
import numpy
import itertools


#==============================================================================
# Reshaper Class
#==============================================================================
class Reshaper(object):
    '''
    This is the manager class for the reshaper.  This class creates and
    controls all of the multiple sub-objects needed to perform the many
    data conversion operations specified by each individual Specifier object.
    '''

    def __init__(self, spec_list, serial=False, debug=1):
        '''
        Constructor

        @param spec_list  An list of Specifier class instantiations that define
                          the information needed for each type of Reshaper
                          operation to be performed.  (That is, a
                          Slice2SeriesSpecifier will invoke the creation of a
                          matching Slice2SeriesReshaper object.)

        @param serial     True or False, indicating whether the Reshaper object
                          should perform its operation in serial (True) or
                          parallel (False).  Default is False (parallel).

        @param debug      Integer value indicating the debug level of output
                          that should be generated during the run.  Zero
                          indicates no output, and a larger value indicates
                          more detailed output.  The default is 1.
        '''

        # Type checking the Specifier list
        if (not isinstance(spec_list, list)):
            err_msg = 'The list of Specifier objects must be of list type.'
            raise TypeError(err_msg)
        for spec in spec_list:
            if (not isinstance(spec, Specifier)):
                err_msg = 'Each object in the Specifier list must be of ' + \
                          'Specifier type.'
                raise TypeError(err_msg)

        ## The list of input specifications
        self._specifier_list = spec_list

        if (type(serial) is not bool):
            err_msg = 'Serial indicator must be True or False.'
            raise TypeError(err_msg)
        if (type(debug) is not int):
            err_msg = 'Debug level indicator must be and integer.'
            raise TypeError(err_msg)

        ## Master messenger utility (for serial or parallel communications)
        self._messenger = create_messenger(serial=serial, debug=debug)
        self._messenger.print_once('Initializing Reshaper', level=1)

        ## Time keeper (for timing/clocking internal operations)
        self._timer = TimeKeeper(messenger=self._messenger)

        ## Limit on the number of output files to write (for debugging)
        self._output_limit = 0

        # Validate the user input data
        self._timer.start('Specifier Validation')
        self._specifier.validate()
        self._timer.stop('Specifier Validation')
        self._messenger.debug_print('Specifier validated')

    def convert(self):
        '''
        Method to perform the Reshaper's designated operation.
        '''
        pass


#==============================================================================
# Slice2SeriesReshaper Class
#==============================================================================
class Slice2SeriesReshaper(Reshaper):
    '''
    This is the main reshaper class.  This is the class that defines how
    the time-slice to time-series reshaping operation is to be performed.
    '''

    def __init__(self, spec, serial=False):
        '''
        Constructor

        @param specifier  An instance of the Specifier class, defining the
                          input specification for this reshaper operation.

        @param serial     True or False, indicating whether the operation
                          should be performed in serial (True) or parallel
                          (False).  The default is to assume parallel operation
                          (but serial will be chosen if the mpi4py cannot be
                          found when trying to initialize decomposition.
        '''
        # Type checking (or double-checking)
        if (not isinstance(spec, Slice2SeriesSpecifier)):
            err_msg = "Slice2SeriesReshaper requires a Slice2SeriesSpecifier" \
                    + " as input."
            raise TypeError(err_msg)

        # Call the base-class constructor
        super(Slice2SeriesReshaper, self).__init__(spec, serial=serial)

        # Setup PyNIO options (including disabling the default PreFill option)
        opt = Nio.options()
        opt.PreFill = False

        # Determine the Format and CompressionLevel options
        # from the NetCDF format string in the Specifier
        if (self._specifier.netcdf_format == 'netcdf'):
            opt.Format = 'Classic'
        elif (self._specifier.netcdf_format == 'netcdf4'):
            opt.Format = 'NetCDF4Classic'
            opt.CompressionLevel = 0
        elif (self._specifier.netcdf_format == 'netcdf4c'):
            opt.Format = 'NetCDF4Classic'
            opt.CompressionLevel = 1
        self._nio_options = opt
        self._messenger.debug_print('PyNIO options set')

        # Open all of the input files
        self._timer.start('Open Input Files')
        self._input_files = []
        for filename in self._specifier.input_file_list:
            self._input_files.append(Nio.open_file(filename, "r"))
        self._timer.stop('Open Input Files')
        self._messenger.debug_print('Input files opened')

                # Validate the input files themselves
        self._timer.start('Input File Validation')
        self._validate_input_files()
        self._timer.stop('Input File Validation')
        self._messenger.debug_print('Input files validated')

        # Sort the input files by time
        self._timer.start('Sort Input Files')
        self._sort_input_files_by_time()
        self._timer.stop('Sort Input Files')
        self._messenger.debug_print('Input files sorted')

        # Retrieve and sort the variables in each time-slice file
        # (To determine if it is time-invariant metadata, time-variant
        # metadata, or if it is a time-series variable)
        self._timer.start('Sort Variables')
        self._sort_variables()
        self._timer.stop('Sort Variables')
        self._messenger.debug_print('Variables sorted')

    def _validate_input_files(self):
        '''
        Perform validation of input data files themselves.  We check
        the file contents here, assuming that the files are already open.
        '''

        # Make a pass through each file to make sure there is a 'time'
        # dimension and a 'time' variable
        for i in range(len(self._input_files)):
            ifile = self._input_files[i]
            if ('time' not in ifile.dimensions):
                err_msg = 'Time dimension not found in file (' \
                        + self._specifier.input_file_list[i] + ')'
                raise LookupError(err_msg)
            if ('time' not in ifile.variables):
                err_msg = 'Time variable not found in file (' \
                        + self._specifier.input_file_list[i] + ')'
                raise LookupError(err_msg)

        # Make sure that the list of variables in each file is the same
        variables = self._input_files[0].variables
        var_names = set(variables.keys())
        for ifile in self._input_files[1:]:
            var_names_next = set(ifile.variables.keys())
            if (len(var_names - var_names_next) != 0):
                err_msg = "All input files do not contain the variables."
                raise ValueError(err_msg)

    def _sort_input_files_by_time(self):
        '''
        Internal method for sorting the input files by time, and check to
        make sure that all of the times spanning across each file do not
        overlap with each other (i.e., that the times across all files are
        monotonicly increasing).  We also store the array of times across
        all input files for future writing.

        @note Currently, this method assumes that all of the input files
              have the same 'time:units' attribute, such that all time variable
              values are measured from the same date-time.  When this is true,
              we do not need to consider the value of the 'time:units'
              attribute itself.  If this assumption is not true, then we need
              to consider the 'time:units" attribute of each file, together
              with that file's time variable values.  This is what the CFTime
              class is intended to do.
        '''

        # Get the time attributes (for convenience) and, for each file,
        # add the times to a list.  (Each file will have an array of times
        # associated with it.  Each array will be added to a list, such
        # that the outer-most list contains an array for each input file)
        time_values = []
        for ifile in self._input_files:
            time_values.append(ifile.variables['time'].get_value())

        # Determine the sort order based on the first time in the time values
        order = range(len(self._input_files))
        new_order = sorted(order, key=lambda i: time_values[i][0])

        # Re-order the list of input files and filenames
        new_file_list = [None] * len(new_order)
        new_filenames = [None] * len(new_order)
        new_values = [None] * len(new_order)
        for i in order:
            new_file_list[i] = self._input_files[new_order[i]]
            new_filenames[i] = self._specifier.input_file_list[new_order[i]]
            new_values[i] = time_values[new_order[i]]

        # Save this data in the new orders
        self._input_files = new_file_list
        self._input_filenames = new_filenames

        # Now, check that the largest time in each file is less than the
        # smallest time in the next file (so that the time spans of each file
        # do not overlap)
        for i in order[:-1]:
            if (new_values[i][-1] >= new_values[i + 1][0]):
                err_msg = 'Times in input files ' + str(new_filenames[i]) \
                        + ' and ' + str(new_filenames[i + 1]) + ' appear to ' \
                        + 'overlap.'
                raise ValueError(err_msg)

        # Now that this is validated, let's string together the numpy array
        # of all times (using the new_values array)
        self._all_time_values = \
            numpy.fromiter(itertools.chain.from_iterable(new_values),
                           dtype='float')

    def _sort_variables(self):
        '''
        Internal method for sorting the variables that exist in each
        time-slice file.  This method determines if each variable is to be
        treated as time-invariant metadata, time-variant metadata (user
        defined), or time-series variables.  All metadata is written to
        every time-series file, and any time-series variable is written to
        its own file.  The time-variant metadata variables are determined
        by user input, and are contained in the Specifier data member,
        Specifier.time_variant_metadata.
        '''

        # Initialize the list of variable names for each category
        self._time_variant_metadata = []
        self._time_invariant_metadata = []
        self._time_series_variables = []

        # Categorize each variable (only looking at first file)
        variables = self._input_files[0].variables
        for var_name in self._input_files[0].variables.keys():
            if ('time' not in variables[var_name].dimensions):
                self._time_invariant_metadata.append(var_name)
            elif (var_name in self._specifier.time_variant_metadata):
                self._time_variant_metadata.append(var_name)
            else:
                self._time_series_variables.append(var_name)

        # Debug output
        self._messenger.debug_print('Time-Invariant Metadata: ' + \
                                  str(self._time_invariant_metadata))
        self._messenger.debug_print('Time-Variant Metadata: ' + \
                                  str(self._time_variant_metadata))
        self._messenger.debug_print('Time-Series Variables: ' + \
                                  str(self._time_series_variables))

    def convert(self):
        '''
        Method to perform the Reshaper's designated operation.  In this case,
        convert a list of time-slice files to time-series.
        '''
        self._messenger.debug_print('Converting time-slices to time-series')

        # For data common to all input files, we reference only the first
        ref_infile = self._input_files[0]

        # Store the common dimensions and attributes for each file
        # (taken from the first input file in the list)
        common_dims = ref_infile.dimensions
        common_atts = ref_infile.attributes

        # Partition the time-series variables across all processors
        tsv_names_local = self._messenger.partition(self._time_series_variables)
        if (self._time_series_limit > 0):
            tsv_names_local = tsv_names_local[0:self._time_series_limit]

        # Print partitions for all ranks
        dbg_msg = 'Local time-series variables are ' \
                + str(tsv_names_local)
        self._messenger.debug_print(dbg_msg, allranks=True)

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
        for out_name in tsv_names_local:
            out_filename = self._specifier.output_file_prefix \
                         + out_name + self._specifier.output_file_suffix
            self._messenger.debug_print('Creating output file: ' + out_filename)

            # Open each output file and create the dimensions and attributes
            self._timer.start('Open Output Files')
            out_file = Nio.open_file(out_filename, 'w',
                                     options=self._nio_options)
            for att_name, att_val in common_atts.items():
                setattr(out_file, att_name, att_val)
            for dim_name, dim_val in common_dims.items():
                if (dim_name == 'time'):
                    out_file.create_dimension(dim_name, None)
                else:
                    out_file.create_dimension(dim_name, dim_val)
            self._timer.stop('Open Output Files')

            # Create the time-invariant metadata variables
            self._timer.start('Create Time-Invariant Metadata')
            for name in self._time_invariant_metadata:
                in_var = ref_infile.variables[name]
                out_var = out_file.create_variable(name,
                    in_var.typecode(), in_var.dimensions)
                for att_name, att_val in in_var.attributes.items():
                    setattr(out_var, att_name, att_val)
            self._timer.stop('Create Time-Invariant Metadata')

            # Create the time-variant metadata variables
            self._timer.start('Create Time-Variant Metadata')
            for name in self._time_variant_metadata:
                in_var = ref_infile.variables[name]
                out_tvm_vars[name] = out_file.create_variable(name,
                    in_var.typecode(), in_var.dimensions)
                for att_name, att_val in in_var.attributes.items():
                    setattr(out_tvm_vars[name], att_name, att_val)
            self._timer.stop('Create Time-Variant Metadata')

            # Create the time-series variable itself
            self._timer.start('Create Time-Series Variables')
            in_var = ref_infile.variables[out_name]
            out_var = out_file.create_variable(out_name,
                in_var.typecode(), in_var.dimensions)
            self._timer.stop('Create Time-Series Variables')

            # Append the output file to list
            out_files[out_name] = out_file

        # Now that each output file has been created, start writing the data
        # (Looping over output file index, which is common in name lists)
        for out_name, out_file in out_files.items():

            dbg_msg = 'Writing output file for variable: ' + out_name
            self._messenger.debug_print(dbg_msg, allranks=True)

            # Create the attributes of the time-series variable
            in_var = ref_infile.variables[out_name]
            out_var = out_file.variables[out_name]
            for att_name, att_val in in_var.attributes.items():
                setattr(out_var, att_name, att_val)

            # Write the time-invariant metadata
            self._timer.start('Write Time-Invariant Metadata')
            for name in self._time_invariant_metadata:
                in_meta = ref_infile.variables[name]
                out_meta = out_file.variables[name]
                if in_meta.rank > 0:
                    out_meta[:] = in_meta[:]
                else:
                    out_meta.assign_value(in_meta.get_value())
            self._timer.stop('Write Time-Invariant Metadata')

            # Write each time-variant variable
            t_index = 0
            for in_file in self._input_files:

                # Write the time-varient metadata
                self._timer.start('Write Time-Variant Metadata')
                for name in self._time_variant_metadata:
                    in_meta = in_file.variables[name]
                    out_meta = out_file.variables[name]
                    out_meta[t_index] = in_meta[:]
                self._timer.stop('Write Time-Variant Metadata')

                # Write the time-series variables
                self._timer.start('Write Time-Series Variables')
                out_var[t_index] = in_file.variables[out_name][:]
                self._timer.stop('Write Time-Series Variables')

                # Increment the time index
                t_index += 1

            # Close the output file
            self._timer.start('Close Output Files')
            out_file.close()
            self._timer.stop('Close Output Files')

        self._messenger.debug_print(
            'Finished converting time-slices to time-series.')

    def print_diagnostics(self):
        '''
        Print out timing and I/O information collected up to this point
        in the Reshaper's operation.
        '''
        self._messenger.debug_print('Printing diagnostics')

        # Before printing anything, synchronize across processors
        self._messenger.sync()

        # Print only on master
        self._messenger.master_print(self._timer.print_times())
