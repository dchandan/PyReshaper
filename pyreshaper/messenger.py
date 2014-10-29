'''
This is the parallel messenger class to provide basic load-balancing and
decomposition routines for the PyReshaper.  This encapsulates all of the
mpi4py calls as needed.  It also has the ability to detect if the mpi4py
module is present.  If not, it throws an exception.

_______________________
Created on Apr 30, 2014

@author: Kevin Paul <kpaul@ucar.edu>
'''


#==============================================================================
# create_messenger factory function
#==============================================================================
def create_messenger(serial=False, debug=1):
    '''
    This is the factory function the creates the necessary messenger
    object for parallel (or serial) operation.  The type must be specified
    by the user with the 'serial' argument, as there is not portable way
    of determining if the run should be assumed to be serial or parallel
    from the environment.

    @param serial  True or False, indicating whether the serial or parallel
                   messenger utility object should be constructed and
                   returned.  DEFAULT: False (parallel operation)

    @param debug   Integer indicating level of debug output desired during
                   operation.  Zero indicates no output, and a higher value
                   indicates more detailed output.  Defaults to 1.

    @return  A messenger object.
    '''

    # Check type
    if (type(serial) is not bool):
        err_msg = "The serial argument must be a bool."
        raise TypeError(err_msg)

    # Construct and return the desired decomp utility
    if (serial):
        return Messenger(debug=debug)
    else:
        return MPIMessenger(debug=debug)


#==============================================================================
# Messenger Base Class
#==============================================================================
class Messenger(object):
    '''
    This is the base class for serial/parallel messenger utility.  This defines
    serial operation, and has no dependencies.  These methods are reimplemented
    in the derived class for parallel decomposition.

    Each Messenger is aware of the 'world' group in which it has been created.
    Each Messenger can also be assigned to a subgroup, which is indicated by
    an integer index.  Operations can be performed across the 'world' group,
    subgroups, or on a specific Messenger.
    '''

    def __init__(self, debug=1):
        '''
        Constructor

        @param debug   Integer indicating level of debug output desired during
                       operation.  Zero indicates no output, and a higher value
                       indicates more detailed output.  Defaults to 1.
        '''
        # Type checking
        if (type(debug) is not int):
            err_msg = "The debug level argument must be an integer."
            raise TypeError(err_msg)

        ## Current debug level for output
        self._debug = debug

        ## The type of decomp utility constructed
        self.messenger_type = 'serial'

        ## The rank of this processor in the world group
        self._world_rank = 0

        ## The number of processors in the world group
        self._world_size = 1

        ## Whether this is the master process/rank for the world group
        self._is_world_master = True

        ## This processor's subgroup identifier
        self._subgroup = 0

        ## The rank of this processor in the subgroup
        self._subgroup_rank = 0

        ## The number of processors in the subgroup
        self._subgroup_size = 1

        ## Whether this is the master process/rank for the subgroup
        self._is_subgroup_master = True

    def assign_to_subgroup(self, subgroup):
        '''
        This method assigns this Messenger to a given subgroup.
        '''
        if (type(subgroup) is not int):
            err_msg = 'Subgroup identifier must be an integer'
            raise TypeError(err_msg)

        # Does nothing in serial
        return

    def partition(self, global_list, weights=None, world=False):
        '''
        Takes a list and returns a subset of the list, corresponding to the
        subset for which the current process is given responsibility.

        For serial execution, this returns the global variables dictionary.

        @param global_list  A list of objects to be partitioned over all of the
                            processors in the current group.

        @param weights      Weights corresponding to each item in the global
                            list.  In serial, this is ignored.

        @param world        True or False, indicating whether this list should
                            be partitioned over the entire world group (True),
                            or over the assigned subgroup (False).

        @return The local list of objects to be the responsibility of the
                current processor.
        '''
        if (not isinstance(global_list, list)):
            err_msg = "List of global variables has wrong type"
            raise TypeError(err_msg)

        return global_list

    def sync(self, world=False):
        '''
        A wrapper on the MPI Barrier method.  Forces execution to wait for
        all other processors to get to this point in the code.

        Does nothing in serial.

        @param world  True or False, indicating whether this barrier should
                      be applied to the entire world group (True), or just
                      to the assigned subgroup (False).
        '''
        return

    def is_master(self, world=False):
        '''
        Returns True or False depending on whether this rank is the master
        rank (e.g., rank 0).  In serial, always returns True.

        @param world  True or False, indicating whether tto return if this
                      processor is the master for the entire world group
                      (True), or just for the assigned subgroup (False).

        @return True or False depending on whether this rank is the master
                rank in this group (e.g., rank 0).
        '''
        if (world):
            return self._is_world_master
        else:
            return self._is_subgroup_master

    def get_rank(self, world=False):
        '''
        Returns the number associated with the local MPI processor/rank of
        this group.  In serial, it always returns 0.

        @param world  True or False, indicating whether to return the world
                      group rank (True) or the assigned subgroup rank (False).

        @return The integer ID associated with the local group or subgroup
        '''
        if (world):
            return self._world_rank
        else:
            return self._subgroup_rank

    def get_size(self, world=False):
        '''
        Returns the number of ranks in this group (i.e., in this MPI
        communicator.)  In serial, it always returns 1.

        @param world  True or False, indicating whether to return the world
                      group size (True) or the assigned subgroup size (False).

        @return The integer number of ranks in the group or subgroup
        '''
        if (world):
            return self._world_rank
        else:
            return self._subgroup_rank

    def sum(self, data, world=False):
        '''
        Sums all of the data across all processors in the group or subgroup.

        In serial, just sums the data passed to this method.

        @param data The data with values to be summed (an iterable)

        @param world  True or False, indicating whether to sum data over the
                      entire world (True) or over only the assigned subgroup
                      (False).

        @return The sum of the data values
        '''
        return reduce(lambda x, y: x + y)

    def find_max(self, data, world=False):
        '''
        Finds the element-by-element maximum of a set of data across all
        processors in the group.  In serial, just returns the data passed in.

        @param data The data in which the maximum is to be found

        @param world  True or False, indicating whether to find the maximum
                      in the entire world (True) or only in the assigned
                      subgroup (False).

        @return A list with the maximum value of each element in the data
        '''
        return data

    def find_min(self, data, world=False):
        '''
        Finds the element-by-element minimum of a set of data across all
        processors in the group.  In serial, just returns the data passed in.

        @param data The data in which the maximum is to be found

        @param world  True or False, indicating whether to find the minimum
                      in the entire world (True) or only in the assigned
                      subgroup (False).

        @return A list with the minimum value of each element in the data
        '''
        return data

    def print_once(self, output, level=0, world=True):
        '''
        This method prints output to stdout, but only from one processor
        per group or subgroup.  The output is prefixed with a tag indicating
        [group:rank:size] of the world group or the assigned subgroup.

        @param output  A string that should be printed to stdout

        @param level  The debug level to apply to this message.  If this level
                      is less than the Messenger's debug level, this message
                      is printed.

        @param world  Whether to print once from the entire world (True), or
                      once from each subgroup (False).
        '''
        if (self.is_master(world=world)):
            self.print_all(output, level=level, world=world)

    def print_all(self, output, level=0, world=False):
        '''
        This prints a message to stdout from every processor.  The output is
        prefixed with a tag indicating [group:rank:size] of the world group or
        the assigned subgroup.

        @param output A string that should be printed to stdout

        @param level A list of ranks from which output should be printed

        @param world  Whether to prefix the output with world labels or
                      subgroup labels
        '''
        if (level < self._debug):
            prefix = '['
            if (world):
                prefix += 'WORLD:' \
                       + str(self._world_rank) + ':' \
                       + str(self._world_size)
            else:
                prefix = '[' + str(self._subgroup) + ':' \
                       + str(self._subgroup_rank) + ':' \
                       + str(self._subgroup_size)
            prefix += '] '
            print prefix + output


#==============================================================================
# MPIMessenger Class
#==============================================================================
class MPIMessenger(Messenger):
    '''
    This is the parallel-operation class for decomposition/parallel utilities.
    This is derived from the Messenger class, which defines the serial
    operation.

    Each Messenger is aware of the 'world' group in which it has been created.
    Each Messenger can also be assigned to a subgroup, which is indicated by
    an integer index.  Operations can be performed across the 'world' group,
    subgroups, or on a specific Messenger.
    '''

    def __init__(self, debug=1):
        '''
        Constructor

        @param debug   Integer indicating level of debug output desired during
                       operation.  Zero indicates no output, and a higher value
                       indicates more detailed output.  Defaults to 1.
        '''
        # Type checking
        if (type(debug) is not int):
            err_msg = "The debug level argument must be an integer."
            raise TypeError(err_msg)

        ## Current debug level for output
        self._debug = debug

        ## The type of decomp utility constructed
        self.messenger_type = 'parallel'

        # Try to import the MPI module and get the world communicator
        try:
            from mpi4py import MPI
            self._world_comm = MPI.COMM_WORLD
            self._subgroup_comm = MPI.COMM_WORLD
        except:
            raise ImportError('Failed to import MPI.')

        ## The rank of this processor in the world group
        self._world_rank = self._world_comm.Get_rank()

        ## The number of processors in the world group
        self._world_size = self._world_comm.Get_size()

        ## Whether this is the master process/rank for the world group
        self._is_world_master = (self._world_rank == 0)

        ## This processor's subgroup identifier
        self._subgroup = 0

        ## The rank of this processor in the subgroup
        self._subgroup_rank = self._world_rank

        ## The number of processors in the subgroup
        self._subgroup_size = self._world_size

        ## Whether this is the master process/rank for the subgroup
        self._is_subgroup_master = self._is_world_master


    def _unweighted_partition(self, main_list, num_sublists, list_idx):
        '''
        This is private helper function designed to implement an unweighted
        partitioning algorithm useful for subdividing a list into any
        number of sublists, returning the 'list_idx' sublist.  Unweighted
        partitioning equates to equal-number partitioning.

        @param main_list     The main list to the partitioned

        @param num_sublists  The number of sublists to be created

        @param list_idx      The index of the list to be returned

        @return  The 'list_idx's corresponding sublist.
        '''

    def _weighted_partition(self, main_list, num_sublists, list_idx, weights):
        '''
        This is private helper function designed to implement a weighted
        partitioning algorithm useful for subdividing a list into any
        number of sublists, returning the 'list_idx' sublist.

        @param main_list     The main list to the partitioned

        @param num_sublists  The number of sublists to be created

        @param list_idx      The index of the list to be returned

        @param weights       The weight of each item in the main list.  The
                             partitioning algorithm creates sublists such that
                             the sum of all of the weights in each sublist
                             should as close to each other as possible.

        @return   The 'list_idx's corresponding sublist.
        '''

    def partition(self, global_vars_names):
        '''
        Takes a list of variable names and extracts a subset of the
        variables to be used by the local processor.

        @param global_vars_names  A list of global time-series variable names

        @return A list of local time-series variable names (a subset of the
                global variables list)
        '''
        if (not isinstance(global_vars_names, list)):
            err_msg = "List of global variables has wrong type"
            raise TypeError(err_msg)

        rank = self._mpi_rank
        size = self._mpi_size

        local_vars_names = global_vars_names.keys()[rank::size]

        return local_vars_names

    def sync(self):
        '''
        A wrapper on the MPI Barrier method.  Forces execution to wait for
        all other processors to get to this point in the code.
        '''
        self._mpi.COMM_WORLD.Barrier()

    def sum(self, data):
        '''
        Sums all of the data across all processors.

        @param data The data with values to be summed (a dictionary)

        @return The sum of the data values
        '''
        total = Messenger.sum(self, data)
        return self._mpi.COMM_WORLD.allreduce(total, op=self._mpi.SUM)

    def max(self, data):
        '''
        Finds the maximum value of each element in a set of data,
        across all processors.

        @param data The list of data in which the maximum is to be found

        @return A list with the maximum value of each element in the data
        '''
        return self._mpi.COMM_WORLD.allreduce(data, op=self._mpi.MAX)
