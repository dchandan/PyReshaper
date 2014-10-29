'''
Class to manage the data conversion process for a given stream of data.
This class is responsible for performing the conversion operation, assuming
all of the data (files) conform to a common specification (i.e., 1 and only 1
specifier).  The operational assumption is that one single specification
will be given one messenger to handle all communication in a single group
dedicated to the single specification.

________________________
Created on May 26, 2014

@author: Kevin Paul <kpaul@ucar.edu>
'''

from specification import Specifier
from messenger import Messenger


class Converter(object):
    '''
    Base class for all data converters.  Defines a simple interface with
    a constructor and a 'convert()' method.  As input it takes a single
    Specifier and a single Messenger.
    '''

    def __init__(self, spcfr, msngr):
        '''
        Constructor
        
        @param spcfr  A Specifier object indicating what operation
        '''
