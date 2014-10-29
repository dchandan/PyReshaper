#!/usr/bin/env python

from distutils.core import setup

setup(name='PyReshaper_n4p',
      version='0.9.0',
      description='Python Time-Slice to Time-Series NetCDF Converter',
      author='Kevin Paul',
      author_email='kpaul@ucar.edu',
      url='https://wiki.ucar.edu/display/~kpaul/PyReshaper',
      packages=['pyreshaper_n4p'],
      scripts=['scripts/slice2series-n4p'],
      requires=['netCDF4', 'mpi4py']
     )
