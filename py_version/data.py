#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# data.py
#
# 2014-01-30 -- created
# 2015-08-22 -- last updated
#
# ~~~~~~~~~
# citation:
# ~~~~~~~~~
# T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J. 
# Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-
# led algorithms for simulating habitats (SPLASH): Modelling radiation evapo-
# transpiration and plant-available moisture, Geoscientific Model Development, 
# 2015 (in progress)

###############################################################################
## IMPORT MODULES:
###############################################################################
import numpy

###############################################################################
## CLASSES
###############################################################################
class DATA:
    """
    Name:     DATA
    Features: This class handles the file IO for reading and writing data.
    """
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Initialization 
    # ////////////////////////////////////////////////////////////////////////
    def __init__(self):
        """
        Name:     DATA.__init__
        Input:    str, input file name (fname)
        Features: Initialize empty class variables
        """
        self.file_name = ""
        self.year = 0
        self.num_lines = 0.
        self.sf_vec = numpy.array([])
        self.tair_vec = numpy.array([])
        self.pn_vec = numpy.array([])
    #
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # Class Function Definitions
    # ////////////////////////////////////////////////////////////////////////
    def read_csv(self, fname, y=-1):
        """
        Name:     DATA.read_csv
        Input:    - str, input CSV filename (fname)
                  - int, year (y)
        Output:   None
        Features: Reads all three daily input variables (sf, tair, and pn) for 
                  a single year from a CSV file that includes a headerline.
        """
        self.file_name = fname
        #
        try:
            data = numpy.loadtxt(fname, 
                                 dtype={'names': ('sf', 'tair', 'pn'),
                                        'formats' : ('f4', 'f4', 'f4')},
                                 delimiter=',',
                                 skiprows=1)
        except IOError:
            print "Could not read input file", fname
        else:
            self.sf_vec = data['sf']
            self.tair_vec = data['tair']
            self.pn_vec = data['pn']
            self.num_lines = data.shape[0]
            #
            if y == -1:
                if data.shape[0] == 366:
                    self.year = 2000
                elif data.shape[0] == 365:
                    self.year = 2001
            else:
                self.year = y
    #
    def read_txt(self, fname, var, y=-1):
        """
        Name:     DATA.read_txt
        Input:    - str, input text file (fname)
                  - str, variable name (i.e., 'pn', 'sf', 'tair')
                  - int, year (y)
        Output:   None.
        Features: Reads plain text file (no header) into one of daily input
                  arrays.
        """
        # Add filename to list:
        if not isinstance(self.file_name, list):
	        self.file_name = []
        self.file_name.append(fname)
        #
        try:
            data = numpy.loadtxt(fname, dtype='f4')
        except IOError:
            print "Could not read input file", fname
        else:
            if var == 'sf':
                self.sf_vec = data
            elif var == 'pn':
                self.pn_vec = data
            elif var == 'tair':
                self.tair_vec = data
            else:
                print 'Variable type not recognized!'
            #
            # Add line numbers to list:
            if not isinstance(self.num_lines, list):
                self.num_lines = []
            self.num_lines.append(data.shape[0])
            #
            if y == -1:
                if data.shape[0] == 366:
                    self.year = 2000
                elif data.shape[0] == 365:
                    self.year = 2001
            else:
                self.year = y
