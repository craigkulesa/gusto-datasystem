"""GUSTO Level 2 Pipeline preprocessor module
"""

__version__ = '0.1'
__date__ = '9/19/2024'
__updated__ = '20240919'   # automatically updated when saving
__author__ = 'V. Tolls, CfA | Harvard & Smithsonian'


import astropy.units as u
import astropy.constants as c
from astropy.io import ascii
import pandas as pd
import numpy as np
from numpy.lib import recfunctions as rfn
from dataclasses import dataclass, asdict
from .GL09PParameters import *
import os

# primary function
class GL09PProcessControl:
    """
    The purpose of the process control is to
    1.) read a script controlling the GL09P
        script synthax is simple:
        module = 0,1 or True/False
    2.) analyze the input script for what to run and what not
    3.) return the array for process control to the pipeline

    filename: GL09P_control_script.txt
    location of script: 


    """

    parlist = []  # Class attribute to track instance names.

    def __init__(self, scriptfile=None, verbose=False):
        """Initializing the process control class

        Parameters
        ----------
        scriptfile : str
            optional line-specific configuration file

        Returns
        -------
        line configuration structure
        """
        # exec("self.{} = '{}'".format(vname, vvalue))

        self.padd('__date__', __date__, '', 'created')
        self.padd('__version__', __version__, '', 'class version')
        self.padd('__updated__', __updated__, '', 'last update')
        #self.padd('__pyfile__', __pyfile__, '', 'file containing this class')
        #self.padd('__path__', __path__, '', 'path to this file')
        self.padd('__author__', __author__, '', 'author')

        if scriptfile is None:
            ifile = '../Data/GL09P_setup_pars.txt'
        else:
            ifile = scriptfile

        if scriptfile is None:
            #wdir = os.getcwd()
            wdir = os.path.dirname(os.path.abspath(__file__))
            ifile = wdir+'/Data/GL09P_control_script.txt'
            print('No GUSTO Pipeline Control Script File file provided, \n' +
                  'using default file: ', ifile, '\n',
                  ' please save in your working directory.\n')

            # open file save window to save the template file somewhere else
        else:
            ifile = parfile
        self.padd('__parfile__', ifile, '', 'file containing parameters')

    def getProcInfo(self):
        """
        Function parsing and returning a structure with the process information 
        to control the pipeline execution, e.g., which lines to be processed and 
        which pipeline function like baseline fit to be executed.
        """
        pcInfo = {}

    print('GL09ProcessControl function called.')

    def padd(self, pname, pvalue, punit='', pinfo='', level=0, levelnames=''):
        """Adding a parameter

        Input:
        ------

        """
        self.parlist.append(pname)
        setattr(self, pname, Par(pname, pvalue, punit, pinfo))

    def cadd(self, cname):
        """Adding an empty ReCclass

        Input:
        ------

        """
        self.parlist.append(cname)
        setattr(self, cname, RecClass())

    def getParList(self):
        # return self.parlist
        return self.__dict__


@dataclass
class Par(GL09PParameters):
    name: str
    value: float
    unit: str
    info: str


class RecClass:
    """Class to create nested parameter classes
    """

    def __init__(self):
        pass

    def addclp(self, clname):
        """Add a parameter class with name clname
        """
        setattr(self, clname, Recpars())

    def padd(self, pname, pvalue, punit='', pinfo=''):
        """Add a parameter to the class
        """
        setattr(self, pname, Par(pname, pvalue, punit, pinfo))


class Recpars:

    def __init__(self):
        self.par1 = 20

    def padd(self, pname, pvalue, punit='', pinfo=''):
        """Add a parameter to the class
        """
        setattr(self, pname, Par(pname, pvalue, punit, pinfo))

    def addclp(self, clname):
        """Add a parameter class with name clname
        """
        setattr(self, clname, Recpars())


if __name__ == '__main__':

    # initialize the process control class
    prcntl = GL09ProcessControl()

    # retrieve the process control information
    prinfo = prcntl.getProcInfo()

    # inspect the process information
    print(prinfo)
