"""
GL2PConfigClass:
GUSTO Pipeline configuration management functions.
The functions provided by this class include:
- setting the input data path
- setting the output data file path
- setting/getting the GUSTO general instrument/observation parameters
etc.
"""

__date__ = '9/19/2024'
__updated__ = '20240919'
__version__ = '0.1'
__author__ = 'V. Tolls, CfA | Harvard & Smithsonian'


import time
import configparser
import os
import sys
import numpy as np
import pkg_resources
from os.path import expanduser, expandvars
#from joblib import Parallel, delayed, cpu_count
import logging
log = logging.getLogger(__name__)

cfg_file0 = pkg_resources.resource_filename('gustoL09P', 'Data/GL09P_setup_pars.txt')


def getRange(icpar, dtype='float', endpoint=True):
    """Function reading an input str and convert ranges to 
    an array of floats. 
    Uses getValues() and returns a numpy float array if <3 parameters are
    provided. 
        
    Parameters
    ----------
    icpar: str
        parameter str like '[start, stop, step]'
    """
    pars = getValues(icpar)
    
    if pars.size < 3:
        return pars
    else:
        # create a grid from the parameters
        if endpoint:
            return np.arange(pars[0], pars[1]+pars[2], pars[2])
        else:
            return np.arange(pars[0], pars[1], pars[2])


def getValues(icpar, dtype='float'):
    """Function converting str data to numpy float array.
        
    Parameters
    ----------
    icpar: str
        parameter str like '[start, stop, step]' or '3.45'
    """
    icpars = icpar.replace('[','').replace(']','').replace(' ','').split(',')
    ipars = np.array((icpars), dtype=dtype)
    
    return ipars


class GL09PConfigClass:
    """GUSTO Level 2 Configuration class
    """

    def __init__(self, configFile=None, verbose=False):
        """Initializing the GUSTO configuration class.
        The configuration file, e.g., GL2P_setup_pars.txt, can be provided via
        command line input, or exists in the standard directory for GUSTO, <user-home>/.gusto, 
        or can be loaded from the package data directory. An example file is provided
        in the package data directory: Data/GL2P_setup_pars.txt

        Parameters
        ----------
        config_file: str
            Full path to an existing configuration file to control the simulator
        verbose : bool
            enable/disable printing output to stdout.

        """

        print('Init configFile: ', configFile)
        print()
        if configFile == '':
            # check if there is a config File in the GUSTO config directory
            uhome = expanduser("~")
            infile = cfg_file0   #uhome+'/.gusto/GL2P_setup_pars.txt'

            # check if standard config file exists
            ret = os.path.exists(infile)
            if ret:
                cffile = infile
            else:
                print('Error: no GUSTO Pipeline configuration file available.\n' +
                      'please create a configuration file in ~/.gusto/GL2Pconfig.txt\n' +
                      'An example file is in the package data directory.' +
                      'Or provide a config with full file path when starting the pipeline.')
                sys.exit(1)
        else:
            ret = os.path.exists(configFile)
            if ret:
                cffile = configFile
            else:
                print('Error: no GUSTO Pipeline configuration file %s does not exist.\n' % (configFile) +
                      'Please provide a valid path to an existing config file!')
                sys.exit(1)
        self.cffile = cffile
        # self.dtime = datetime.now().strftime("%Y%m%d_%H%M%S.%f").rstrip('0')
        self.dtstr = datetime.now().strftime("%Y%m%d_%H%M%S")

    def getConfigInfo(self, verbose=False):
        """Function returning the configuration information dictionary.

        Parameters
        ----------
        verbose : bool
            enable/disable printing output to stdout.

        """
        # read the configuration file, but parsing is done in another function
        self.cfg = self.readConfigFile(self.cffile, verbose=verbose)
        if verbose:
            print('\ngetConfigInfo: \n', self.cfg)

        return self.cfg

    def readConfigFile(self, cffile=None, verbose=False):
        """Function parsing the GUSTO pipeline configuration file.
        Note, the synthax of the configuration file is fixed. If the file
        cannot be read, the program exits with an error message.
        
        Parameters
        ----------
        cffile: str
            Full path to an existing configuration file to control the simulator
        verbose : bool
            enable/disable printing output to stdout.


        """
        print('readConfigFile', verbose)
        config = configparser.ConfigParser()

        cfg = config.read(cffile)
        if verbose:
            print('readConfigFile: ', cfg)

        cfg_sect = config.sections()

        # GUSTO_Directories
        # gpars = config.options("GUSTO_Parameters")
        defpars = config["DEFAULT"]
        defdict = {}
        for key in defpars.keys():
            defdict[key] = defpars[key]

        # GUSTO_Directories
        # gpars = config.options("GUSTO_Parameters")
        gdirs = config["GUSTO_Directories"]
        gddict = {}
        for key in gdirs.keys():
            if 'work' in key:
                # add the date/time stamp to the working directory
                print('work dir (before): ', gdirs[key])
                gddict[key] = gdirs[key]+'_'+self.dtstr
                print('work dir (after): ', gdict[key])
            else:
                print('other: ', gdirs[key])
                gddict[key] = gdirs[key]

        # GUSTO_Parameters
        # gpars = config.options("GUSTO_Parameters")
        gpars = config["GUSTO_Parameters"]
        gpdict = {}
        for key in gpars.keys():
            gpdict[key] = gpars[key]

        # create a result dictionary
        sdict = {}
        sdict['default'] = defdict
        sdict['gdirs'] = gddict
        sdict['gparams'] = gpdict
        self.sdict = sdict

        return sdict

    def applyConfig(self, verbose=False):
        """Function checking and setting the directories needed for the pipeline 
        to work. For data access, the directory information is retrieved on the fly.
        This function performs check if the directories are available and, e.g., in 
        case of the work directory, it creates it if it does not exists. 
        
        Parameters
        ----------
        verbose : bool
            enable/disable printing output to stdout.


        """
        sdict = self.sdict

        # work directory
        if not os.path.exists(sdict['work_dir']):
            if verbose:
                print('GUSTO Pipeline working directory does not exist. Creating it!')
                os.path.mkdir()


def GL09P_CreateConfigFile(configFileName=None):
    """Convenience function to create a standard GUSTO Pipeline configuration file.
    
    WILL BE DEPRECATED => example file in package Data directory

    Parameters
    ----------
    configFileName : str
        Full path to the to be created configuration file.
    """

    if configFileName is None:
        print('Error: No configFileName provided. Cannot create file!')
        sys.exit(1)

    config = configparser.ConfigParser(allow_no_value=True)

    config['DEFAULTS'] = {'number_of_cpus': 4}

    config['GUSTO_Directories'] = \
        {'# with full path!': None,
         'configDir': "",
         'gustoIncomingDataDir': "",
         'gustoCIIDataDir': "",
         'gustoNIIDataDir': "",
         'gustoOIDataDir': "",
         'logDir': ""}

    config['GUSTO_Parameters'] = \
        {'# Full path for GParameterFile!': None,
         'GParameterFile': ''}

    configStr = '# GUSTO Level 2 Pipeline Configuration File \
                [DEFAULT] \n\
                number_of_cpus = 4 \n\
                 \n\
                [GUSTO_Directories] \n\
                configDir = "" \n\
                gustoIncomingDataDir = "" \n\
                gustoCIIDataDir = "" \n\
                gustoNIIDataDir = "" \n\
                gustoOIDataDir = "" \n\
                logDir = "" \n\
                 \n\
                [GUSTO_Parameters] \n\
                # with full path! \n\
                GParameterFile = '' \n\
                '
    with open(configFileName, "w") as f:
        config.write(f)
