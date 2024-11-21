"""Gusto Level 2 Pipeline class
"""

__version__ = '0.1'
__date__ = '20240919'
__updated__ = '20240919'
__author__ = 'V. Tolls, CfA | Harvard & Smithsonian'


import time
import configparser
import os
import sys
from os.path import expanduser, expandvars
from datetime import datetime
import pkg_resources
from pprint import pprint
import numpy as np
#from joblib import Parallel, delayed, cpu_count
import logging
log = logging.getLogger(__name__)

ProcFile0 = pkg_resources.resource_filename(
    'gustoL09P', 'data/CIIconfig.txt')

spectralLines = ['CII', 'NII', 'OI']


class GL09PipelineSetupClass:
    """GUSTO Level 1 Pipeline setup: reading configuration and processing data
    from config file.
    """

    def __init__(self, configFile=None, scriptFile=None, verbose=False):
        """Initializing the Level 1 pipeline setup


        Parameters
        ----------
        configFile : str
            Optional pipeline configuration file including data directories, etc.
        scriptFile: str
            the script file contains parameters to control the processing of
            the 2 spectral lines, i.e., line specific parameters
        """

        # read the configuration file
        #status = self.initializePipeline(verbose=verbose, configFile=configFile)
        self.dtstr = datetime.now().strftime("%Y%m%d_%H%M%S")

    def getProcessInfo(self, procFile=None, verbose=False):
        """Function reading the pipeline processing information file, e.g. ~/.gusto/CIIconfig.txt.
        The file should contain only processing information needed to process the selected 
        molecular line properly.
        
        The process information is stored in dictionaries for each of the processing parts. This 
        enables adding more processing parts later without changing the existing ones.
        
        The information is stored in the file pointed at by procFile. This is a user provide file.
        An example file is given for CII, located in the package Data directory: CIIconfig.txt
        
        Example accessing the data:
            gp2 = GL09PipelineClass(verbose=True)
            
            pctl = gp2.getProcessInfo(procFile='/Users/volkertolls/.gusto/CIIconfig.txt')
            
            keys = pctl.keys()
            for key in keys:
                print('%s: %s'%(key, pctl[key]))
            # print(pctl)
            # print(pctl.keys())
            print()
            print(pctl['Pars'].keys())
            print(pctl['PDefs']['line'])
        
        or from within the pipeline:
            print(self.GL09Pproc['Pars'])
        
        Parameters
        ----------
        procFile: str
                    Full path to the processing configuration file.

        """
        # if procFile:
        #     procPath = "/".join(procFile.split('/')[:-1])
        # else:
        #     procFile = ProcFile0
        #     print('Using built-in procFile: ', procFile)

        # if not os.path.isdir(procPath):
        #     print('Error: GUSTO configuration directory does not exist: %s \n'%(configPath))
        #     print('creating directory!\n')
        #     os.makedirs(procPath, exist_ok = True)

        # if not os.path.exists(procFile):
        #     print('Error: processing config file does not exist: ', procFile, '\n')
        #     sys.exit('Error: No processing control file provided.')

        pconfig = configparser.ConfigParser(interpolation=EnvInterpolation())
        pconfig.optionxform = str

        # read the configuration file
        pconfig.read(procFile)
        self.procFile = procFile

        pdict = {}
        psecs = pconfig.sections()

        # get all defaults
        defaults = pconfig.defaults()
        pcdict = {}
        for key in defaults.keys():
            pcdict[key] = defaults[key]

        pipepars = pconfig["Parameters"]
        ppdict = {}
        for key in pipepars.keys():
            ppdict[key] = pipepars[key]

        pipectl = pconfig["PipelineProcs"]
        ctldict = {}
        for key in pipectl.keys():
            ctldict[key] = pipectl[key]

        pipearea = pconfig['InputData']
        regdict = {}
        for key in pipearea.keys():
            regdict[key] = pipearea[key]

        self.GL09Pproc = {}
        self.GL09Pproc['PDefs'] = pcdict
        self.GL09Pproc['Pars'] = ppdict
        self.GL09Pproc['PProcs'] = ctldict

        if verbose:
            print('getProcessInfo: procFile: ', self.procFile)
            print(self.GL09Pproc)
            print()

        return self.GL09Pproc

    def initializePipeline(self, verbose=False, configFile=None):
        """Function reading the pipeline configuration file and 
        setting up the pipeline. The pipeline configuration file is different from the 
        process configuration. The configuration file contains only setup information
        like directories, number of cores to use, etc.

        Parameters
        ----------
        configfile : str
                    optional line-specific configuration file

        Returns
        -------
        line configuration structure
        """
        print('Executing managePipeline function')
        # read configuration file

        if configFile == None:
            uhome = expanduser("~")
            self.configFile = uhome+'/.gusto/GL09P_setup_pars.txt'
        else:
            self.configFile = configFile

        self.configFileExists()

        if verbose:
            print('Reading pipeline configuration file: %s\n' %
                  (self.configFile))

        config = configparser.ConfigParser(interpolation=EnvInterpolation())
        config.optionxform = str

        # read the configuration file
        config.read(self.configFile)

        cfg_sect = config.sections()

        sdict = {}

        # get all defaults
        defaults = config.defaults()
        defdict = {}
        for key in defaults.keys():
            defdict[key] = defaults[key]

        # GUSTO_Directories
        # gpars = config.options("GUSTO_Parameters")
        gdirs = config["GUSTO_Directories"]
        gddict = {}
        for key in gdirs.keys():
            if 'work' in key:
                # add the date/time stamp to the working directory
                gddict[key] = gdirs[key]+'_'+self.dtstr
            else:
                gddict[key] = gdirs[key]

        # GUSTO_Parameters
        # gpars = config.options("GUSTO_Parameters")
        gpars = config["GUSTO_Parameters"]
        gpdict = {}
        for key in gpars.keys():
            gpdict[key] = gpars[key]

        # GUSTO_Processing
        gpars = config["GUSTO_Processing"]
        gcdict = {}
        for key in gpars.keys():
            gcdict[key] = gpars[key]

        self.GL09Pconf = {}
        self.GL09Pconf['default'] = defdict
        self.GL09Pconf['gdirs'] = gddict
        self.GL09Pconf['gparams'] = gpdict
        self.GL09Pconf['gprocs'] = gcdict

        if verbose:
            print('loaded configuration file: %s\n' % (self.configFile))
            print('self.GL09Pconf: ', self.GL09Pconf, '\n')

        # checking if data directories are in place
        status = self.GL09PDirsSetup(verbose=verbose)

        # read master file log
        #status = self.GL09PreadMasterDataFileLog()

        return 1

    def getConfigInfo(self, verbose=False):
        # read the configuration file, but parsing is done in another function
        if verbose:
            print('getConfigInfo: ', self.GL09Pconf)

        return self.GL09Pconf

    def configFileExists(self):
        configPath = "/".join(self.configFile.split('/')[:-1])

        if not os.path.isdir(configPath):
            if self.configFile is not None:
                # check if the file is in the module data directory
                self.configFile = pkg_resources.resource_filename(
                        'gustoL09P', os.path.join('Data',self.configFile))
            else:
                print(
                    'Error: GUSTO configuration directory does not exist: %s \n' % (configPath))
                # print('creating directory!\n')
                # os.makedirs(configPath, exist_ok=True)
                sys.exit()

        if not os.path.exists(self.configFile):
            print('Error: config file does not exist: ', self.configFile)
            print('Creating a template file.')
            print('Directories must be filled in!\n')

            ostr = '# GUSTO Level 2 Pipeline Configuration File \n' + \
                '[DEFAULT] \n' + \
                '# The home directory can be either ${ghome} or be provided \n' + \
                '# with an explicit path. \n' + \
                '# home = ${ghome}/GUSTO/DATA/Pipeline \n' + \
                'ghome = /rdat/Projects/GUSTO/Data/Pipeline \n' + \
                ' \n' + \
                '[GUSTO_Processing] \n' + \
                'number_of_cpus = 4 \n' + \
                ' \n' + \
                '[GUSTO_Directories] \n' + \
                '# The directories should be relative to the home directory. \n' + \
                'workDir = ${ghome}/work \n' + \
                'configDir = ${ghome}/.gusto \n' + \
                'gustoIncomingDataDir = ${ghome}/incoming \n' + \
                'gustoCIIDataDir = ${ghome}/CIIData \n' + \
                'gustoNIIDataDir = ${ghome}/NIIData \n' + \
                'gustoOIDataDir = ${ghome}/OIData \n' + \
                'logDir = ${ghome}/work/logs \n' + \
                'tmpDir = ${ghome}/work/tmp \n' + \
                ' \n' + \
                '[GUSTO_Parameters] \n' + \
                '# with full path! \n' + \
                'GParameterFile = empty \n'

            with open(self.configFile, 'x', encoding='utf-8') as f:
                f.write(ostr)

    def GL09PreadMasterDataFileLog(self, verbose=False):
        if verbose:
            print('Reading the master data file log.\n')

        ifile = os.path.join(self.configDir, "GUSTO_masterDataFileLog.xls")

    def GL09PDirsSetup(self, verbose=False):
        if verbose:
            print('Checking if data and pipeline work directories exist.\n')

        # get the directories from config info
        gl4pDirs = self.GL09P_getDirs()
        #print(gl4pDirs)


        # 'workDir', 'L08DataDir', 'L09DataDir', 'L10DataDir', 'logDir', 'tmpDir'
        self.workDir = gl4pDirs['workDir']
        if not os.path.isdir(self.workDir):
            print('Creating working directory: %s \n' % (self.workDir))
            os.makedirs(self.workDir, exist_ok=True)
        self.logDir = gl4pDirs['logDir']
        if not os.path.isdir(self.logDir):
            print('Creating log directory: %s \n' % (self.logDir))
            os.makedirs(self.logDir, exist_ok=True)
        self.tmpDir = gl4pDirs['tmpDir']
        if not os.path.isdir(self.tmpDir):
            print('Error: tmp directory does not exist: %s \n' % (self.tmpDir))
            os.makedirs(self.tmpDir, exist_ok=True)

        # the next set of directories should exist!
        self.L08Dir = gl4pDirs['L08DataDir']
        if not os.path.isdir(self.L08Dir):
            print('Error: Level 0.8 data directory does not exist. \nCreating Level 0.8 data directory: %s \n' %
                  (self.L08Dir))
        self.L09Dir = gl4pDirs['L09DataDir']
        if not os.path.isdir(self.L09Dir):
            print('Error: Level 0.9 data directory does not exist. \nCreating Level 0.9 data directory: %s \n' %
                  (self.L09Dir))
        self.L10Dir = gl4pDirs['L10DataDir']
        if not os.path.isdir(self.L10Dir):
            print('Warning: Level 1.0 data directory does not exist. \nCreating Level 1.0 data directory: %s \n' %
                  (self.L10Dir))

        if verbose:
            print('work dir: ', self.workDir)
            print('tmp dir: ', self.tmpDir)
            print('log dir: ', self.logDir)
            print('Level 0.8 Data Dir: ', self.L08Dir)
            print('Level 0.9 Data Dir: ', self.L09Dir)
            print('Level 1.0 Data Dir: ', self.L10Dir)
            print()

        return 1

    def getGL09PConfigInfo(self):
        return self.GL09Pconf

    def readLineConfiguration(self, line, lineconfigfile=None, verbose=False):
        """Function to read the configuration file for the selected line.

        The files for the line configurations should be located in the 
        directory .gusto in the user's home directory. However, the function 
        can also be called using the lineconfigfile keyword in case the 
        configuration files are stroed elsewhere. The appropriate file names
        are CIIconfig.txt, NIIconfig.txt, and OIconfig.txt.

        Parameters
        ----------
        line : str
            selected emission line for processing: CII, NII, or OI
        lineconfigfile : str
            optional full path to line-specific configuration file to be used

        Returns
        -------
        line configuration structure
        """

        if line is None:
            print('Error: Please select one emission line: CII, NII, or OI')
            print('Aborting.')

        if line not in spectralLines:
            print('Error: valid lines are only: CII, NII, or OI!')
            print('Aborting.')

        line = line.upper()

        if lineconfigfile is None:
            if line == 'CII':
                lcfg = readLCfile('CIIconfig.txt', verbose=verbose)
            elif line == 'NII':
                lcfg = readLCfile('NIIconfig.txt', verbose=verbose)
            elif line == 'OI':
                lcfg = readLCfile('OIconfig.txt', verbose=verbose)

    def readLCfile(self, lcfFile, verbose=False):
        """Function reading the line configuration file.
        """
        lcp = configparser.ConfigParser()
        lcp.optionxform = str

        # read the configuration file
        lcp.read(lcfFile)
        self.lcfFile = lcfFile

        cfg_sect = config.sections()

        sdict = {}

        # get all defaults
        defaults = lcp.defaults()
        defdict = {}
        for key in defaults.keys():
            defdict[key] = defaults[key]

        # GUSTO_Directories
        # gpars = config.options("GUSTO_Parameters")
        lpars = lcp["GUSTO_Parameters"]
        lpdict = {}
        for key in lpars.keys():
            lpdict[key] = lpars[key]

        self.lineConf = {}
        self.lineConf['default'] = defdict
        self.lineConf['params'] = lpdict

        if verbose:
            print('loaded configuration file: ', self.lcfFile)
            print(self.lineConf)

    def getLineConfigInfo(self):
        return self.lineConf

    def GL09P_getDirs(self, verbose=False):
        return self.GL09Pconf['gdirs']

    def GL09P_loadData(self, verbose=False):
        """function loading the GUSTO Level 1 data
        """
        print('loadData module')

    def getGL09PVersion(self):
        """Function returning the GUSTO Level 2 pipeline version.
        """
        return __version__

    def getGL09PDate(self):
        """Function returning the GUSTO Level 2 pipeline version.

        Returns
        -------
        creation date and date of last update
        """
        return __date__, __updated__

    def analyzeSpectrumMask(self):
        """Function analyzing the full spectrum mask.
        It returns just a boolean 0 or 1 for not usable or usable the 
        spectrum, respectively.
        """
        print('analyzeSpectrumMask: to be implemented')
        return 1

    def analyzePixelMask(self):
        """Function analyzing flags set for each pixel in a spectrum.
        """
        print('analyzeSpectrumMask: to be implemented')
        return 1


class EnvInterpolation(configparser.ExtendedInterpolation):
    """Interpolation which expands environment variables in values."""

    def before_read(self, parser, section, option, value):
        value = super().before_read(parser, section, option, value)
        return expandvars(value)



def getFloatRange(crange):
    """Function converting string range to float range.
    
   
    Parameters
    ----------
    crange : string
            string version of the range, two numbers separated by comma
    
    
    Return
    ______
    numpy array of 2 floating point number for range.
    
    """
    cra = crange.split(',')
    range = [float(cra[0]), float(cra[1])]
    return np.array(range)


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from GL09PipelineClass import *
    from GL09PParameters import *

    test = 3

    if test == 0:
        print('Test run of GUSTO Level 2 Pipeline')
        gp1 = GL09PipelineClass(verbose=True)
        print()
        print(gp1.getGL09PDate())
        print(gp1.getGL09PVersion())

    if test == 1:
        print('Test run of GUSTO Level 2 Pipeline')
        gp1 = GL09PipelineClass(verbose=True)
        print()
        print(gp1.getGL09PConfigInfo())

    if test == 2:
        print('Test run of GUSTO Level 2 Pipeline')
        gp1 = GL09PipelineClass(verbose=True)

        gl4pDirs = gp1.GL09P_getDirs(verbose=True)

        print()
        print('GL09P directories:')
        print(gl4pDirs.keys())
        print()
        print(gl4pDirs['configDir'])

    if test == 3:
        print('Loading the process control file')
        gp1 = GL09PipelineClass(verbose=True)
        
        #pctl = gp2.getProcessInfo(procFile='/Users/volkertolls/.gusto/CIIconfig.txt')
        pctl = gp1.getProcessInfo(procFile='')
        
        keys = pctl.keys()
        for key in keys:
            print('%s: %s'%(key, pctl[key]))
        # print(pctl)
        # print(pctl.keys())
        print()
        print(pctl['Pars'].keys())
        print(pctl['PDefs']['line'])
        

