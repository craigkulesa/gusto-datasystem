#! /usr/bin/env python
"""
Support program creating a GUSTO Level 2 Pipeline configuration file.
Newly create file must be modified to fit inidividual setup.
"""
import os
import sys
from os.path import expanduser, expandvars

__date__ = '9/19/2024'
__updated__ = '20240919'
__version__ = '0.1'
__author__ = 'V. Tolls, CfA | Harvard & Smithsonian'

uhome = expanduser("~")


def GUSTOCreateConfigFile(configFile='GL2Pconfig.txt', configPath=uhome + '/.gusto/', overwrite=False):
    """
    """

    if (configFile == None) | (configFile == ''):
        configFile = 'GL2Pconfig.txt'
    else:
        # check if full path is provided (not necessary!)
        head, tail = os.path.split(configFile)
        if head != '':
            # path or something provided
            if os.path.exists(head):
                # if the path exists, we set it as configPath
                configPath = head
            else:
                print('Error in provided configFile variable: ', configFile)
                print('Terminating ...')
                sys.exit(1)

    if not os.path.exists(configPath):
        print('Error: Directory of configuration file does not exist!')
        print('Do you want to create the directory?')
        inp = input('(yes/NO: ')
        if inp.lower() == 'yes':
            print('Trying to create directory')
            if os.makedirs(configPath):
                print('Directory successfully created!')
            else:
                print('Error: Could not create directory. Please check!')

    ofile = os.path.join(configPath, configFile)

    # last test: does file exist!
    if os.path.exists(ofile) & (not overwrite):
        print('Error: file already exists: ', ofile)
        print('If a new file is desired, either delete the existing file of set the overwrite keyword to True.')
        sys.exit(1)

    ostr = '# GUSTO Level 2 Pipeline Configuration File \n' + \
        '[DEFAULT] \n' + \
        'number_of_cpus = 4 \n' + \
        ' \n' + \
        '[GUSTO_Directories] \n' + \
        'configDir = ~/.gusto/GL2Pconfig.txt \n' + \
        'workDir = /rdat/Projects/GUSTO/DATA/GL2P/ \n' + \
        'gustoIncomingDataDir = /rdat/Projects/GUSTO/DATA/GL2P/Incoming/ \n' + \
        'gustoCIIDataDir = /rdat/Projects/GUSTO/DATA/GL2P/CIIdata \n' + \
        'gustoNIIDataDir = /rdat/Projects/GUSTO/DATA/GL2P/NIIdata \n' + \
        'gustoOIDataDir = /rdat/Projects/GUSTO/DATA/GL2P/OIdata/ \n' + \
        'logDir = /rdat/Projects/GUSTO/DATA/GL2P/log/ \n' + \
        ' \n' + \
        '[GUSTO_Parameters] \n' + \
        '# with full path! \n' + \
        'GParameterFile = ~/git/gustoL09/Data/GUSTO_Baseline_Data_draft3.xlsx'

    if overwrite:
        mode = 'w'
    else:
        mode = 'x'

    with open(ofile, mode, encoding='utf-8') as f:
        f.write(ostr)

    print('Successfully created a template GUSTO L2 configuration File:\n\n  %s' % (ofile))
    print('\nNow customize the file with your preferred text editor.')


if __name__ == '__main__':

    # regular mode is 0 !!!
    test = 0

    if test == 0:
        GUSTOCreateConfigFile()

    elif (test == 1):
        overwrite = True
        GUSTOCreateConfigFile(overwrite=overwrite)

