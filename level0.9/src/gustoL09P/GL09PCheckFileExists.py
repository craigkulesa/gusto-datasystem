"""Function checks if GUSTO related config/parameter file exists


History:
Version 0.1 
"""
import os
from os.path import expanduser
import ntpath

__version__ = '0.1'
__date__ = '9/19/2024'
__updated__ = '20240919'
__author__ = 'V. Tolls, CfA | Harvard & Smithsonian'



def GUSTOCheckFileExists(ifile, noexceptionflag=False, returnboolean=True, 
                         replacementflag=False, verbose=False):
    """Function checking if ifile exists in the provided directroy, 
    or if there is a replacement in the class base directory or 
    the directory ~/.gusto.
    On Unix, the "~" for the home directory will be expanded to the full 
    path of the home directory.

    INPUT:
    ------
    ifile: string
        full path to file
    returnboolean: boolean
        returns True/False or full path
    replacementflag: boolean
        checking if also checking for replacement files in class base direcotry or 
        the GUSTO pipeline setup directory ~/.gusto

    OUTPUT:
    -------
        returns 0 or full path to file
    """

    if verbose:
        print(f'checking ifile: {ifile} \n{os.path.exists(ifile)}\n')

    # remove potential use of "~" for home directory
    ifile = os.path.expanduser(ifile)
    uhome = expanduser("~")
    head, tail = os.path.split(ifile)
    froot, ext = os.path.splitext(tail)

    if os.path.exists(ifile):
        # check standard configuration file exists
        if returnboolean:
            return True
        else:
            return ifile
    elif (os.path.exists('./' + tail)) & replacementflag:
        if returnboolean:
            return True
        else:
            return os.getcwd() + '/.gusto/' + tail
    elif (os.path.exists(uhome + '/.gusto/' + tail)) & replacementflag:
        if returnboolean:
            return True
        else:
            return uhome+'/.gusto/' + tail
    else:
        if noexceptionflag:
            print(f'Configuration file does not exist.\n'+
                            '  Provide a valid configuration file.\n'+
                            '  The standard configuration file should '+
                            'be: %s\n'%(uhome+'/.gusto/' + tail)+
                            '  or provide the path to a valid alternate location.')
        else:
            raise NoneException(f'Configuration file does not exist.\n'+
                                '  Provide a valid configuration file.\n'+
                                '  The standard configuration file should '+
                                'be: %s\n'%(uhome+'/.gusto/' + tail)+
                                '  or provide the path to a valid alternate location.')


class NoneException(Exception):
    pass


if __name__ == "__main__":

    test = 1

    if test==1:
        uhome = expanduser("~")
        
        # check if standard GUSTO Level 2 pipeline configuration file exists
        tfile = uhome+'/.gusto/GL2Pconfig.txt'
        res = GUSTOCheckFileExists(tfile, noexceptionflag=True, returnboolean=True, verbose=True)
        print(f'Checking if file {tfile} exist (boolean response): {res}')
        print()
        
        tfile = '~/.gusto/GL2Pconfig.txt'
        res = GUSTOCheckFileExists(tfile, noexceptionflag=True, returnboolean=False)
        print(f'Checking if file {tfile} exist: {res}')
        print()
        
        tfile = './GL2Pconfig.txt'
        res = GUSTOCheckFileExists(tfile, noexceptionflag=True)
        print(f'Checking if file {tfile} exist: {res}')
        print()
        
        tfile = '~/.gusto/SL09Pconfig.txt'
        res = GUSTOCheckFileExists(tfile, noexceptionflag=True)
        print(f'Checking if file {tfile} exist: {res}')
        print()

