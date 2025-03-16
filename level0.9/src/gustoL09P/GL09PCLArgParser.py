"""Command line argument parser for GUSTO Level 0.9 pipeline
"""
import argparse
import textwrap
import importlib
import logging
from pprint import pprint
from .GL09PUtils import *

cfg_file0 = importlib.resources.files('gustoL09P') / 'Data/GL09P_setup_pars.txt'


def GL09PCLArgParser(verbose=False):
    r"""Function parsing the command line input.

    Parameters
    ----------
    verbose : bool
        Print more information when processing.

    Returns
    -------
    args : dict
        Command line arguments

    Examples
    --------
    """
    
    
    parser = argparse.ArgumentParser(
        prog='execGL09P',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
            GUSTO Level 1 Pipeline CLI
            --------------------------------
                I have indented it
                exactly the way
                I want it
            '''),
        epilog='')
    parser.add_argument('--configFile', '-c', nargs='?', type=str,
                        default=cfg_file0,
                        help='GUSTO L1 Pipeline configuration file. This file contains the data directories, etc.')
    parser.add_argument('--startLevel', '-s', nargs='?', type=str,
                        default='0.8',
                        help='GUSTO Data Level the pipeline should start processing (default=0.8); \npossible entries are 0.8, 0.9 and 0.95')
    parser.add_argument('--lines', nargs='?', type=str,
                        default='CII',
                        help='Line to be processed: either CII or NII, if both lines are desired, do not use this option but the config file')
    parser.add_argument('--endLevel', '-e', nargs='?', type=str,
                        default='0.9',
                        help='GUSTO Data Level produced by pipeline (default=1.0); possible entries are 0.9 and 1.0')
    parser.add_argument('--scanRange', '-r', nargs=2, type=int,
                        default=[2000, 30000],
                        help='Range of scans to be processed by pipeline.')
    parser.add_argument('--loglevel', '-l', type=str,
                        #default='INFO',
                        help='sets the log level of the {tpipe}')
    parser.add_argument('--calmethod', '-m', nargs='?', type=int,
                        #default='2',
                        help='select GUSTO spectrum calibration method')
    parser.add_argument('--verbose', '-v', action=argparse.BooleanOptionalAction,
                        help='sets verbosity of the {tpipe}')
    parser.add_argument('--debug', '-d', action=argparse.BooleanOptionalAction,
                        help='sets settings for debugging pipeline')
    args = parser.parse_args()
    if verbose:
        print('commandline arguments:\n', args, '\n')
        print('configFile: ', args.configFile)
        print()

    # this overrides the verbosity from above in case it is enabled
    if args.verbose is not None:
        verbose = args.verbose

    return args



def manageArgs(cfi, args, verbose=False):
    r"""Function applying the command line arguments to the configuration record array.

    Parameters
    ----------
    cfg : dict
        Configuration record array.
    args : dict
        Command line arguments.

    Returns
    -------
    cfg : dict
        Updated configuration record array.

    Examples
    --------
    """


    if args.debug is not None:
        print('args.debug: ', args.debug, type(args.debug))
        cfi['gprocs']['debug'] = args.debug

    if args.lines is not None:
        print('args.lines: ', args.lines)
        cfi['gprocs']['lines'] = [args.lines]
    
    if args.calmethod is not None:
        print('args.calmethod: ', args.calmethod, type(args.calmethod), cfi['gprocs']['drmethod'], type(cfi['gprocs']['drmethod']))
        cfi['gprocs']['drmethod'] = args.calmethod
    else:
        cfi['gprocs']['drmethod'] = int(cfi['gprocs']['drmethod'])
        print('calmethod: ', cfi['gprocs']['drmethod'], type(cfi['gprocs']['drmethod']))

    if args.loglevel is not None:
        numeric_level = getattr(logging, args.loglevel.upper(), None)
    elif cfi['gprocs']['loglevel'] != '':
        print('cfi: ', cfi['gprocs']['loglevel'].upper())
        numeric_level = getattr(logging, cfi['gprocs']['loglevel'].upper(), None)
    else:
        numeric_level = getattr(logging, 'INFO', None)
    cfi['gprocs']['loglevel'] = numeric_level
    
    # check the rowflagfilter: numeric or string
    cfi['gprocs']['rowflagfilter'] = anaFlagString(cfi['gprocs']['rowflagfilter'])
    
    return cfi

