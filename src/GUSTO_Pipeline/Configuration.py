"""Command line argument parser for GUSTO pipeline
"""
import argparse
import configargparse
import logging
from importlib.resources import files

cfg_file0 = files('GUSTO_Pipeline') / 'config.gusto'

def getConfiguration(verbose=False):
    """Function parsing the command line input.
    """    
    parser = configargparse.ArgParser(default_config_files=[cfg_file0, '~/.config.gusto'], prog='runGUSTO', ignore_unknown_config_file_keys=True)
    parser.add('-c', '--config', required=False, is_config_file=True, help='config file path')
    parser.add('-e', '--erase', required=False, action=argparse.BooleanOptionalAction, help='erase contents of output folder before starting', default=False)
    parser.add('-v', '--verbose', required=False, action=argparse.BooleanOptionalAction, help='sets verbose output', default=False)
    parser.add('-d', '--debug', required=False, action=argparse.BooleanOptionalAction, help='sets debugging output', default=False)
    parser.add('-p', '--path', required=False, help="\tOverarching data path")
    parser.add('-b', '--band', required=True, help='GUSTO band range: 1, 2, or 1 2 for both', nargs="+")
    parser.add('-l', '--level', required=False, help='Pipeline step range to process: single number or start stop for multiple steps', nargs="+")
    parser.add('-j', '--cpus', required=False, help='set number of CPUs to use')
    parser.add('-s', '--scanid', required=True, help='scanID range', nargs=2)
    parser.add('--polyorder', required=False, help='Baseline polynomial order', default='1')
    parser.add('--calmethod', required=False, help='Calibration method to use', default='cal_scaledGainHOTs')
    parser.add('--despurmethod', required=False, help='despur method to use', default='polyRes')
    parser.add('--spurchannelfilter', required=False, action=argparse.BooleanOptionalAction, help='apply filter for spur masks', default=False)
    
    args = parser.parse_args()
    print(parser.format_values())

    if not args.path.endswith('/'):
        args.path=args.path+'/'
    return args
