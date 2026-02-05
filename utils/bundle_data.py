#!/usr/bin/env python3.11
import os
import glob
import shutil
import argparse
import configargparse
from importlib.resources import files
from astropy.io import fits

cfg_file0 = files('GUSTO_Pipeline') / 'config.gusto'


def doCopy(scan, path):
    with fits.open(scan) as hdu:
        header = hdu[0].header    
        try:
            object = header['OBJECT']
            line   = header['LINE']
        except Exception as e:
            print("BAD header .. skipping file", scan)
            return
    dirOut = path + '/' + line + '/' + object + '/'
    if not os.path.exists(dirOut):
        os.makedirs(dirOut)
        print(f"Directory {dirOut} created.")
    shutil.copy(scan, dirOut)    
    return None


if __name__ == "__main__":
    parser = configargparse.ArgParser(default_config_files=[cfg_file0, '~/.config.gusto'], prog='runGUSTO', ignore_unknown_config_file_keys=True)
    parser.add('-p', '--path', required=False, help="\tOverarching data path")
    args = parser.parse_args()
    print(parser.format_values())

    if not args.path.endswith('/'):
        args.path=args.path+'/'
        
    # Read config file for Data Paths
    files = sorted(glob.glob(f"{args.path}level1/*.fits"))

    for file in files:
        print('file = ',file)
        doCopy(file, args.path + 'bundles/')
