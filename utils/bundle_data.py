#!/usr/bin/env python3.11
import os
import glob
import shutil
import configparser
import argparse
from astropy.io import fits


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
    config = configparser.ConfigParser()
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", help="\tFilename partial", default="*")
    args = parser.parse_args()

    # Read config file for Data Paths
    config.read('../common/config.ini')
    path = config.get('Paths', 'L10_path')
    outpath = config.get('Paths', 'bundle_path')
    partial_filename = args.file
    files = sorted(glob.glob(f"{path}/{partial_filename}.fits"))

    for file in files:
        print('file = ',file)
        doCopy(file, outpath)
