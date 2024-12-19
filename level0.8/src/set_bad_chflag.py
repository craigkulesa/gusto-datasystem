import os
import sys
import glob
import shutil
import configparser
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from pybaselines import Baseline, utils
from astropy import constants as const

import argparse
import multiprocessing
from functools import partial

from tqdm import tqdm

def doStuff(scan, args):
    # open fits file
    hdu    = fits.open(scan, mode='update')
    # read header
    header = hdu[0].header
    npix   = header['NPIX']
    line   = header['LINE']
    # read data_table
    data   = hdu[1].data
    spec   = data['spec'] 
    nrow   = len(spec)
    CHANNEL_FLAG = data['CHANNEL_FLAG']

    if 'HISTORY' not in header:
        tqdm.write("{:s}, No HISTORY in header .. flagging bad channels".format(scan))
    else:
        history_list = header['HISTORY']

        if (f'known bad {line} channels flagged' not in history_list):
            tqdm.write("{:s}, No bad channel flags in header .. flagging bad channels".format(scan))

        elif (f'known bad {line} channels flagged' in history_list) and args.force:
            tqdm.write("{:s}, Flagging already done .. but forcing".format(scan))

        elif (f'known bad {line} channels flagged' in history_list) and not args.force:
            tqdm.write("{:s}, Flagging already done .. stopping".format(scan))
            return None

    # Known spikes
    # LO Interference
    # MHz       Channel
    # 317-342   31-36   66-70
    # 649-659   66-69   134-136
    # 972       100-101 201-203
    #
    # Iridium
    # 1314      133-136 266-272
    # 1458      147-152 294-303
    # 1610      165-168 331-335
    # Look through all of the spectra to set row_flags for noisy data

    #RESET
    for i in range(nrow-1):
        CHANNEL_FLAG[i] = 0

    if(line=='NII'):
        for i in range(nrow-1):
            CHANNEL_FLAG[i][31:36]   |= 1<<4    # LO 1 330 MHz Bit 4 SPUR
            CHANNEL_FLAG[i][66:69]   |= 1<<4    # LO 2 656 MHz Bit 4 SPUR
            CHANNEL_FLAG[i][100:101] |= 1<<4    # LO 3 980 MHz Bit 4 SPUR
            CHANNEL_FLAG[i][133:136] |= 1<<7    # Iridium 1314 MHz *Variable Bit 7 VARIABLE SPUR
            CHANNEL_FLAG[i][147:152] |= 1<<7    # Iridium 1458 MHz *Variable Bit 7 VARIABLE SPUR
            CHANNEL_FLAG[i][165:168] |= 1<<4    # Iridium 1 1616-1625 MHz Bit 4 SPUR

        hdu[0].header.add_history('known bad NII channels flagged')

    elif(line=='CII'):
        for i in range(nrow-1):
            CHANNEL_FLAG[i][66:70]   |= 1<<4    # LO 1 330 MHz Bit 4 SPUR
            CHANNEL_FLAG[i][134:136] |= 1<<4    # LO 2 656 MHz Bit 4 SPUR
            CHANNEL_FLAG[i][201:203] |= 1<<4    # LO 3 980 MHz Bit 4 SPUR
            CHANNEL_FLAG[i][266:272] |= 1<<7    # Iridium 1314 MHz *Variable Bit 7 VARIABLE SPUR
            CHANNEL_FLAG[i][294:303] |= 1<<7    # Iridium 1458 MHz *Variable Bit 7 VARIABLE SPUR
            CHANNEL_FLAG[i][331:335] |= 1<<4    # Iridium 1 1616-1625 MHz Bit 4 SPUR

        hdu[0].header.add_history('known bad CII channels flagged')

    else:
        hdu[0].header.add_history('error, no NII or CII line found')

    # Write the change back to the fits file
    hdu.flush()

    return None


if __name__ == "__main__":
    # ConfigParser Object
    config = configparser.ConfigParser()

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", help="\tFilename partial", default="ACS*")
    parser.add_argument("--force", help="\tForce update", action=argparse.BooleanOptionalAction, default=False)
    args = parser.parse_args()

    # Read config file for Data Paths
    config.read('../../common/config.ini')
    path = config.get('Paths', 'L08_path')

    partial_filename = args.file
    files = sorted(glob.glob(f"{path}/{partial_filename}.fits"))

    with multiprocessing.Pool() as pool:
        pool.map(partial(doStuff, args=args), files)
            




