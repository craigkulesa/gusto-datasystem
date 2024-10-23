import os
import sys
import glob
import shutil
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from pybaselines import Baseline, utils
from astropy import constants as const

def doStuff(scan):
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
        print("No history entries found in the header .. continuing")
    else:
        history_list = header['HISTORY']
        if f'known bad {line} channels flagged' in history_list:
            print('flagging already done .. stopping')
            exit

    # Known spikes
    # LO Interference
    # MHz       Channel
    # 317-342   31-36   66-70
    # 649-659   66-69   134-136
    # 972       100-101 201-203
    #
    # Iridium
    # 1610      165-168 331-335
    # Look through all of the spectra to set row_flags for noisy data
    if(line=='NII'):
        for i in range(nrow-1):
            CHANNEL_FLAG[i][31:36]   = 1    # LO 1 330 MHz
            CHANNEL_FLAG[i][66:69]   = 1    # LO 2 656 MHz
            CHANNEL_FLAG[i][100:101] = 1    # LO 3 980 MHz
            CHANNEL_FLAG[i][165:168] = 1    # Iridium 1 1616-1625 MHz

        hdu[0].header.add_history('known bad NII channels flagged')

    elif(line=='CII'):
        for i in range(nrow-1):
            CHANNEL_FLAG[i][66:70]   = 1    # LO 1 330 MHz
            CHANNEL_FLAG[i][134:136] = 1    # LO 2 656 MHz
            CHANNEL_FLAG[i][201:203] = 1    # LO 3 980 MHz
            CHANNEL_FLAG[i][331:335] = 1    # Iridium 1 1616-1625 MHz

        hdu[0].header.add_history('known bad CII channels flagged')

    else:
        hdu[0].header.add_history('error, no NII or CII line found')

    # Write the change back to the fits file
    hdu.flush()


directory = "/home/young/Projects/GUSTO/level0.8/ACS3"
partial = sys.argv[1]
search_files = sorted(glob.glob(f"{directory}/{partial}.fits"))
for file in (search_files):
    doStuff(file)


