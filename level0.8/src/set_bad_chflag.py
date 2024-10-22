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
    hdu    = fits.open(scan, mode='update')
    data   = hdu[1].data
    spec   = data['spec'] 
    CHANNEL_FLAG = data['CHANNEL_FLAG']
    npix   = hdu[0].header['NPIX']
    line   = hdu[0].header['LINE']
    nrow   = len(spec)
    spec   = data['spec']

    # Known spikes
    # LO Interference
    # MHz   Channel
    # 317-342   65, 66, 67, 68
    # 649-659   133, 134
    # 972       198, 199
    #
    # Iridium
    # 1610      329, 330
    # Look through all of the spectra to set row_flags for noisy data
    if(line=='NII'):
        for i in range(nrow-1):
            #for j in range(npix-1):
            #    CHANNEL_FLAG[i][j]   = 0
            CHANNEL_FLAG[i][31:36]   = 1    # LO 1 330 MHz
            CHANNEL_FLAG[i][66:69]   = 1    # LO 2 656 MHz
            CHANNEL_FLAG[i][100:101] = 1    # LO 3 980 MHz
            CHANNEL_FLAG[i][165:168] = 1    # Iridium 1 1616-1625 MHz
    elif(line=='CII'):
        for i in range(nrow-1):
            #for j in range(npix-1):
            #    CHANNEL_FLAG[i][j]   = 0
            CHANNEL_FLAG[i][66:70]   = 1    # LO 1 330 MHz
            CHANNEL_FLAG[i][134:136] = 1    # LO 2 656 MHz
            CHANNEL_FLAG[i][201:203] = 1    # LO 3 980 MHz
            CHANNEL_FLAG[i][331:335] = 1    # Iridium 1 1616-1625 MHz

    # Write the change back to the fits file
    print(scan)
    hdu.flush()


directory = "/home/young/Projects/GUSTO/level0.8/ACS5"
partial = sys.argv[1]
search_files = sorted(glob.glob(f"{directory}/{partial}.fits"))
for file in (search_files):
    doStuff(file)


