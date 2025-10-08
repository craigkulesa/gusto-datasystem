import os
import sys
import glob
from datetime import datetime
import shutil
import configargparse
import argparse
from astropy.io import fits
import numpy as np
from multiprocessing import Pool
from functools import partial
import subprocess
from .flagdefs import *
from .DataIO import *


def runGitLog():
    try:
        result = subprocess.run(['git', 'log', '-1', '--format=%cd', '--date=format-local:%Y-%m-%d %H:%M:%S %Z', '--pretty=format:Level 0.8 commit %h by %an %ad', '--', 'L08_chanflags.py'], capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        return f"Error: {e}"


def find_mixer(line, mixer):
    if(line=='NII'):
        mixers=[2,3,4,6]
    elif(line=='CII'):
        mixers=[2,3,5,8]

    # return index in mixer list of given mixer number
    if mixer in mixers:
        return mixers.index(mixer)
    else:
        return 5


def L08_Pipeline(args, scanRange):
    global commit_info
    dirDataOut = args.path+'level0.8/'
    prefix = ['NII_', 'CII_'] 
    os.makedirs(dirDataOut, exist_ok=True)
    if args.erase:
        clear_folder(dirDataOut)

    commit_info = runGitLog()  # do this only once     
    path = args.path+'level0.7/'
    sum_files = 0
    
    for band in args.band:
        files = makeFileGlob(path, prefix[int(band)-1], 'fits', scanRange)
        sum_files += len(files)
        if(args.cpus):
            pool = Pool(processes=int(args.cpus))
        else:
            pool = Pool()
        pool.map(partial(doStuff, options=args), files)

    return sum_files
    

def doStuff(scan, options):
    global commit_info
    # open fits file, header, then data
    hdu    = fits.open(scan)
    header = hdu[0].header
    npix   = header['NPIX']
    line   = header['LINE']
    data   = hdu[1].data
    spec   = data['DATA']
    mixer  = data['MIXER']
    nrow   = len(spec)
    ROW_FLAG  = data['ROW_FLAG']
    CHANNEL_FLAG = data['CHANNEL_FLAG']

    if 'HISTORY' in header:
        history_list = header['HISTORY']

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

    # RESET CHANNEL FLAGS -- no prior processing has set any of these yet
    for i in range(nrow):
        CHANNEL_FLAG[i] = 0

    if(line=='NII'):
        x0=12
        x1=24
        l0=12
        mixers=[2,3,4,6]
        thresh=[.1, .1, .1, .1]
        
        for i in range(nrow):
            CHANNEL_FLAG[i][31:36]   |= ChanFlags.VARIABLE_SPUR    # LO 1 330 MHz 
            CHANNEL_FLAG[i][66:69]   |= ChanFlags.VARIABLE_SPUR    # LO 2 656 MHz 
            CHANNEL_FLAG[i][100:101] |= ChanFlags.VARIABLE_SPUR    # LO 3 980 MHz 
            CHANNEL_FLAG[i][133:136] |= ChanFlags.VARIABLE_SPUR    # Iridium 1314 MHz 
            CHANNEL_FLAG[i][147:152] |= ChanFlags.VARIABLE_SPUR    # Iridium 1458 MHz 
            CHANNEL_FLAG[i][165:168] |= ChanFlags.VARIABLE_SPUR    # Iridium 1 1616-1625 MHz Bit 
            CHANNEL_FLAG[i][0:31]    |= ChanFlags.OOB    # Out of band 0-300MHz (lower)
            CHANNEL_FLAG[i][410:511] |= ChanFlags.OOB    # Out of band 4000-5000MHz (upper)

        hdu[0].header.add_history('known bad NII channels flagged')

    elif(line=='CII'):
        x0=25
        x1=50
        l0=25
        mixers=[2,3,5,8]
        thresh=[.01, .1, .1, .0004]
        
        for i in range(nrow):
            CHANNEL_FLAG[i][66:70]   |= ChanFlags.VARIABLE_SPUR    # LO 1 330 MHz Bit 
            CHANNEL_FLAG[i][134:136] |= ChanFlags.VARIABLE_SPUR    # LO 2 656 MHz Bit 
            CHANNEL_FLAG[i][201:203] |= ChanFlags.VARIABLE_SPUR    # LO 3 980 MHz Bit 
            CHANNEL_FLAG[i][266:272] |= ChanFlags.VARIABLE_SPUR    # Iridium 1314 MHz 
            CHANNEL_FLAG[i][294:303] |= ChanFlags.VARIABLE_SPUR    # Iridium 1458 MHz 
            CHANNEL_FLAG[i][331:335] |= ChanFlags.VARIABLE_SPUR    # Iridium 1 1616-1625 MHz 
            CHANNEL_FLAG[i][0:62]    |= ChanFlags.OOB    # Out of band 0-300MHz (lower)
            CHANNEL_FLAG[i][820:1023]|= ChanFlags.OOB    # Out of band 4000-5000MHz (upper)

        hdu[0].header.add_history('known suspect CII channels flagged')

    else:
        hdu[0].header.add_history('error, no NII or CII line found')

    # compute fringing rowflag
    xdata = np.arange(0,l0)
    data  = np.zeros(l0)
    for i in range(nrow):
        # compute standard deviation
        ydata = spec[i][x0:x1]
        z = np.polyfit(xdata, ydata, 5)
        p = np.poly1d(z)
        for j in range(l0):
            data[j] = ydata[j] - p(j)

        if(np.std(data) > thresh[find_mixer(line, mixer[i])]):
            ROW_FLAG[i] |=  (RowFlags.RINGING_BIT1 | RowFlags.RINGING_BIT0)  # set ringing bits
        else:
            ROW_FLAG[i] &= ~(RowFlags.RINGING_BIT1 | RowFlags.RINGING_BIT0)  # clear ringing bits
      
    hdu[0].header.add_history('badly fringing rows flagged')
    header['DLEVEL'] = 0.8
    now = datetime.now()
    header['PROCTIME'] = now.strftime("%Y%m%d_%H%M%S")
    header['COMMENT'] = commit_info
    
    # Write the changes back to fits as new file
    hdu.writeto(scan.replace("level0.7", "level0.8").replace("L07", "L08"), overwrite=True)
    return None


