#!/usr/bin/env python3.11

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

sys.path.append("../../common/")
import flagdefs as flags

commit_info = ''


def clear_folder(folder_path):
    print('Erasing contents of '+folder_path)
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")

            
def flatten(xss):
    return [x for xs in xss for x in xs]


def runGitLog():
    try:
        result = subprocess.run(['git', 'log', '-1', '--format=%cd', '--date=format-local:%Y-%m-%d %H:%M:%S %Z', '--pretty=format:Level 0.8 commit %h by %an %ad', '--', './L08-pipeline.py'], capture_output=True, text=True, check=True)
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

    
def makeFileGlob(path, options):
    fileList = []
    if int(options.scanid[1]) - int(options.scanid[0]) > 28422: # do them all
        fileList.append(glob.glob(path+'/*.fits'))
    else:
        for scanID in range(int(options.scanid[0]), int(options.scanid[1])+1):
            fileList.append(glob.glob(path+'/*_*_'+str(scanID).zfill(5)+'.fits'))
    return(flatten(fileList))


def doStuff(scan, options):
    # open fits file, header, then data
    if options.update == True:
        hdu    = fits.open(scan, mode='update')
    else:
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

        if (f'known suspect {line} channels flagged' in history_list):
            if options.force:
                print("{:s}, Flagging already done .. but forcing".format(scan))

            else:
                print("{:s}, Flagging already done .. stopping".format(scan))
                return None

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
            CHANNEL_FLAG[i][31:36]   |= flags.ChanFlags.VARIABLE_SPUR    # LO 1 330 MHz 
            CHANNEL_FLAG[i][66:69]   |= flags.ChanFlags.VARIABLE_SPUR    # LO 2 656 MHz 
            CHANNEL_FLAG[i][100:101] |= flags.ChanFlags.VARIABLE_SPUR    # LO 3 980 MHz 
            CHANNEL_FLAG[i][133:136] |= flags.ChanFlags.VARIABLE_SPUR    # Iridium 1314 MHz 
            CHANNEL_FLAG[i][147:152] |= flags.ChanFlags.VARIABLE_SPUR    # Iridium 1458 MHz 
            CHANNEL_FLAG[i][165:168] |= flags.ChanFlags.VARIABLE_SPUR    # Iridium 1 1616-1625 MHz Bit 
            CHANNEL_FLAG[i][0:31]    |= flags.ChanFlags.OOB    # Out of band 0-300MHz (lower)
            CHANNEL_FLAG[i][410:511] |= flags.ChanFlags.OOB    # Out of band 4000-5000MHz (upper)

        hdu[0].header.add_history('known bad NII channels flagged')

    elif(line=='CII'):
        x0=25
        x1=50
        l0=25
        mixers=[2,3,5,8]
        thresh=[.01, .1, .1, .0004]
        
        for i in range(nrow):
            CHANNEL_FLAG[i][66:70]   |= flags.ChanFlags.VARIABLE_SPUR    # LO 1 330 MHz Bit 
            CHANNEL_FLAG[i][134:136] |= flags.ChanFlags.VARIABLE_SPUR    # LO 2 656 MHz Bit 
            CHANNEL_FLAG[i][201:203] |= flags.ChanFlags.VARIABLE_SPUR    # LO 3 980 MHz Bit 
            CHANNEL_FLAG[i][266:272] |= flags.ChanFlags.VARIABLE_SPUR    # Iridium 1314 MHz 
            CHANNEL_FLAG[i][294:303] |= flags.ChanFlags.VARIABLE_SPUR    # Iridium 1458 MHz 
            CHANNEL_FLAG[i][331:335] |= flags.ChanFlags.VARIABLE_SPUR    # Iridium 1 1616-1625 MHz 
            CHANNEL_FLAG[i][0:62]    |= flags.ChanFlags.OOB    # Out of band 0-300MHz (lower)
            CHANNEL_FLAG[i][820:1023]|= flags.ChanFlags.OOB    # Out of band 4000-5000MHz (upper)

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
            ROW_FLAG[i] |=  (flags.RowFlags.RINGING_BIT1 | flags.RowFlags.RINGING_BIT0)  # set ringing bits
        else:
            ROW_FLAG[i] &= ~(flags.RowFlags.RINGING_BIT1 | flags.RowFlags.RINGING_BIT0)  # clear ringing bits
      
    hdu[0].header.add_history('badly fringing rows flagged')
    header['DLEVEL'] = 0.8
    now = datetime.now()
    header['PROCTIME'] = now.strftime("%Y%m%d_%H%M%S")
    header['COMMENT'] = commit_info
    
    # Write the changes back to fits in update mode or as new file
    if(options.update == True):
        hdu.flush()
    else:
        hdu.writeto(scan.replace("level0.7", "level0.8"))
    return None


if __name__ == "__main__":
    parser = configargparse.ArgParser(default_config_files=['../../config.gusto', '~/.config.gusto'], ignore_unknown_config_file_keys=True)
    parser.add('-c', '--config', required=False, is_config_file=True, help='config file path')
    parser.add('-e', '--erase', required=False, action=argparse.BooleanOptionalAction, help='erase contents of output folder before starting')
    parser.add('-p', '--path', required=False, help="\tData path")
    parser.add('-u', '--update', required=False, help="\tUpdate fits file in place, vs making a new file", action=argparse.BooleanOptionalAction, default=False)
    parser.add('-j', '--cpus', required=False, help='set number of CPUs to use')
    parser.add('-f', '--force', required=False, help="\tForce update even if step appears to have been done already", action=argparse.BooleanOptionalAction, default=False)
    parser.add('-s', '--scanid', required=True, help='scanID range', nargs=2)

    options = parser.parse_args()

    print(options)
    print(parser.format_values())

    commit_info = runGitLog()  # lookup git commit info only once
    
    dirDataOut = options.path+'level0.8/'
    if not os.path.exists(dirDataOut):
        os.makedirs(dirDataOut)
        print(f"Directory {dirDataOut} created.")
    else:
        print(f"Reusing directory {dirDataOut}")
        if options.erase:
            clear_folder(dirDataOut)

    if options.update == True:  # update within level 0.8 directory
        path = dirDataOut
        print("Using level 0.8 folder for input and updating files in place")
    else: # read from level 0.7 directory
        path = options.path+'level0.7'
        print("Using level 0.7 folder for input and creating new level 0.8 files")
    files = makeFileGlob(path, options)
    print(f"Found {len(files)} files to process")
    
    if(options.cpus):
        pool = Pool(processes=int(options.cpus))
    else:
        pool = Pool()
    pool.map(partial(doStuff, options=options), files)
            
