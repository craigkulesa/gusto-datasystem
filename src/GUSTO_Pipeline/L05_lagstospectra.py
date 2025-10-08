#!/usr/bin/env python3.11

import os
#import sys
import glob
import subprocess
from importlib.resources import files
from .Logger import *

logger = logging.getLogger('pipelineLogger')

def makeFileGlob(inDir, prefix, suffix, scanRange):
    ignore = [10086, 13638, 17751, 27083, 28089, 4564, 7165, 7167]
    filter=prefix+'*.'+suffix
    sdirs = sorted(glob.glob(os.path.join(inDir,filter)))
    dsc = [int(os.path.split(sdir)[1].split('_')[2].split('.')[0]) for sdir in sdirs]
    dfiles = []
    for i,ds in enumerate(dsc):
        if (ds >= scanRange[0]) & (ds <= scanRange[1]) & (ds not in ignore):
            dfiles.append(sdirs[i])            
    return dfiles


def L05_Pipeline(args, scanRange):
    L05_basedir = str(files('GUSTO_Pipeline') / 'level0.5')
    prefix = ['ACS5_', 'ACS3_']
    dirDataOut = args.path+'level0.5'
    os.makedirs(dirDataOut, exist_ok=True)
    inDir = args.path+'lags/'
    sum_files = 0
    if args.erase:
        eraseStr = '-e'
    else:
        eraseStr = ''
    if args.cpus:
        cpuStr = '-j '+args.cpus
        print('Number of cores used for processing: %i\n'%int(args.cpus))
    else:
        cpuStr = ''
    outPath = args.path+'level0.5'
    
    for band in args.band:
        print("Processing Band ", band)
        fileList = makeFileGlob(inDir, prefix[int(band)-1], 'dat', scanRange)
        with open(L05_basedir+"/filelist.txt", "w") as file:
            file.write("\n".join(fileList))
        print("Processing ", len(fileList), " files, please wait...")
        try:
            result=subprocess.run(['./run_level05.sh', eraseStr, cpuStr, '-f '+L05_basedir+'/filelist.txt', '-o '+outPath], cwd=L05_basedir, capture_output=True, text=True, check=True)
            logger.info(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error: {e}")
        sum_files += len(fileList)

    return sum_files
