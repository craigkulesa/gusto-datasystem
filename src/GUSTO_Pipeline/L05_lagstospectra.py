import os
import subprocess
from importlib.resources import files
from .Logger import *
from .DataIO import *

logger = logging.getLogger('pipelineLogger')

def L05_Pipeline(args, scanRange, verbose=False):
    global logger
    L05_basedir = str(files('GUSTO_Pipeline') / 'level0.5')
    prefix = ['ACS5_', 'ACS3_']
    dirDataOut = args.path+'level0.5'
    os.makedirs(dirDataOut, exist_ok=True)
    inDir = args.path+'lags/'
    sum_files = 0
    if args.erase:
        clear_folder(dirDataOut)
    if args.cpus:
        cpuStr = '-j '+args.cpus
        logger.info('Number of cores used for processing: %i'%int(args.cpus))
    else:
        cpuStr = ''
    outPath = args.path+'level0.5'
    
    for band in args.band:
        logger.info(f"Processing Band {band}")
        fileList = makeFileGlob(inDir, prefix[int(band)-1], 'dat', scanRange)
        with open(L05_basedir+"/filelist.txt", "w") as file:
            file.write("\n".join(fileList))
        logger.info(f"Processing {len(fileList)} files, please wait...")
        try:
            result=subprocess.run(['./run_level05.sh', cpuStr, '-f '+L05_basedir+'/filelist.txt', '-o '+outPath], cwd=L05_basedir, capture_output=True, text=True, check=True)
            logger.debug(result.stdout)
        except subprocess.CalledProcessError as e:
            logger.info(f"Error: {e}")
        sum_files += len(fileList)

    return sum_files
