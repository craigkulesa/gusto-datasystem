import os
import glob
import subprocess
from importlib.resources import files
from .Logger import *

logger = logging.getLogger('pipelineLogger')

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
        clear_folder(dirDataOut)
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
            result=subprocess.run(['./run_level05.sh', cpuStr, '-f '+L05_basedir+'/filelist.txt', '-o '+outPath], cwd=L05_basedir, capture_output=True, text=True, check=True)
            logger.info(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error: {e}")
        sum_files += len(fileList)

    return sum_files
