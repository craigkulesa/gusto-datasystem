"""
GUSTO Pipeline data file handling class for GUSTO SDFITS
"""
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from importlib.resources import files
from .Logger import *
import os
import shutil
import glob
import numpy as np
import numpy.ma as ma
import subprocess

warnings.filterwarnings('ignore', category=Warning,
                        message=' FITSFixedWarning: ', append=True)

logger = logging.getLogger('pipelineLogger')

def clear_folder(folder_path):
    logger.warning('Erasing '+folder_path)
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            logger.warning(f"Failed to delete {file_path}. Reason: {e}")

            
def makeFileGlob(inDir, prefix, suffix, scanRange):
    #ignore = [10086, 13638, 17751, 27083, 28089, 4564, 7165, 7167]
    ignore = [0]
    filter=prefix+'*.'+suffix
    sdirs = sorted(glob.glob(os.path.join(inDir,filter)))
    dsc = [int(os.path.split(sdir)[1].split('_')[2].split('.')[0]) for sdir in sdirs]
    dfiles = []
    for i,ds in enumerate(dsc):
        if (ds >= scanRange[0]) & (ds <= scanRange[1]) & (ds not in ignore):
            dfiles.append(sdirs[i])            
    return dfiles



def runGitLog(level, file):
    gitdir = files('GUSTO_Pipeline')
    try:
        result = subprocess.run(['git', 'log', '-1', '--format=%cd', '--date=format-local:%Y-%m-%d %H:%M:%S %Z', '--pretty=format:Level '+level+' commit %h by %an %ad', '--', file], cwd=gitdir, capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        return f"Error: {e}"



def loadSDFITS(ifile, verbose=False, usemask=False):
    """Function loading GUSTO SDFITS files

    Parameters
    ----------
    ifile : string
        full path to data file
    verbose : boolean
        if TRUE, print info to STDOUT

    Returns
    -------
    Returns spectrum, data, and header arrays. Spectrum is the data['spec'] array 
    converted into a numpy masked array.
    """
    if verbose:
        logging.info('Loading SDFITS data ...')
    
    with fits.open(ifile) as hdu:
        data   = hdu[1].data
        hdr    = hdu[0].header
        hdr1    = hdu[1].header

    keys = data.dtype.names
        
    ss = np.argsort(data['UNIXTIME'])
    data = data[ss]
    
    if 'MIXER' not in keys:
        data['MIXER'] =  hdr1['MIXER']
    
    spec = ma.MaskedArray(data['DATA'], mask=np.zeros(data['DATA'].size))
    return spec, data, hdr, hdr1

