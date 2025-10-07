"""
GUSTO Pipeline data file handling class.
This class handles the data file at the output of the Level 1 
pipeline and the Level 0.7 data files.

V. Tolls, SAO

created: 9/19/2024
"""

import datetime
import inspect
import os
import sys
import time
import warnings

from astropy.coordinates import SkyCoord
from astropy.io import ascii, fits
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning

import numpy as np
import numpy.ma as ma

import logging
log = logging.getLogger(__name__)

__updated___ = '20240919'
__version__ = 0.1
__date__ = '20240919'
fheader, tail = os.path.split(inspect.stack()[0][1])
__pyfile__ = tail
__path__ = fheader

warnings.filterwarnings('ignore', category=Warning,
                        message=' FITSFixedWarning: ', append=True)


def loadL08Data(ifile, verbose=False, usemask=False):
    """Function loading Level 0.8 data.

    Parameters
    ----------
    ifile : string
        full path to data L0.8 file
    verbose : boolean
        if TRUE, print info to STDOUT

    Returns
    -------
    Returns spectrum, data, and header arrays. Spectrum is the data['spec'] array 
    converted into a numpy masked array.
    """
    if verbose:
        print('Loading L0.8 data ...')
    
    with fits.open(ifile) as hdu:
        data   = hdu[1].data
        hdr    = hdu[0].header
        hdr1    = hdu[1].header

    keys = data.dtype.names
    if 'spec' in keys:
        dkey = 'spec'
    elif 'DATA' in keys:
        dkey = 'DATA'
        
    ss = np.argsort(data['UNIXTIME'])
    data = data[ss]
    
    if 'MIXER' not in keys:
        data['MIXER'] =  hdr1['MIXER']
    
    # combine spectrum and channel_flag mask
    if usemask == True:
        spec = ma.MaskedArray(data[dkey], mask=data['CHANNEL_FLAG'])
    else:
        spec = ma.MaskedArray(data[dkey], mask=np.zeros(data[dkey].size))
    return spec, data, hdr, hdr1


if __name__ == "__main__":

    # test summary:
    # test = 1:   reading Level 0.8 data files
    #

    test = 1

    if test == 1:
        print('No test yet')
        pass

