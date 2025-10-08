"""
GUSTO Pipeline data file handling class for GUSTO L1 SDFITS
"""

import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
import numpy as np
import numpy.ma as ma

warnings.filterwarnings('ignore', category=Warning,
                        message=' FITSFixedWarning: ', append=True)


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
        print('Loading SDFITS data ...')
    
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

