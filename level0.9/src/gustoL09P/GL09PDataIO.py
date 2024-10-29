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
from pprint import pprint
import re
import sys
import time
import warnings

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii, fits
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS, WCSSUB_LONGITUDE, WCSSUB_LATITUDE, WCSSUB_SPECTRAL
from astropy.wcs import validate as WCS_validate
import astropy.wcs
import pkg_resources
from tqdm import tqdm

import astropy.constants as c
import astropy.constants as const
import astropy.units as u
from gustoL09P.GL09PParameters import GL09PParameters
import numpy as np
import numpy.ma as ma
import parsl
from parsl.app.app import python_app, bash_app
from parsl.channels import LocalChannel
from parsl.config import Config
from parsl.configs.local_threads import config
from parsl.executors import HighThroughputExecutor
from parsl.providers import LocalProvider

__updated___ = '20240919'
__version__ = 0.1
__date__ = '20240919'
fheader, tail = os.path.split(inspect.stack()[0][1])
__pyfile__ = tail
__path__ = fheader

# warnings.filterwarnings('ignore', category=Warning, message=' FITSFixedWarning: The WCS transformation has more axes (3) than the image it is associated with (2) [astropy.wcs.wcs]', append=True)
warnings.filterwarnings('ignore', category=Warning,
                        message=' FITSFixedWarning: ', append=True)


def loadL08Data(ifile, verbose=False):
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

    # header keys
    # ['SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 'COMMENT', 'COMMENT', 'CALID', 'CUNIT1', 'CRPIX1', 'CRVAL1', 'CDELT1', 'HKSCANID', 'TELESCOP', 
    #  'LINE', 'LINEFREQ', 'BAND', 'NPIX', 'DLEVEL', 'PROCTIME', 'SER_FLAG', 'COMMENT', 'CRADLE02', 'CRYCSEBK', 'CRYOPORT', 'CALMOTOR', 'CRADLE03', 
    #  'QAVCCTRL', 'COOLRTRN', 'FERADIAT', 'CRYCSEFT', 'CRADLE04', 'OAVCCTRL', 'COOLSUPL', 'CRADLE01', 'EQUILREF', 'SECONDRY', 'COMMENT', 'B1_SYNTH', 
    #  'B1_PWR_3', 'B1_PWR_4', 'B1M5_AMP', 'B1_PWR_1', 'B1_PWR_2', 'B2_UCTRL', 'B2MLTDRV', 'B2_PWR_3', 'B2_PWR_4', 'B2AVA183', 'B1M5MULT', 'B2M5_AMP', 
    #  'B2_PWR_1', 'B2_PWR_2', 'COMMENT', 'T_IS', 'T_IVCS', 'T_LNA', 'T_MIXER', 'T_OS', 'T_OVCS', 'T_QCL', 'T_TANK', 'COMMENT', 'GOND_ALT', 'GOND_LAT', 
    #  'GOND_LON', 'ELEVATON', 'OBJECT', 'IF0', 'SYNTFREQ', 'SYNTMULT', 'VLSR', 'COMMENT']

    # data table
    # (numpy.record, [('MIXER', '>i4'), ('NINT', '>i4'), ('UNIXTIME', '>f8'), ('NBYTES', '>i4'), ('CORRTIME', '>i4'), ('INTTIME', '>f4'), 
    # ('ROW_FLAG', '>i4'), ('Ihigh', '>i4'), ('Qhigh', '>i4'), ('Ilow', '>i4'), ('Qlow', '>i4'), ('Ierr', '>i4'), ('Qerr', '>i4'), 
    # ('VIhi', '>f4'), ('VQhi', '>f4'), ('VIlo', '>f4'), ('VQlo', '>f4'), ('scanID', '>i4'), ('subScan', '>i4'), ('scan_type', 'S6'), 
    # ('THOT', '>f4'), ('RA', '>f4'), ('DEC', '>f4'), ('filename', 'S48'), ('PSat', '>f4'), ('Imon', '>f4'), ('Gmon', '>f4'), 
    # ('spec', '>f4', (1024,)), ('CHANNEL_FLAG', '>i2', (1024,))])
    
    if verbose:
        print(list(hdr.keys()))
        print(type(data))
        print(data.dtype)
        print('spectra: ', data['spec'].shape)
        print('row flag: ', data['CHANNEL_FLAG'].shape)
        print('channel mask: ', data['CHANNEL_FLAG'].shape)

    # combine spectrum and channel_flag mask
    spec = ma.MaskedArray(data['spec'], mask=data['CHANNEL_FLAG'])
    spec = ma.masked_invalid(spec)
    spec.mask[:,612:1023] = 1

    return spec, data, hdr




if __name__ == "__main__":

    # test summary:
    # test = 1:   reading Level 0.8 data files
    #

    test = 1

    if test == 1:
        pass

