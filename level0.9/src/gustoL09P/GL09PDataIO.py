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
# import parsl
# from parsl.app.app import python_app, bash_app
# from parsl.channels import LocalChannel
# from parsl.config import Config
# from parsl.configs.local_threads import config
# from parsl.executors import HighThroughputExecutor
# from parsl.providers import LocalProvider

import logging
log = logging.getLogger(__name__)

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
        hdr1    = hdu[1].header

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
    
    keys = data.dtype.names
    if 'spec' in keys:
        dkey = 'spec'
    elif 'DATA' in keys:
        dkey = 'DATA'
        
    ss = np.argsort(data['UNIXTIME'])
    data = data[ss]
    
    if 'MIXER' not in keys:
        data['MIXER'] =  hdr1['MIXER']
    
    # if verbose:
    #     print(list(hdr.keys()))
    #     print(type(data))
    #     print(data.dtype)
    #     print('spectra: ', data[dkey].shape)
    #     print('row flag: ', data['CHANNEL_FLAG'].shape)
    #     print('channel mask: ', data['CHANNEL_FLAG'].shape)

    # combine spectrum and channel_flag mask
    spec = ma.MaskedArray(data[dkey], mask=data['CHANNEL_FLAG'])
    #spec = ma.masked_invalid(spec)
    #spec.mask[:,612:1023] = 1

    return spec, data, hdr, hdr1


def saveProcData():
    """This is not working code yet. Old code is just being parked until 
    made usable.
    """

    # Table for data extension:
    # ra = data['RA'][osel]
    # dec = data['DEC'][osel]
    # ra = data['RA']
    # dec = data['DEC']
    # rflag = data['ROW_FLAG'][osel] 
    #
    # c0 = SkyCoord(ra[0], dec[0], unit=(u.deg, u.deg), frame='icrs')
    # c1 = SkyCoord(ra[osel]*u.degree, dec[osel]*u.degree, frame='icrs')
    # rc1 = c1.transform_to(SkyOffsetFrame(origin=c0)) 
    # raoff = rc1.lon.deg
    # decoff = rc1.lat.deg 
    
    tsyseff_avg = np.nanmean(tsyseff[:,200:400], axis=1)
    # timeobs = Time(data['UNIXTIME'][osel], format='unix').fits
    # dobs = [tt[0].split('T')[0] for tt in timeobs]
    # tobs = [tt[0].split('T')[1] for tt in timeobs]
    tred = Time(datetime.datetime.now()).fits
    
    # aa = np.zeros(n_OTF)
    # ao = np.ones(n_OTF)
    
    
    # changing the output data format to
    # primarily passing through the header and the data table
    # and only few variables and the header keys are updated or added
    
    # this is the offset from the longitude in CRVAL2 -> RA
    col1 = Column(raoff, name='CDELT2', description='longitude')
    # this is the offset from the longitude in CRVAL3 -> Dec
    col2 = Column(decoff, name='CDELT3', description='latitude')
    col3 = Column(ta, name='DATA', description='spectrum')
    col4 = Column(tsyseff_avg, name='TSYS', description='system temperature (K)')
    col5 = Column(aa, name='IMAGFREQ', description='image band frequency (Hz)')
    col6 = Column(dobs, name='DATE-OBS', description='date of obs')
    col7 = Column(tobs, name='UT', description='time obs started')
    col8 = Column(ao*float(hdr['ELEVATON']), name='ELEVATIO', description='')
    col9 = Column(ao*(Thot[0]+Thot[1])/2., name='THOT', description='hot load temperature (K)')
    col10 = Column(ao*Tsky, name='TCOLD', description='Sky temperature (K)')
    col11 = Column(ta.mask, name='CHANNEL_FLAG', description='pixel mask')
    col12 = Column(rflag, name='Row_flag', description='flagging compromised spectra')
    fT = Table([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12])
        
    # now, create the header
    

    # NAXIS1: number of channels
    # NAXIS2: pos coord1
    # NAXIS3: pos coord2
    # NAXIS4: Stokes
    wcs_dict = {
        "CTYPE1": "Velocity    ",
        "CUNIT1": "m/s",
        "CDELT1": 0.2,
        "CRPIX1": 0.0,
        "CRVAL1": 0.0,
        "NAXIS1": n_pix,
        "CTYPE2": "RA---GLS",
        "CUNIT2": "deg",
        "CDELT2": 0.1,
        "CRPIX2": 0.0,
        "CRVAL2": ra[0],
        "NAXIS2": 1,
        "CTYPE3": "DEC--GLS",
        "CUNIT3": "deg",
        "CDELT3": 0.1,
        "CRPIX3": 0.0,
        "CRVAL3": dec[0],
        "NAXIS3": 1,
        # "CTYPE4": "STOKES  ",
        # "CUNIT4": "",
        # "CDELT4": 0.0,
        # "CRPIX4": 0.0,
        # "CRVAL4": 1,
        # "NAXIS4": 1,
    }
    wcs = astropy.wcs.WCS(wcs_dict)
    header = wcs.to_header()
    
    # additional header keywords
    header['TELESCOP'] = hdr['TELESCOP']
    header['BAND'] = hdr['BAND']
    header['LINE'] = line
    header['OBJECT'] = 'GUSTO_GPS'
    header['FOFFSET'] = 0.0
    header['RESTFREQ'] = rfreq
    header['LINEFREQ'] = hdr['LINEFREQ']
    header['IMAGFREQ'] = rfreq
    header['IF0'] = hdr['IF0']
    header['SYNTFREQ'] = hdr['SYNTFREQ']
    header['VELO-LSR'] = hdr['VLSR']    # m/s
    header['VELDEF'] = 'RADI-LSR'
    header['DELTAV'] = 0.2
    header['GAINIMAG'] = 0.99999E-2
    header['BEAMEFF'] = 0.95
    header['FORWEFF'] = 0.95
    header['EPOCH'] = 2000
    header['DATE-OBS'] = tobs[0][0]
    header['DATE-RED'] = tred
    header['UT'] = 0.0
    header['OBSTIME'] = 0.1
    header['CREATOR'] = 'GUSTO Level 0.9 Pipeline'
    header['HISTORY'] = 'pipeline version: %s'%(__version__)
    header['HISTORY'] = 'last pipeline update: %s'%(__updated__)
    header['HISTORY'] = 'file created: %s'%(time.strftime("%c"))
    header['CALID'] = hdr['CALID']
    header['DLEVEL'] = '0.9'
    header['PROCTIME'] = time.strftime("%c")
    header['SER_FLAG'] = hdr['SER_FLAG']
    header['T_MIXER'] = hdr['T_MIXER']
    header['GOND_ALT'] = hdr['GOND_ALT']
    header['GOND_LAT'] = hdr['GOND_LAT']
    header['GOND_LON'] = hdr['GOND_LON']
    header['ELEVATON'] = hdr['ELEVATON']
    # header[''] = hdr['']
    # header[''] = hdr['']
    
    
    # now write the fits file
    hdu = fits.PrimaryHDU()
    hdu2 = fits.BinTableHDU(data=fT, header=header, name='MATRIX', ver=True)
    hdulist = fits.HDUList([hdu, hdu2])
    os.makedirs(outDir, exist_ok=True)
    ofile = os.path.join(outDir, os.path.split(dfile)[1].replace('.fits','_%s_L09.fits'%(mix)))
    hdulist.writeto(ofile, overwrite=True)
    
    #if verbose:
    print('saved file: ', ofile)
    # logger.info('saved file: ', ofile)


    # bunit = 'K (Ta)'
    # RADESYS = 'FK5'
    # EQUINOX = 2000
    # RA = 
    # DEC = 
    # LINE = 'CII'
    # RESTFREQ = 
    # IMAGFREQ = 
    # SPECSYS = 'LSRK'
    # OBJECT = 


if __name__ == "__main__":

    # test summary:
    # test = 1:   reading Level 0.8 data files
    #

    test = 1

    if test == 1:
        print('No test yet')
        pass

