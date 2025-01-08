#!/usr/bin/env python
"""
This is the executable for the GUSTO L2 Pipeline.
"""

__date__ = '12/12/2024'
__updated__ = '20241212'
__version__ = '0.1'
__author__ = 'V. Tolls, CfA | Harvard & Smithsonian'

from joblib import Parallel, delayed
from joblib import Memory
import glob
import numpy as np
import time
import pkg_resources
import parsl
import sys
import os
import logging
import matplotlib.pyplot as plt
import multiprocessing
from matplotlib.patches import Rectangle
from pathlib import Path
from pprint import pprint
from multiprocessing.pool import Pool
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.linalg import cholesky
from scipy.stats import norm

from astropy import constants as const
from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time

from .GL09PDataIO import loadL08Data
from .GL09PUtils import *
from .GL09PLogger import *

import logging
log = logging.getLogger(__name__)

offsetsfile0 = pkg_resources.resource_filename('gustoL09P', 'Data/offsets.txt')


def GL10Pipeline(cfi, scanRange, verbose=False):
    """Function processing the Level 0.95 data injecting 
    coordinate corrections. 

    1: read in the data from file
    2: apply coordinate offset correction to each spectrum
    3: save back to files with only setting a flag about the baseline correction,
    but all other parameters are passed through


    Parameters
    ----------
    scanRange : int array

    Returns
    -------

    """
    
        #logger = logging.getLogger('GL09PLogger')
    
    if cfi['gprocs']['debug']==True:
        print('\nExecuting debug mode.\n')
        logger = multiprocessing.log_to_stderr()
        logger.setLevel(multiprocessing.SUBDEBUG)
        n_procs = 1
    else:
        n_procs = multiprocessing.cpu_count() - 2
        
    print('Number of cores used for processing: %i\n'%(n_procs))
    
    # get lines for processing
    lines = cfi['gprocs']['lines'].replace('[','').replace(']','').replace(' ','').split(',')
    print('Lines: ', lines[0])
    lines= ['CII', 'NII']
    
    ignore = [0]
    
    # read offsets file
    offsetsfile = cfi['gprocs']['offsetsfile']
    if offsetsfile == '':
        offsetsfile = offsetsfile0
    if not os.path.exists(offsetsfile):
        print('offsetsfile: ', offsetsfile)
        print('No file with coordinate offsets available. Terminating pipeline run!')
        sys.exit(1)
        
    offs = np.genfromtxt(offsetsfile, delimiter='\t', skip_header=2, 
                     dtype=[('mxpix', 'U4'), ('az', '<f8'), ('el', '<f8'), ('comment', 'U16')])

    
    for line in lines:
        if verbose:
            print(line)
        # identify the files for processing
        inDir = cfi['gdirs']['L095DataDir']
        outDir = cfi['gdirs']['L10DataDir']
        if line=='CII':
            filter = 'ACS5*.fits'
        else:
            filter = 'ACS3*.fits'
        print('outDir: ', outDir)
        print('filter: ', os.path.join(inDir,filter))
        
        # sdirs = sorted(glob.glob(os.path.join(inDir,filter), root_dir=inDir))
        #print(glob.glob(os.path.join(inDir,filter)))
        sdirs = sorted(glob.glob(os.path.join(inDir,filter)))
        #print('single result: ', sdirs[0], os.path.split(sdirs[0]))
        dsc = [int(os.path.split(sdir)[1].split('_')[1].split('.')[0]) for sdir in sdirs]
        
        sdirs.sort(key=lambda sdirs: dsc)
        
        dfiles = []
        for i,ds in enumerate(dsc):
            if (ds >= scanRange[0]) & (ds <= scanRange[1]) & (ds not in ignore):
                dfiles.append(sdirs[i])
                        
        n_ds = len(dfiles)
        if int(cfi['gprocs']['max_files']) > 0:
            n_ds = int(cfi['gprocs']['max_files'])
            dfiles = dfiles[:n_ds]
                    
        paramlist = [[a, b, c, d, e, f] for a in [line] for b in [inDir] for c in [outDir] for d in dfiles for e in [offs] for f in [bool(cfi['gprocs']['debug'])]]
        # paramlist = [[a, b, c, d, e] for a in [line] for b in [inDir] for c in [outDir] for d in dfiles for e in worker_configurer]
        #print(paramlist)
        if verbose:
            print('Number of data files: ', n_ds, len(sdirs))
            #print('Selected data files: ', dfiles)
        
        
        # setup multiprocessing loop here to process each file in list
        with Pool(n_procs) as pool:
            # execute tasks in order
            for result in pool.imap(processL10, paramlist):
                print(f'Processed: {result}', flush=True)
        
    return n_ds



def processL10(params, verbose=True):
    """Function applying the actual baseline fit


    Parameters
    ----------
    params : list
            list of parameters used in process

    Returns
    -------

    """
    
    #loadL08Data(dfile, verbose=True)
    line, inDir, outDir, dfile, offsets, debug = params[0], params[1], params[2], params[3], params[4], params[5]
    
    # define some processing data first (maybe relocat to function later?)

    if 'ACS3' in dfile:
        pixel_cut = 600
        band = 2
        add = 'B2'
    else:
        pixel_cut = 300
        band = 1
        add = 'B1'

    datavalid = True

    #logger.info('loading: ', os.path.join(inDir,dfile), ' for line: ', line)
    spec, data, hdr, hdr1 = loadL08Data(os.path.join(inDir,dfile), verbose=False)
    print(os.path.join(inDir,dfile))
    rowFlag = data['ROW_FLAG']
    
    n_spec, n_pix = spec.shape
    
    prange = [int(hdr['pgpixst']), int(hdr['pgpixen'])]
    
    # osel = np.argwhere((otfID == data['scanID']) & (otfID.size>=1) & (rfsID.size>2) & (rhsID.size>2) & (hotID.size>2) & (mix == data['MIXER']) & (data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0))
    # osel = np.argwhere((data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0)).flatten()
    
    # if len(osel) <= 0:
    #     print('WARNING: No OTF spectra available.')
    #     # logger.warning('No OTF spectra available.')
    #     datavalid = False
    #     return 0

    osel = np.argwhere((data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0)).flatten()
    
    umixers = np.unique(data['MIXER'])
    
    
    print('processing data of %s...'%(dfile))
    # insert the coordinate corrections
    # Note: the coordinate correction is not yet final
    # and will be (iteratively) improved
    
    mxoffs = getMixerOffsets(band, umixers, verbose=verbose)
    
    for i, mix in enumerate(umixers):
        azoff = mxoffs['az'][i]
        eloff = mxoffs['el'][i]
        
        msel = np.argwhere((data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0) & (data['MIXER']==mix)).flatten()
        if msel.size > 0:
            ras = data['RA'][msel]
            decs = data['DEC'][msel]
            utime = data['UNIXTIME'][msel]
            
            glat = hdr['GOND_LAT']
            glon = hdr['GOND_LON']
            galt = hdr['ELEVATON']
            
            cc = SkyCoord(ra=ras*u.deg, dec=decs*u.deg, frame='icrs')
            
            balloon = EarthLocation(lat=glat, lon=glon, height=galt*u.m)
            
            otime = Time(utime, format='unix')
            
            aa = AltAz(location=EarthLocation(lat=glat, lon=glon, height=galt*u.m), obstime=otime)
            altaz = cc.transform_to(aa)
            
            # apply beam offsets
            naz = altaz.az.deg*u.deg + azoff*u.deg
            nalt = altaz.alt.deg*u.deg + eloff*u.deg
            
            ncc = SkyCoord(AltAz(az=naz, alt=nalt, obstime=otime, location=balloon))
            nradec = ncc.transform_to('icrs')
            
            data['RA'][msel] = nradec.ra.deg
            data['DEC'][msel] = nradec.dec.deg
        

    # now we have to save the data in a FITS file
    
    # insert delta_ra and delta_dec in table for compatibility with WCS
    # insert WCS information in header
    # all will be relative to the first ra/dec pair
    cc0 = SkyCoord(data['ra'][osel[0]]*u.deg, data['dec'][osel[0]]*u.deg, frame='icrs')
    cc = SkyCoord(data['ra']*u.deg, data['dec']*u.deg, frame='icrs')    
    dra, ddec = cc0.spherical_offsets_to(cc)
    # add the columns to the data
    cols = data.columns
    new_cols = fits.ColDefs([
        fits.Column(name='dra', format='D', array=dra.deg),
        fits.Column(name='ddec', format='D', array=ddec.deg)])
    nhdu = fits.BinTableHDU.from_columns(cols + new_cols)
    data = nhdu.data

    # select only the processed good data
    odata = data[osel]
    
    # add 
    hdr.set('CDELT2', value=0.000000001, comment=(''), after='CDELT1')
    hdr.set('CRVAL2', value=(data['ra'][osel[0]]*u.deg).value, comment=(''), after='CDELT1')
    hdr.set('CRPIX2', value=0, comment=(''), after='CDELT1')
    hdr.set('CUNIT2', value='deg', comment=(''), after='CDELT1')
    hdr.set('CTYPE2', value='RA---GLS', comment=(''), after='CDELT1')
    hdr.set('CDELT3', value=0.000000001, comment=(''), after='CDELT2')
    hdr.set('CRVAL3', value=(data['dec'][osel[0]]*u.deg).value, comment=('CDELT2'))
    hdr.set('CRPIX3', value=0, comment=(''), after='CDELT2')
    hdr.set('CUNIT3', value='deg', comment=(''), after='CDELT2')
    hdr.set('CTYPE3', value='DEC--GLS', comment=(''), after='CDELT2')
    
    tred = Time(datetime.datetime.now()).fits
    
    # hdr.insert('VLSR', ('PROC_LEV', 0.95, 'pipeline processing level'), after=True)
    # hdr.add_comment('Pipeline Processing', before='PROC_LEV')
    # hdr.insert('PROC_LEV', ('PROCDATE', tred.split('T')[0], 'Date of processing'))
    # hdr.insert('PROCDATE', ('PROCTIME', tred.split('T')[1], 'Time of processing'))
    hdr['DLEVEL'] = 1.0
    hdr['PROCTIME'] = tred
    
    hdr.set('', value='', after='BS_ITERM')
    hdr.set('', value='          Level 0.95 Pipeline Processing', after='BS_ITERM')
    hdr.set('', value='', after='BS_ITERM')
    hdr.set('L10PTIME', value=tred, comment=('L1.0 pipeline processing time'))
    
    os.makedirs(outDir, exist_ok=True)
    ofile = os.path.join(outDir, os.path.split(dfile)[1].replace('_L095.fits','_L10.fits'))
    fits.writeto(ofile, data=None, header=hdr, overwrite=True)
    fits.append(ofile, data=odata, header=hdr1)

    print('saved file: ', ofile)
    # logger.info('saved file: ', ofile)
    
    return dfile

   

def getMixerOffsets(band, mixers, type='THEORY', offsetfile=None, verbose=False):
    """Function retrieving the GUSTO on-sky mixer offsets from file.

    usage:
    ------
    
    aa = getMixerOffsets(1, [3, 5, 8], verbose=True)
    print(aa['az'])     # prints: [0.06079  0.062584 0.093356]
    print(aa[0]['az'])  # prints: 0.06079
    print(aa['az'][0])  # prints: 0.06079
    
    
    Parameters
    ----------
    band : str
        GUSTO band: 'B1' or 'B2'
    mixers : int array
        array of mixer indices
    type : str
        type of offset determination: 'Theory' or 'AS_MEASURED'
    offsetfile : str
        path to file with mixer offsets
    verbose : bool
        flag for output to STDOUT
        

    Returns
    -------
    function returns rec array with mixer offsets

    """

    if band==1:
        # [NII]
        cband = 'B1'
    elif band==2:
        # [CII]
        cband = 'B2'

    if offsetfile is None:
        offsetfile = offsetsfile0
    # if verbose:
    #     print('\nCoordinate offset file: ', offsetfile)

    data = np.genfromtxt(offsetfile, delimiter='\t', skip_header=2, 
                         dtype=[('mxpix', 'U4'), ('az', '<f8'), ('el', '<f8'), ('type', 'U16')])
    
    cmixers = ['B%iM%i'%(band, i) for i in mixers]
    
    ar = [np.argwhere((cmixer == data['mxpix'])&((type==data['type'])|(data['type']=='FIDUCIAL'))).flatten() for cmixer in cmixers]

    # if verbose:
    #     print('cmixers: ', cmixers)
    #     print('args: ', ar)
    #     print('data: \n', data[ar].flatten())
        
    return data[ar].flatten()



