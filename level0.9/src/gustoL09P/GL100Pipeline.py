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
from .GL09PDataIO import loadL08Data
from .GL09PUtils import *
from .GL09PLogger import *
import logging
log = logging.getLogger(__name__)

offsetsfile0 = pkg_resources.resource_filename('gustoL09P', 'Data/offsets.txt')


def GL100Pipeline(cfi, scanRange, verbose=False):
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
    offsetsfile = cfi['gdirs']['L09DataDir']
    if offsetsfile is None:
        offsetsfile = offsetsfile0
    if not os.path.exists(offsetsfile):
        print('No file with coordinate offsets available. Terminating pipeline run!')
        os.exit(1)
        
    offs = np.genfromtxt(offsetsfile, delimiter='\t', skip_header=2, 
                     dtype=[('mxpix', 'U4'), ('az', '<f8'), ('el', '<f8'), ('comment', 'U16')])

    
    for line in lines:
        if verbose:
            print(line)
        # identify the files for processing
        inDir = cfi['gdirs']['L095DataDir']
        outDir = cfi['gdirs']['L1DataDir']
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
                    
        paramlist = [[a, b, c, d, e, f] for a in [line] for b in [inDir] for c in [outDir] for d in dfiles for e in [offsets] for f in [bool(cfi['gprocs']['debug'])]]
        # paramlist = [[a, b, c, d, e] for a in [line] for b in [inDir] for c in [outDir] for d in dfiles for e in worker_configurer]
        #print(paramlist)
        if verbose:
            print('Number of data files: ', n_ds, len(sdirs))
            #print('Selected data files: ', dfiles)
        
        
        # setup multiprocessing loop here to process each file in list
        with Pool(n_procs) as pool:
            # execute tasks in order
            for result in pool.imap(processL100, paramlist):
                print(f'Processed: {result}', flush=True)
        
    return n_ds



def processL100(params, verbose=True):
    """Function applying the actual baseline fit


    Parameters
    ----------
    params : list

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
    rowFlag = data['ROW_FLAG']
    
    n_spec, n_pix = spec.shape
    
    prange = [int(hdr['pgpixst']), int(hdr['pgpixen'])]
    
    # osel = np.argwhere((otfID == data['scanID']) & (otfID.size>=1) & (rfsID.size>2) & (rhsID.size>2) & (hotID.size>2) & (mix == data['MIXER']) & (data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0))
    # osel = np.argwhere((data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0)).flatten()
    
    if len(osel) <= 0:
        print('WARNING: No OTF spectra available.')
        # logger.warning('No OTF spectra available.')
        datavalid = False
        return 0

    spec_OTF = np.squeeze(spec[osel,:])
    n_OTF, n_otfpix = spec_OTF.shape 
    
    umixers = np.unique(data['MIXER'])
    
    
    print('processing data ...')
    # create the calibrated spectra
    for i0 in range(n_OTF):
        # insert the coordinate corrections
        # Note: the coordinate correction is not yet final
        # and will be (iteratively) improved
        
        mxoffs = getMixerOffsets(band, umixers, verbose=verbose)
        
        for i, mix in enumerate(mixers):
            azoff = mxoffs['az'][i]
            eloff = mxoffs['el'][i]
            
            osel = np.argwhere((data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0) & (data['MIXER']==mix)).flatten()
            ras = data['RA'][osel]
            decs = data['DEC'][osel]
            utime = data1['UNIXTIME']
            
            glat = hdr['GOND_LAT']
            glon = hdr['GOND_LON']
            galt = hdr['ELEVATON']
            
            cc = SkyCoord(ra=ras*u.deg, dec=decs*u.deg, frame='icrs')
            
            balloon = EarthLocation(lat=glat, lon=glon, height=galt*u.m)
            
            otime = Time(utime, format='unix')
            
            aa = AltAz(location=EarthLocation(lat=glat, lon=glon, height=galt*u.m), obstime=otime)
            altaz = cc.transform_to(aa)
               
            # get beam offsets
            az_bm = 0.030104  #azoffset[band, mix]
            alt_bm = 0.000023  #eloffset[band, mix]
            
            # apply beam offsets
            naz = altaz.az.deg*u.deg + azoff*u.deg
            nalt = altaz.alt.deg*u.deg + altoff*u.deg
            
            ncc = SkyCoord(AltAz(az=naz, alt=nalt, obstime=otime, location=balloon))
            nradec = ncc.transform_to('icrs')
            
            data['RA'][osel] = nradec.ra.deg
            data['DEC'][osel] = nradec.dec.deg
        
        pass

    # now we have to save the data in a FITS file

    data['CHANNEL_FLAG'] [osel,:] = spec_OTF.mask
    
    osel = np.argwhere((data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0)).flatten()
    odata = data[osel]
    
        
    tred = Time(datetime.datetime.now()).fits
    
    # hdr.insert('VLSR', ('PROC_LEV', 0.95, 'pipeline processing level'), after=True)
    # hdr.add_comment('Pipeline Processing', before='PROC_LEV')
    # hdr.insert('PROC_LEV', ('PROCDATE', tred.split('T')[0], 'Date of processing'))
    # hdr.insert('PROCDATE', ('PROCTIME', tred.split('T')[1], 'Time of processing'))
    hdr['PROC_LEV'] = 1.0
    hdr['PROCDATE'] = tred.split('T')[0]
    hdr['PROCTIME'] = tred.split('T')[1]
    
    os.makedirs(outDir, exist_ok=True)
    ofile = os.path.join(outDir, os.path.split(dfile)[1].replace('.fits','_L095.fits'))
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

    data = np.genfromtxt(offsetfile, delimiter='\t', skip_header=2, 
                         dtype=[('mxpix', 'U4'), ('az', '<f8'), ('el', '<f8'), ('type', 'U16')])
    
    cmixers = ['B%iM%i'%(band, i) for i in mixers]
    
    ar = [np.argwhere((cmixer == data['mxpix'])&(type==data['type'])).flatten() for cmixer in cmixers]

    if verbose:
        print('cmixers: ', cmixers)
        print('args: ', ar)
        print('data: \n', data[ar].flatten())
        
    return data[ar].flatten()



