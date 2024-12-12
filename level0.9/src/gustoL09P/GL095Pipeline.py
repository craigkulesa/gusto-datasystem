#!/usr/bin/env python
"""
This is the executable for the GUSTO L2 Pipeline.
"""

__date__ = '12/06/2024'
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
from pybaselines import Baseline, utils
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.linalg import cholesky
from scipy.stats import norm
from .GL09PDataIO import loadL08Data
from .GL09PUtils import *
from .GL09PLogger import *
import logging
log = logging.getLogger(__name__)


def GL095Pipeline(cfi, scanRange, verbose=False):
    """Function processing the Level 0.9 data and perform a
    baseline correction. 

    1: read in the data from file
    2: apply bseline correction to each spectrum
    3: save back to files with only setting a flag about the baseline correction,
    but all other parameters are passed through


    Parameters
    ----------
    scanRange : int

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
    
    
    for line in lines:
        if verbose:
            print(line)
        # identify the files for processing
        inDir = cfi['gdirs']['L09DataDir']
        outDir = cfi['gdirs']['L095DataDir']
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
            
        
        paramlist = [[a, b, c, d, e, f] for a in [line] for b in [inDir] for c in [outDir] for d in dfiles for e in [int(cfi['gprocs']['drmethod'])] for f in [bool(cfi['gprocs']['debug'])]]
        # paramlist = [[a, b, c, d, e] for a in [line] for b in [inDir] for c in [outDir] for d in dfiles for e in worker_configurer]
        #print(paramlist)
        if verbose:
            print('Number of data files: ', n_ds, len(sdirs))
            #print('Selected data files: ', dfiles)
        
        
        # setup multiprocessing loop here to process each file in list
        with Pool(n_procs) as pool:
            # execute tasks in order
            for result in pool.imap(processL09, paramlist):
                print(f'Processed: {result}', flush=True)
        
    return n_ds



def processL09(params, verbose=True):
    """Function applying the actual baseline fit


    Parameters
    ----------
    params : list

    Returns
    -------

    """
    
    #loadL08Data(dfile, verbose=True)
    line, inDir, outDir, dfile, drmethod, debug = params[0], params[1], params[2], params[3], params[4], params[5]
    
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
    osel = np.argwhere((data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0)).flatten()
    
    if len(osel) <= 0:
        print('WARNING: No OTF spectra available.')
        # logger.warning('No OTF spectra available.')
        datavalid = False
        return 0

    spec_OTF = np.squeeze(spec[osel,:])
    n_OTF, n_otfpix = spec_OTF.shape 
    
    basecorr = np.zeros(spec_OTF.shape)
    
    print('processing data ...')
    # create the calibrated spectra
    for i0 in range(n_OTF):
        # perform a baseline correction on all OTF spectra
        
        # this call might require to address the baseline correction method
        # possible methods are 1); polunomial and 2) 
        #basecorr[i0,:] = baseCorrection()
        # only correct the good part of the spectrum
        # try to mask bad spectral pixels using weights set to zero: TBI
        
        bs, ws = arplsw(spec_OTF[i0,prange[0]:prange[1]], lam=1e2, ratio=0.02, itermax=100)
        
        basecorr[i0,prange[0]:prange[1]] = bs
        
        spec_OTF[i0,prange[0]:prange[1]] -= basecorr[i0,prange[0]:prange[1]]

    # now we have to save the data in a FITS file

    keys = data.dtype.names
    if 'spec' in keys:
        dkey = 'spec'
    elif 'DATA' in keys:
        dkey = 'DATA'
        
    data[dkey][osel,:] = spec_OTF.data
    data['CHANNEL_FLAG'] [osel,:] = spec_OTF.mask
    
        
    tred = Time(datetime.datetime.now()).fits
    
    # hdr.insert('VLSR', ('PROC_LEV', 0.95, 'pipeline processing level'), after=True)
    # hdr.add_comment('Pipeline Processing', before='PROC_LEV')
    # hdr.insert('PROC_LEV', ('PROCDATE', tred.split('T')[0], 'Date of processing'))
    # hdr.insert('PROCDATE', ('PROCTIME', tred.split('T')[1], 'Time of processing'))
    hdr['PROC_LEV'] = 0.95
    hdr['PROCDATE'] = tred.split('T')[0]
    hdr['PROCTIME'] = tred.split('T')[1]
    
    os.makedirs(outDir, exist_ok=True)
    ofile = os.path.join(outDir, os.path.split(dfile)[1].replace('.fits','_L095.fits'))
    fits.writeto(ofile, data=None, header=hdr, overwrite=True)
    fits.append(ofile, data=data, header=hdr1)
    
    # if in debug mode, add more extensions
    # if debug & datavalid:
    #     # data needed
    #     # Tsys1, Tsys2, fraca, fracb, Tsyseff, hcorr, spref
    #     # we have spectra/data for each otf spectrum: Tsyseff, hcorr, spref, fraca, fracb
    #     # single spectra: Tsys1, Tsys2, 
    #     # small groupf of spectra: hots
    #     col1 = Column(aTsyseff, name='Tsyseff', description='effective Tsys per OTF spectrum')
    #     col2 = Column(ahcorr, name='hcorr', description='effective hot spectrum')
    #     col3 = Column(aspref, name='spref', description='effective hot for ref spectrum')
    #     col4 = Column(aspref2, name='spref2', description='ref spectrum')
    #     col5 = Column(afrac, name='frac', description='fraction of Tsys1/2')
    #     col6 = Column(spec, name='spec', description='original spectra')
    #     col7 = Column(aTa2, name='tant', description='spectrum without hot correction')
    #     fT = Table([col1, col2, col3, col4, col5, col6, col7])
    #     fits.append(ofile, data=fT.as_array())
    #
    #     col21 = Column(np.vstack(aghots), name='hots', description='averaged hot spectra')
    #     col22 = Column(np.hstack(aghtim), name='htime', description='utime of avg. hots')
    #     col23 = Column(np.hstack(aghmix), name='hmixer', description='mixer of avg. hots')
    #     col24 = Column(np.hstack(aghtint), name='htint', description='integration time of avg. hots')
    #     fgh = Table([col21, col22, col23, col24])
    #     fgh.write(ofile, append=True)
    #
    #     col31 = Column(np.array(atsys), name='Tsys', description='Tsys before/after OTF')
    #     col32 = Column(np.array(atsmix), name='tsmix', description='mixer of Tsys')
    #     fits.append(ofile, data=Table([col31,col32]).as_array())
    #
    #     col41 = Column(aTa, name='Ta', description='single mixer antenna temperature')
    #     fits.append(ofile, data=Table([col41]).as_array())
    #
    #     # this is a crutch to properly name the extensions!
    #     with fits.open(ofile) as hdu:
    #         hdu[2].header['EXTNAME'] = 'DEBUG1'
    #         hdu[3].header['EXTNAME'] = 'HOTS'
    #         hdu[4].header['EXTNAME'] = 'Tsys'
    #         hdu[5].header['EXTNAME'] = 'Ta'
    #         hdu.writeto(ofile, overwrite=True)

    
    print('saved file: ', ofile)
    # logger.info('saved file: ', ofile)
        
        
    return dfile



def arplsw(iy, lam=1e4, ratio=0.05, itermax=100):
    """
    Inputs:
        y:
            input data (1-dim spectrum)
        lam:
            parameter that can be adjusted by user. The larger lambda is,
            the smoother the resulting background, z
        ratio:
            wheighting deviations: 0 < ratio < 1, smaller values allow less negative values
        itermax:
            number of iterations to perform
    Output:
        z: the fitted background vector 
        
        w: return the weights vector 

    """
    nsel = np.isfinite(iy)
    y = iy[nsel]
    
    N = len(y)
    D = sparse.eye(N, format='csc')
    D = D[1:] - D[:-1]  # numpy.diff( ,2) does not work with sparse matrix. This is a workaround.

    H = lam * D.T * D
    w = np.ones(N)
    for i in range(itermax):
        W = sparse.diags(w, 0, shape=(N, N))
        WH = sparse.csc_matrix(W + H)
        C = sparse.csc_matrix(cholesky(WH.todense()))
        z = spsolve(C, spsolve(C.T, w * y))
        d = y - z
        dn = d[d < 0]
        m = np.mean(dn)
        s = np.std(dn)
        wt = 1. / (1 + np.exp(2 * (d - (2 * s - m)) / s))
        if np.linalg.norm(w - wt) / np.linalg.norm(w) < ratio:
            break
        w = wt
    return z, w

