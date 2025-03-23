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
from astropy.time import Time
from astropy.io import fits
from astropy.table import QTable, Table, Column
from astropy import units as u
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.linalg import cholesky
from scipy.stats import norm
from .GL09PDataIO import loadL08Data
from .GL09PUtils import *
from .GL09PLogger import *
import logging

logger = logging.getLogger('GL09PLogger')


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
        
    # logger
    logger = logging.getLogger('GL09PLogger')
    #logger.addHandler(logging.StreamHandler())
    logger.setLevel(cfi['gprocs']['loglevel'])
    logger.info('Started Level 0.95 Pipeline')
    
    debug = cfi['gprocs']['debug']

    if debug==True:
        print('\nExecuting debug mode.\n')
        # logger = multiprocessing.log_to_stderr()
        # logger.setLevel(multiprocessing.SUBDEBUG)
        n_procs = 1
    else:
        n_procs = multiprocessing.cpu_count() - 2
        
    print('Number of cores used for processing: %i\n'%(n_procs))
    
    # get lines for processing
    print(cfi['gprocs']['lines'], type(cfi['gprocs']['lines']), type([0]), type(cfi['gprocs']['lines'])!=type([0]))
    if type(cfi['gprocs']['lines'])!=type([0]):
        lines = cfi['gprocs']['lines'].replace('[','').replace(']','').split(' ')
    else:
        lines = cfi['gprocs']['lines']
    print('Lines: ', lines)

    # print('Lines: ', lines[0])
    # lines= ['CII', 'NII']
    
    ignore = [0]
    
    
    for line in lines:
        if verbose:
            print('\nProcessing line: ', line)
        # identify the files for processing
        inDir = cfi['gdirs']['L09DataDir']
        outDir = cfi['gdirs']['L095DataDir']
        os.makedirs(outDir, exist_ok=True)
        if line=='NII':
            filter = 'NII*.fits'
        else:
            filter = 'CII*.fits'
        print('outDir: ', outDir)
        print('filter: ', os.path.join(inDir,filter))
        
        # sdirs = sorted(glob.glob(os.path.join(inDir,filter), root_dir=inDir))
        #print(glob.glob(os.path.join(inDir,filter)))
        sdirs = sorted(glob.glob(os.path.join(inDir,filter)))
        #print('single result: ', sdirs[0], os.path.split(sdirs[0]))
        dsc = [int(os.path.split(sdir)[1].split('_')[2].split('.')[0]) for sdir in sdirs]
        
        sdirs.sort(key=lambda sdirs: dsc)
        
        dfiles = []
        for i,ds in enumerate(dsc):
            if (ds >= scanRange[0]) & (ds <= scanRange[1]) & (ds not in ignore):
                dfiles.append(sdirs[i])
                        
        n_ds = len(dfiles)
        if int(cfi['gprocs']['max_files']) > 0:
            n_ds = int(cfi['gprocs']['max_files'])
            dfiles = dfiles[:n_ds]
            
        params = {'line': line, 'inDir': inDir, 'outDir': outDir, 
                  'basemethod': int(cfi['gprocs']['basemethod']),
                  'polyorder': int(cfi['gprocs']['polyorder']),
                  'debug': cfi['gprocs']['debug'], 'verbose': verbose,
                  'loglevel':cfi['gprocs']['loglevel'], 
                  'rowflagfilter': int(cfi['gprocs']['rowflagfilter'])}
        paramlist = [[a, b] for a in dfiles for b in [params]]
        # paramlist = [[a, b, c, e, f, g, d] for a in [line] for b in [inDir] for c in [outDir] for e in [int(cfi['gprocs']['drmethod'])] for f in [debug] for g in [cfi['gprocs']['loglevel']] for d in dfiles]
        # paramlist = [[a, b, c, d, e] for a in [line] for b in [inDir] for c in [outDir] for d in dfiles for e in worker_configurer]
        # pprint(paramlist)
        if verbose:
            print('Number of data files: ', n_ds, len(sdirs))
            #print('Selected data files: ', dfiles)
        
        
        # setup multiprocessing loop here to process each file in list
        with Pool(n_procs) as pool:
            # execute tasks in order
            for result in pool.imap(processL09, paramlist):
                #print(f'Processed: {result}', flush=True)
                logger.info(f'Processed: {result}')
        
    logger.info('Level 0.95 Pipeline finished')
    return n_ds



def processL09(paramlist, verbose=True):
    """Function applying the actual baseline fit


    Parameters
    ----------
    params : list

    Returns
    -------

    """
    import logging
    
    dfile = paramlist[0]
    params = paramlist[1]
    # line, inDir, outDir, dfile, drmethod, debug = params[0], params[1], params[2], params[3], params[4], params[5]
    line, inDir, outDir, drmethod = params['line'], params['inDir'], params['outDir'], params['basemethod']
    debug, verbose, loglevel, rowflagfilter = params['debug'], params['verbose'], params['loglevel'], params['rowflagfilter']
    polyorder = params['polyorder']

    logger = logging.getLogger('GL09PLogger')
    #logger.addHandler(logging.StreamHandler())
    #print('loglevel: ', loglevel)
    logger.setLevel(loglevel)
    if debug:
        logger.info('processing in debug mode')
    #print('debug proc 0.95: ', debug)
    
    # define some processing data first (maybe relocat to function later?)

    if 'ACS3' in dfile:
        # line: CII
        pixel_cut = 600
        band = 2
        add = 'B2'
    else:
        # line: NII
        pixel_cut = 300
        band = 1
        add = 'B1'

    # values needed for baseline correction 
    # will be moved to config file eventually
    bs_lam = 5e4
    bs_lam2 = 1e2
    bs_ratio = 0.02
    bs_itermax = 100
    datavalid = True

    lstr = 'loading file: %s for line: %s'%(os.path.join(inDir,dfile), line)
    logger.info(lstr)
    spec, data, hdr, hdr1 = loadL08Data(os.path.join(inDir,dfile), verbose=False)
    chFlag = data['CHANNEL_FLAG']
    rowFlag = data['ROW_FLAG']
    rflag = checkRowflag(data['ROW_FLAG'], rowflagfilter=rowflagfilter)
    
    n_spec, n_pix = spec.shape
    
    rms = np.zeros(n_spec)
    basecorrf = np.zeros((n_spec,n_pix))
    
    prange = [int(hdr['goodpxst']), int(hdr['goodpxen'])]
    
    # osel = np.argwhere((otfID == data['scanID']) & (otfID.size>=1) & (rfsID.size>2) & (rhsID.size>2) & (hotID.size>2) & (mix == data['MIXER']) & (data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0))
    osel = np.argwhere((data['scan_type'] == 'OTF') & (rflag)).flatten()
    
    if len(osel) <= 0:
        print('File: ', os.path.join(inDir,dfile))
        print(rowFlag.flatten())
        print('WARNING: No OTF spectra available.')
        # logger.warning('No OTF spectra available.')
        datavalid = False
        return 0

    spec_OTF = np.squeeze(spec[osel,:])
    if osel.size==1:
        spec_OTF = spec_OTF[np.newaxis, :]
    n_OTF, n_otfpix = spec_OTF.shape 
    
    basecorr = np.zeros(spec_OTF.shape)
    rmsotf = np.zeros(n_OTF)
    rf = np.zeros(n_OTF)
    
    print('processing data ...')
    # create the calibrated spectra
    for i0 in range(n_OTF):
        # perform a baseline correction on all OTF spectra
        
        # this call might require to address the baseline correction method
        # possible methods are 1); polunomial and 2) 
        #basecorr[i0,:] = baseCorrection()
        # only correct the good part of the spectrum
        # try to mask bad spectral pixels using weights set to zero: TBI
        yywn = spec_OTF[i0,prange[0]:prange[1]]
        chf = chFlag[i0,prange[0]:prange[1]]
        # fsel = np.argwhere(np.isfinite(yywn))
        # nsel = np.argwhere(np.isfinite(yywn) == False)
        # yy = yywn[fsel]
        pargs = np.argwhere(chf == 0).flatten()
        margs = np.argwhere(chf > 0).flatten()

        yy = yywn.copy()
        yy[margs] = np.interp(margs, pargs, yywn[pargs])
        
        if yy.size > 0:
            xx = np.arange(yy.size)
            if basemethod == 0:
                # no baseline fit!
                basecorr[i0,:] = np.zero(basecorr[i0,:].size)
            elif basemethod == 1:
                # polynomial fitter
                
                pass
            elif basemethod == 2:
                baseline_fitter = Baseline(x_data=xx)
                bs, pars = baseline_fitter.aspls(yy, bs_lam)
                rbs, rpars = baseline_fitter.aspls(yy, bs_lam2)
                # rbs, rws = arplsw(yy, lam=bs_lam2, ratio=bs_ratio, itermax=bs_itermax)
                rmsotf[i0] = np.std(yy-rbs)
                
                # bswn = np.zeros(yywn.size)
                # bswn[:] = np.nan
                # bswn[np.squeeze(fsel)] = bs
                basecorr[i0,prange[0]:prange[1]] = bs   #bswn
            else:
                print('Error: Baseline fitting method %i does not exist. \n')
                print('Treated as no baseline fit.')
                basecorr[i0,:] = np.zero(basecorr[i0,:].size)
            
            spec_OTF[i0,prange[0]:prange[1]] -= basecorr[i0,prange[0]:prange[1]]
    
        else:
            # we have to set the rowflag since no valid data
            rf[i0] = 9
            
    # now we have to save the data in a FITS file

    keys = data.dtype.names
    if 'spec' in keys:
        dkey = 'spec'
    elif 'DATA' in keys:
        dkey = 'DATA'
        
    data[dkey][osel,:] = spec_OTF.data
    data['CHANNEL_FLAG'][osel,:] = spec_OTF.mask
    data['ROW_FLAG'][osel] = rf
    rms[osel] = rmsotf
    basecorrf[osel,:] = basecorr
    
    # add the rms column to the data array
    if 'rms' not in data.names:
        rms_col = fits.ColDefs([fits.Column(name='rms', format='E', array=rms)])
        hdu = fits.BinTableHDU.from_columns(data.columns + rms_col)
        data = hdu.data
    else:
        data['rms'] = rms
        
    tred = Time(datetime.datetime.now()).fits
    
    # updating header keywords
    hdr['DLEVEL'] = 0.95
    hdr['PROCTIME'] = tred
    
#    hdr.insert('VLSR', ('PROC_LEV', 0.9, 'pipeline processing level'), after=True)
#    hdr.add_comment('Level 0.9 Pipeline Processing', before='PROC_LEV')
    #hdr.add_comment('Level 0.9 Pipeline Processing', after='VLSR')
    hdr.set('L095PTIM', value=tred, comment=('L0.95 pipeline processing time'), after='RHID2')
    hdr.set('', value='', after='RHID2')
    hdr.set('', value='          Level 0.95 Pipeline Processing', after='RHID2')
    hdr.set('', value='', after='RHID2')
    hdr.set('bs_lam', value=bs_lam, comment='lambda value for ASPLS baseline correction')
    hdr.set('bs_lam2', value=bs_lam2, comment='lambda value for arplsw rms baseline correction')
    hdr.set('bs_ratio', value=bs_ratio, comment='ratio value for arplsw rms baseline correction')
    hdr.set('bs_iterm', value=bs_itermax, comment='max. iterations value for arplsw rms baseline correction')
    hdr.add_comment('L0.95 processing time: %s'%(tred))
    hdr.add_comment('L0.95 version: %s'%(__version__))
    hdr.add_comment('L0.95 last pipeline update: %s'%(__updated__))
    hdr.add_comment('L0.95 developer: %s'%(__author__))
    hdr.set('', value='', after='BS_ITERM')
    
    os.makedirs(outDir, exist_ok=True)
    ofile = os.path.join(outDir, os.path.split(dfile)[1].replace('_L09.fits','_L095.fits'))
    fits.writeto(ofile, data=None, header=hdr, overwrite=True)
    fits.append(ofile, data=data, header=hdr1)
    
    # if in debug mode, add more extensions
    if debug & datavalid:
        col51 = Column(basecorrf, name='basecorr', description='baseline correction')
        col52 = Column(rms, name='rms', description='baseline RMS')
        col53 = Column(spec.data, name='spec', description='spectrum before baseline correction')
        fits.append(ofile, data=Table([col51,col52, col53]).as_array())
        # this is a crutch to properly name the extensions!
        with fits.open(ofile) as hdu:
            hdu[2].header['EXTNAME'] = 'Baseline'
            hdu.writeto(ofile, overwrite=True)

    
    #print('saved file: ', ofile)
    logger.info('saved file: %s'%(ofile))
        
        
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

