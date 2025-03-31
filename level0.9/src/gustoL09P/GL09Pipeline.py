#!/usr/bin/env python
"""
This is the GUSTO L09 Pipeline.
"""

__date__ = '9/19/2024'
__updated__ = '20250210'
__version__ = '0.3.2'
__author__ = 'V. Tolls, CfA | Harvard & Smithsonian'

from joblib import Parallel, delayed
from joblib import Memory
import glob
import numpy as np
import numpy.ma as ma
import time
import argparse
import textwrap
import importlib
import parsl
import sys
import os
import logging
import matplotlib.pyplot as plt
import itertools
import multiprocessing
import astropy.wcs
import zipfile
import lmfit
#import glidertools.cleaning as gc
from matplotlib.patches import Rectangle
from lmfit import Model, Parameters
from pathlib import Path
from pprint import pprint, pformat
from datetime import datetime
from astropy.table import QTable, Table, Column
from astropy import units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame, ICRS
from astropy.time import Time
from astropy.io import fits
from multiprocessing.pool import Pool
from itertools import repeat
from scipy.optimize import minimize
from scipy.signal import savgol_filter

from .GL09PipelineSetupClass import GL09PipelineSetupClass
#from .GL09PConfigClass import GL09PConfigClass, getValues, getRange
from .GL09PProcessControl import GL09PProcessControl
from .GL09PDataIO import loadL08Data
from .GL095Pipeline import GL095Pipeline
from .GL09PUtils import *
from .GL09PLogger import *
from .GL10Pipeline import *
from .GL09PCLArgParser import *

#logger = logging.getLogger(__name__)

spectralLines = ['CII', 'NII']
n_sL = len(spectralLines)

cfg_file0 = importlib.resources.files('gustoL09P') / 'Data/GL09P_setup_pars.txt'
par_file0 = importlib.resources.files('gustoL09P') / 'Data/GUSTO_BaselineData_draft3.xlsx'
pmfileB1 = importlib.resources.files('gustoL09P') / 'Data/ACS5_B1_pixel_masks.txt'
pmfileB2 = importlib.resources.files('gustoL09P') / 'Data/ACS3_B2_pixel_masks.txt'
#cii_file0 = pkg_resources.resource_filename('gustoL09P', 'Data/CIIconfig.txt')

tpipe = 'GUSTO L1 Pipeline'

runtime = time.strftime('%Y%m%d%H%M%S')


def runGL09P(cfi_file=None, verbose=False):
    r"""Function executing the GUSTO Level 2 pipeline.
    The process is:
    1) read the command line parameters
    2) read the configuration parameters
    3) initialize the pipeline
       - analyze which parts of the pipeline are executed
    4) execute the selected pipeline parts
       - calling the appropriate sub-pipeline functions
    4) saving the data if not saved in sub-pipeline

    Parameters
    ----------
    verbose : boolean
        True: plenty of outputs

    Returns
    -------
    None

    Examples
    --------
    on command line execute:
    execGL09P
    """

    # initialize the pipeline
    # this also sets all the directories
    if verbose: 
        print('\n%s: Initializing GUSTO L1 Pipeline'%(time.strftime("%c")))
        print()

    args = GL09PCLArgParser(verbose=verbose)
    if args.configFile is not None:
        cfi_file = args.configFile

    gL09P = GL09PipelineSetupClass()
    status = gL09P.initializePipeline(verbose=verbose, configFile=cfi_file)
    
    #########################################################
    # get the pipeline configuration
    if verbose:
        print('\n%s: Reading pipeline configuration file ...\n'%(time.strftime("%c")))
    cfi = gL09P.getConfigInfo(verbose=verbose)


    cfi = manageArgs(cfi, args, verbose=verbose)

    if verbose:
        print('\nProcessing settings:')
    pprint(cfi)
    print()
    print('debug: ', cfi['gprocs']['debug'], type(cfi['gprocs']['debug']))


    # initialize logging:
    logDir = cfi['gdirs']['logDir']
    os.makedirs(logDir, exist_ok=True)
    
    numeric_level = cfi['gprocs']['loglevel']
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s\nValid options are: DEBUG, INFO, WARNING, ERROR, CRITICAL.' % numeric_level)
    logfile = os.path.join(logDir,'gusto_pipeline_%s.log'%(time.strftime("%Y%m%d%H%M%S")))
    logger = init_logger(loglevel=numeric_level, logfile=logfile, loggername='GL09PLogger')
    logger.addHandler(logging.StreamHandler())
    
    logger.info('Started logging.')
    logger.info('Pipeline configuration file: %s'%(args.configFile))
    logger.info('Pipeline configuration:')
    logger.info(pformat(cfi))


    #########################################################
    # determine the pipeline parts to be executed
    # acceptable ranges are '0.9', '0.95' and '1.0'
    levels = ['0.8', '0.9', '0.95', '1.0']
    isl = args.startLevel
    if isl not in levels:
        print('Please check start data level ...')
    iarg = levels.index(isl)
    elevels = levels[iarg:]
    iel = args.endLevel
    if iel not in levels:
        print('Please check end data level!\nPossible levels are: %s'%(elevels))
    earg = levels.index(iel)
    exlevels = levels[iarg:earg+1]
    print('\nExecuting pipeline levels: ', isl, iel, exlevels)
    print()

    #########################################################
    # determine the pipeline parts to be executed
    # acceptable ranges are '0.8', '0.9', '0.95' and '1.0'
    scanRange = args.scanRange
    print('Scan range for data processing: ', scanRange)
    
    
    print()
    print()
    #########################################################
    # this is the part that forkes to the various pipeline levels    
    
    if '0.9' in exlevels:
        stt = time.time()
        print('\n%s: Executing 0.9 pipeline: calibrating spectra'%(time.strftime("%c")))
        n_pf = GL09Pipeline(cfi, scanRange, verbose=verbose)
        ent = time.time()
        print('Level 0.9 pipeline done.\n')
        print('Execution time: %.2fh  %.2fm   %.2fs'%((ent-stt)/3600.,(ent-stt)/60.,ent-stt))
        print('Execution time: %s'%(time.strftime("%Hh %Mm %Ss",time.gmtime(ent-stt))))
        print('Processed %i files.'%(n_pf))
        print()
        print('The logfile is: %s', logfile)
        print()
    
    if '0.95' in exlevels:
        print('Executing Level 0.95 pipeline: baseline fitter')
        res = GL095Pipeline(cfi, scanRange, verbose=verbose)
        print('Level 0.9 to 0.95 done.\n')
    
    if '1.0' in exlevels:
        print('Executing Level 1.0 pipeline: coordinate corrections')
        res = GL10Pipeline(cfi, scanRange, verbose=verbose)
        print('Level 0.95 to 1.0 done.\n')
    

    print()
    print('Pipeline log file: ', logfile)
    print()
    
    
    #########################################################
    # setup the processing loop for each file
    # the loop should use multiprocessing (parameters like # of cores, should be in setup file)
    # - each file is loaded
    # - determine the data type: OTF or PS to be handled separately
    # - for each file process the refs/hots -> y-factor -> Tsys
    # - for each file process the OTFs
    # - re-evaluate the row flags for each OTF spectrum
    # save spectra to FITS file
        
    



def GL09Pipeline(cfi, scanRange, verbose=False):
    """Function processing the Level 0.8 data. Input are uncalibrated 
    REF, HOT, and OTF spectra and output are calibrated OTF spectra


    Parameters
    ----------
    param1 : int

    Returns
    -------

    """
    #logger = logging.getLogger('GL09PLogger')
    
    if cfi['gprocs']['debug']==True:
        print('\nExecuting debug mode.\n')
        # logger = multiprocessing.log_to_stderr()
        # logger.setLevel(multiprocessing.SUBDEBUG)
        logger = logging.getLogger()
        logger.setLevel(cfi['gprocs']['loglevel'])
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
    # lines= ['CII', 'NII']
    
    # these are the old values for the Oct 24 data sets
    # ignore = [8194, 8338, 9182, 9306, 9314, 9342, 10246, 10250, 10254, 10446, 10606, 
    #           10854, 11006, 11026, 11074, 11102, 11106, 11126, 11134, 24531, 24803, 24881, 26294, 26296]
    
    # these are for the Feb 25 data sets
    ignore = [10086, 13638, 17751, 27083, 28089, 4564, 7165, 7167]
    
    # load ranges for 2nd pixel masking
    
    inDir = cfi['gdirs']['L08DataDir']
    outDir = cfi['gdirs']['L09DataDir']
    os.makedirs(outDir, exist_ok=True)
    
    for line in lines:
        if verbose:
            print('\nProcessing line: ', line)
        # identify the files for processing
        
        # load spikes masking data from files
        if line=='NII':
            filter = 'NII*.fits'
            mranges = readMaskRanges(pmfileB1)
        else:
            filter = 'CII*.fits'
            mranges = readMaskRanges(pmfileB2)
        print('outDir: ', outDir)
        print('filter: ', os.path.join(inDir,filter))
        
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
            
        if line == 'NII':
            vcut = float(cfi['gprocs']['vcut1'])
            pvrange = np.array(getValues(cfi['gprocs']['niivrange']), dtype=int)
            pxrange = np.array(getValues(cfi['gprocs']['niiprange']), dtype=int)
            # pxrange = getRange(cfi['gprocs']['niiprange'], dtype=int)
            # pvrange = getRange(cfi['gprocs']['niivprange'], dtype=int)
        elif line == 'CII':
            vcut = float(cfi['gprocs']['vcut2'])
            pvrange = np.array(getValues(cfi['gprocs']['ciivrange']), dtype=int)
            pxrange = np.array(getValues(cfi['gprocs']['ciiprange']), dtype=int)
            # pxrange = getRange(cfi['gprocs']['ciiprange'], dtype=int)
            # pvrange = getRange(cfi['gprocs']['ciivrange'], dtype=int)
            
        
        params = {'line': line, 'inDir': inDir, 'outDir': outDir, 
                  'drmethod': int(cfi['gprocs']['drmethod']),
                  'debug': cfi['gprocs']['debug'], 'verbose': verbose,
                  'vcut': vcut, 'pxrange': pxrange, 'mranges': mranges,
                  'pvrange': pvrange, 'rowflagfilter': int(cfi['gprocs']['rowflagfilter']),
                  'addpixelflag': cfi['gprocs']['addpixelflag'],
                  'checkringflag': cfi['gprocs']['checkringflag'],
                  'applychannelfilter': cfi['gprocs']['applychannelfilter']}
        paramlist = [[a, b] for a in dfiles for b in [params]]

        if verbose:
            print('Number of data files: ', n_ds, len(sdirs))
            print('Selected data files: ', dfiles)
            print('checkringflag:', cfi['gprocs']['checkringflag'])
        
        
        # setup multiprocessing loop here to process each file in list
        with Pool(n_procs) as pool:
            # execute tasks in order
            for result in pool.imap(processL08, paramlist):
            # for result in pool.imap(processL08, zip(dfiles, repeat(params))):
                print(f'Processed: {result}', flush=True)
        
    return n_ds


def processL08(paramlist):
    """Function processing the Level 0.8 data. Input are uncalibrated 
    REF, HOT, and OTF spectra and output are calibrated OTF spectra


    Parameters
    ----------
    param1 : int

    Returns
    -------

    """
    
    #loadL08Data(dfile, verbose=True)
    dfile = paramlist[0]
    params = paramlist[1]
    # line, inDir, outDir, dfile, drmethod, debug = params[0], params[1], params[2], params[3], params[4], params[5]
    line, inDir, outDir, drmethod, debug = params['line'], params['inDir'], params['outDir'], params['drmethod'], params['debug']
    verbose = params['verbose']
    vcut = params['vcut']
    mranges = params['mranges']
    addpixelflag = params['addpixelflag']
    checkringflag = params['checkringflag']
    rowflagfilter = params['rowflagfilter']
    applychannelfilter = params['applychannelfilter']
    # good pixel ranges
    pxrange = (int(params['pxrange'][0]), int(params['pxrange'][1]))
    # ringing check ranges
    pvrange = (int(params['pvrange'][0]), int(params['pvrange'][1]))

    if debug:
        fadd = '_m%i'%(drmethod)
    
    # define some additional processing data (maybe relocat to function later?)
    if 'CII' in dfile:
        # line: CII
        # all pixels before and after these values are masked as bad
        band = 2
        add = 'B2'
        Tsky = 45  # Kelvin
        rfreq = 1900.5369  # GHz
    else:
        # line: NII
        # all pixels before and after these values are masked as bad
        band = 1
        add = 'B1'
        Tsky = 33.5  # Kelvin
        rfreq = 1461.131406

    datavalid = True

    #logger.info('loading: ', os.path.join(inDir,dfile), ' for line: ', line)
    print('Loading file: ', dfile)
    spec, data, hdr, hdr1 = loadL08Data(os.path.join(inDir,dfile), verbose=False)
    spec_org = spec.copy()
    # we perform a global setting of the data mask to avoid math problems
    # like division by zero below for Tsys and the calibration
    # this defines the global good pixel (or channel) range of the spectra
    # and is indicated by setting bit 8 of the channel_flag
    # the values in these pixels will be set to nan outside the good range and 
    # interpolated within the good range later in the pipeline processing after 
    # the baseline fit.
    # The channel_flag cannot be set in the mask array, but must be set directly:
    # spec.mask[:,pxrange[1]:] = np.bitwise_or(spec.mask[:,pxrange[1]:] ,(1<<7))
    # spec.mask[:,:pxrange[0]+1] = np.bitwise_or(spec.mask[:,:pxrange[0]+1], (1<<7))
    
    data['CHANNEL_FLAG'][:,pxrange[1]:] = np.bitwise_or(data['CHANNEL_FLAG'][:,pxrange[1]:] ,(1<<7))
    data['CHANNEL_FLAG'][:,:pxrange[0]+1] = np.bitwise_or(data['CHANNEL_FLAG'][:,:pxrange[0]+1], (1<<7))
    if applychannelfilter:
        if debug:
            fadd += '_ac'
        spec.mask = data['CHANNEL_FLAG']
        if addpixelflag:
            # more agressive spur flagging
            for i in range(mranges.shape[0]):
                # set the spur flag for the mranges
                data['CHANNEL_FLAG'][:,mranges[i,0]:mranges[i,1]+1] = np.bitwise_or(data['CHANNEL_FLAG'][:,mranges[i,0]:mranges[i,1]+1] ,(1<<7))
                spec.mask = data['CHANNEL_FLAG']
    else:
        spec.mask = np.zeros(spec.shape)
    
    n_spec, n_pix = spec.shape
    
    # for now, process all mixers
    umixers = np.unique(data['MIXER'])
    n_umix = umixers.size
    
    tsyseff_avg = np.zeros(umixers.size)
    aTsyseff = np.zeros([n_spec,n_pix])
    aTsys_mean = np.zeros([n_umix])
    amixer = np.zeros(n_umix)
    ahcorr = np.zeros([n_spec,n_pix])
    ascorr = np.zeros([n_spec,n_pix])
    aspref = np.zeros([n_spec,n_pix])
    aspref2 = np.zeros([n_spec,n_pix])
    ahcoef = np.zeros([n_spec,30])
    afrac = np.zeros([n_spec,2])
    aTam = np.zeros([n_spec])
    aTrms = np.zeros([n_spec])
    asscl = np.zeros([n_spec])
    aTa = np.zeros([n_spec,n_pix])
    aTa2 = np.zeros([n_spec,n_pix])
    if drmethod==3:
        aTsyseff_sm = np.zeros([n_spec,n_pix])
        ahcorr_sm = np.zeros([n_spec,n_pix])
        aspref_sm = np.zeros([n_spec,n_pix])
        aspref2_sm = np.zeros([n_spec,n_pix])
    aghots = []   # hots
    aghtim = []   # unixtime of hots
    aghtint = []  # integration time of hots
    aghmix = []   # mixer of hots
    atsys = []
    ayfac = []
    atsmix = []
    datavalid = np.ones(n_umix, dtype=bool)
    rfsflag = 0

    for k, mix in enumerate(umixers):
        # first check crudely if we have enough data of various scan_types
        msel = np.nonzero(data['MIXER']==mix)
        rflag = checkRowflag(data['ROW_FLAG'], rowflagfilter=rowflagfilter)
        amixer[k] = mix
        
        otfID, rfsID, rhsID, hotID = getSpecScanTypes(mix, spec, data, hdr, rowflagfilter=rowflagfilter)
        check = (np.argwhere(data['scan_type']=='REF').size > 3) & \
                (np.argwhere(data['scan_type']=='HOT').size > 3) & \
                (np.argwhere(data['scan_type']=='REFHOT').size > 3) & \
                (np.argwhere(data['scan_type']=='OTF').size > 5) & \
                (otfID.size>0) & (rfsID.size>0) & (rhsID.size>0) & (hotID.size>0) & np.any(rflag[msel])
                # (otfID.size>0) & (rfsID.size>0) & (rhsID.size>0) & (hotID.size>0) & np.any(data['ROW_FLAG'][msel]==rowflagfilter)
        if not check:
            print('mix, dfile')
            print('check: ', check)
            print('specs: ', spec.shape)
            print('REFs: ', np.argwhere(data['scan_type']=='REF').size, (np.argwhere(data['scan_type']=='REF').size > 3))
            print('HOTs: ', np.argwhere(data['scan_type']=='HOT').size, (np.argwhere(data['scan_type']=='HOT').size > 3))
            print('REFHOTs: ', np.argwhere(data['scan_type']=='REFHOT').size, (np.argwhere(data['scan_type']=='REFHOT').size > 3))
            print('OTFs: ', np.argwhere(data['scan_type']=='OTF').size, (np.argwhere(data['scan_type']=='OTF').size > 5))
            print('other: ', (otfID.size>0), (rfsID.size>0), (rhsID.size>0), (hotID.size>0))
            print('IDs: ', otfID, rfsID, rhsID, hotID)
            print('rowflagfilter: ', rowflagfilter, np.any(checkRowflag(data['ROW_FLAG'][msel], rowflagfilter=rowflagfilter)))
            print('GL09Pipeline: Not enough data available for processing. ROW_FLAGs are set appropriately. ')
            data['ROW_FLAG'][msel] |= 4   # flagged as missing data
            datavalid[k] = False
            return 0
        
        tsys, refs, rhots, rtime, htime, Thot, Tsky, rhIDs, rfsflags, yfac = getCalSpectra(mix, spec, data, hdr, rowflagfilter, verbose=True)
        rfsflag += rfsflags*10**k
        # tsys is a masked array if valid or an int if no good
        if type(tsys)==type(0):
            print('No Tsys available! Stop processing %i of dfile %s'%(mix, dfile), tsys)
            # logger.error('No Tsys available! Stop processing mix of dfile ', mix, dfile, tsys)
            # logger.info('Tsys: ', tsys)
            datavalid[k] = False
            continue
        tsys_mean = np.nanmean(tsys[:,pxrange[0]:pxrange[1]])
        aTsys_mean[k] = tsys_mean
        print('<Tsys_%i>: %.2f'%(mix, tsys_mean))
        tsys.fill_value = 0.0
        #print('tsys shape: ', tsys.shape)
        #print(list(tsys))
        prange = [40, 350]
        #pxs = np.arange(n_pix)
        
        
        otfID, rfsID, rhsID, hotID = getSpecScanTypes(mix, spec, data, hdr, rowflagfilter=rowflagfilter, verbose=verbose)
        # osel = np.argwhere((otfID == data['scanID']) & (otfID.size>=1) & (rfsID.size>2) & (rhsID.size>2) & (hotID.size>2) & (mix == data['MIXER']) & (data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0))
        # 3/8/25: changed this check to exclude the requirement for HOTs since we have sequences without HOT measurements!
        osel = np.argwhere((otfID == data['scanID']) & (rfsID.size>=1) & (rhsID.size>=1) & #(hotID.size>=1) & 
                           (otfID.size>=1) & (mix == data['MIXER']) & (data['scan_type'] == 'OTF') & 
                           (rflag)).flatten()
        #print('otfID.size: ', otfID.size)
        #print('osel: ', len(osel))
        if len(osel) > 2:
            # print('processing OTFs')
            # print('OTFs: ', otfID)
            # print('REFs: ', rfsID)
            # print('REFHOTs: ', rhsID)
            # print('HOTs: ', hotID)
            # print('mixer: ', mix)
            pass
        else:
            print('WARNING: No OTF spectra available for mixer: ', mix)
            # logger.warning('No OTF spectra available.')
            # print('failed processing OTFs')
            # print('scanID: ', (otfID == data['scanID']), '\nrfsID: ', (rfsID.size>=1), '\nrhsID: ', (rhsID.size>=1), '\notfID: ', (otfID.size>=1), 
            #                '\nmix: ', (mix == data['MIXER']), '\nscan type: ', (data['scan_type'] == 'OTF'), '\nrowflag: ', (data['ROW_FLAG'][mix]==0))
            # print('rowflag: ', (data['ROW_FLAG']))
            # print(osel)
            # print('OTFs: ', otfID)
            # print('REFs: ', rfsID)
            # print('REFHOTs: ', rhsID)
            # print('HOTs: ', hotID)
            # print('rfsflag: ', rfsflag)
            # print('mixer: ', mix)
            #print(k, mix, umixers)
            datavalid[k] = False
            data['ROW_FLAG'][msel] |= 1<<19   # flagged as missing data
            continue

        spec_OTF = np.squeeze(spec[osel,:])
        stime = data['UNIXTIME'][osel]
        btime = (rtime[0] + htime[0]) / 2. # before OTFs
        atime = (rtime[1] + htime[1]) / 2. # after OTFs
        if atime > btime:
            fracb = (stime - btime) / (atime - btime)
            fraca = (atime - stime) / (atime - btime)
        else:
            # this is for the case of a single REF/REFHOT pair
            fraca = np.ones(stime.size) / 2.
            fracb = np.ones(stime.size) / 2.
        
        n_OTF, n_opix = spec_OTF.shape
        # antenna temperature is a masked array
        ta = ma.zeros([n_OTF, n_opix])
        ta.mask = spec.mask
        ta2 = ma.zeros([n_OTF, n_opix])
        ta2.mask = spec.mask
        scl = np.ones([n_OTF, n_opix])
        sscl = np.ones(n_OTF)
        tam = np.zeros(n_OTF)
        trms = np.zeros(n_OTF)
        tsyseff = np.zeros([n_OTF, n_opix])
        hcorr = np.zeros([n_OTF, n_opix])
        scorr = np.zeros([n_OTF, n_opix])
        hcoef = np.zeros([n_OTF, 30])
        spref = np.zeros([n_OTF, n_opix])
        spref2 = np.zeros([n_OTF, n_opix])
        if drmethod==3:
            tsyseff_sm = np.zeros([n_OTF, n_opix])
            hcorr_sm = np.zeros([n_OTF, n_opix])
            spref_sm = np.zeros([n_OTF, n_opix])
            spref2_sm = np.zeros([n_OTF, n_opix])

        # this call returns the results for all spectra, but we need only everything for the OTF spectra
        # ahgroup is the assignment of hots to all spectra
        # aghots are the grouped hots
        # aghtim is the unixtime associated with the grouped hots
        # aglast is a flag indicating that there is a hot at the end of the OTF strip
        ahgroup, ghots, ghtim, ghtint, glast, htflag = getHotInfo(spec, data, mix, dfile=dfile, rowflagfilter=rowflagfilter, verbose=True)
        if type(ahgroup)==type(0):
            print('Encountered problem with HOT groups. Flaging rows.')
            data['ROW_FLAG'][msel] |= 1<<19   # flagged as missing data
            continue
        # reduce the assignment to the OTF spectra only
        hgroup = ahgroup[osel]
        gsz = ghots.shape
        n_ghots = gsz[0]
        if n_ghots > 0:
            aghots.append(ghots)
            aghtim.append(ghtim)
            aghtint.append(ghtint)
            aghmix.append(np.zeros(n_ghots)+mix)
        else:
            print('No hots available...')
            data['ROW_FLAG'][msel] |= 1<<19   # flagged as missing data
            continue
        

        atsys.append(tsys)
        ayfac.append(yfac)
        atsmix.append(mix)
        
        if drmethod==3:
            # Define the model function
            def scalefunc(params, spec, sRn):
                model = params['a0']
                nRs, _ = sRn.shape
                for i in range(1,nRs):
                    model += params['a%i'%(i)].value * sRn[i-1,:]
                return ((spec - model)/model)**2
            
            def getsRf(params, sRn):
                model = params['a0']
                nRs, _ = sRn.shape
                for i in range(1,nRs):
                    model += params['a%i'%(i)].value * sRn[i-1,:]
                return model
            
            seq_hots = np.vstack([rhots,ghots])
            n_shots = seq_hots.shape[0]
        
        if drmethod==4:
            def scalefunc(x, c1, c2):
                x1=x[0]
                x2=x[1]
                y = c1 - (c2 * x1 +x2)
                # Return square difference 
                return(np.sum(y*y))
            # rhots and refs should be same number of scans, but might not be
            nhot = rhots.shape[0]
            nref = refs.shape[0]
            nghot = ghots.shape[0]

            # we got the y-factor above
            if yfac.shape[0] == 2:
                y1 = yfac[0,:]
                y2 = yfac[1,:]
            else:
                y1 = yfac[0,:]
                y2 = yfac[0,:]

            # combine all HOTs
            seq_hots = np.vstack([rhots,ghots])

        if drmethod==5:
            # Define the model function
            def scalefunc(params, spec, sRn):
                model = params['a0']
                nRs, _ = sRn.shape
                for i in range(1,nRs):
                    model += params['a%i'%(i)].value * sRn[i-1,:]
                return ((spec - model)/model)**2
            
            def getsRf(params, sRn):
                model = params['a0']
                nRs, _ = sRn.shape
                for i in range(1,nRs):
                    model += params['a%i'%(i)].value * sRn[i-1,:]
                return model
            
            # seq_hots = np.vstack([rhots,ghots])
            # n_shots = seq_hots.shape[0]
        
        # create the calibrated spectra
        for i0 in range(n_OTF):
            tsyseff[i0,:] = fracb[i0] * tsys[0,:] + fraca[i0] * tsys[1,:]
            # we might have to replace flagged pixels in tsyseff to not cause a problem in the spectra
            # => skipped for now since pixel flags should be very similar to flagged pixels in spectra
            if drmethod==1:
                # method 1
                spref[i0,:] = fracb[i0] * refs[0,:] + fraca[i0] * refs[1,:]
                ta[i0,:] = 2.*tsyseff[i0,:] * (spec_OTF[i0,:] - spref[i0,:])/spref[i0,:]
            elif drmethod==2:
                # method 2: using HOTS to mitigate drifts
                # apply the REFHOTS to the refs
                # ToDo: check if the inttime of refs/hots/otfs matters
                # calculate the fraction of hot used for the spectrum
                ht1 = ghtim[hgroup[i0]]
                if hgroup[i0]+1 <= hgroup.max():
                    ht2 = ghtim[hgroup[i0]+1]
                else:
                    ht2 = ghtim[hgroup[i0]]
                # ht2 = ghtim[hgroup[i0]+1]
                if ht2 > ht1:
                    hfrac = (stime[i0]-ht1)/(ht2-ht1)
                else:
                    hfrac = 1.0
                # determine the hots for the individual OTF spectra
                if hgroup[i0]+1 <= hgroup.max():
                    hcorr[i0,:] = ghots[hgroup[i0],:]*hfrac + (1-hfrac) * ghots[hgroup[i0]+1,:]
                else:
                    hcorr[i0,:] = ghots[hgroup[i0],:]*hfrac + (1-hfrac) * ghots[hgroup[i0],:]
                #hcorr_sm[i0,:] = smoothSpectrum(hcorr[i0,:], cflag, window_length=5, polyorder=2)
                # hcorr[i0,:] = ghots[hgroup[i0],:]*hfrac + (1-hfrac) * ghots[hgroup[i0]+1,:]
                # determine the hots-reduced REF spectra
                # spref[i0,:] = (fracb[i0] * refs[0,:] / ghots[0,:] + fraca[i0] * refs[1,:] / ghots[-1,:]) * hcorr[i0,:]
                scorr[i0,:] = hcorr[i0,:] / (fracb[i0] * rhots[0,:] + fraca[i0] * rhots[1,:])
                spref[i0,:] = (fracb[i0] * refs[0,:] + fraca[i0] * refs[1,:]) * scorr[i0,:]
                spref2[i0,:] = fracb[i0] * refs[0,:] + fraca[i0] * refs[1,:]
                #spref_sm[i0,:] = smoothSpectrum(spref_sm[i0,:], cflag, window_length=5, polyorder=2)
                
                # calculate a scaling factor
                scl[i0,:] = spec_OTF[i0,:] / (spref[i0,:])
                sscl[i0] = ma.median(scl[i0,:])
                #print('Scaling factor: %.4f'%(sscl[i0]))
                
                # put everything together. issue: divide by zero -> catch in masks
                ta[i0,:] = 2.*tsyseff[i0,:] * (spec_OTF[i0,:] - spref[i0,:])/spref[i0,:]
                ta2[i0,:] = 2.*tsyseff[i0,:] * (spec_OTF[i0,:] - spref[i0,:])/spref2[i0,:]
                
                ta[i0,:] = ta[i0,:] - ma.median(ta[i0,:])
                ta2[i0,:] = ta2[i0,:] - ma.median(ta2[i0,:])
                
                # calculate an average Ta
                if np.any(np.isfinite(ta[i0,pxrange[0]:pxrange[1]])):
                    tam[i0] = np.nanmean(ta[i0,pxrange[0]:pxrange[1]])
                else:
                    tam[i0] = np.nan
                # print('Ta mean (%i): %.4f'%(i0, Ta_mean), hgroup[i0], hfrac, fraca[i0], fracb[i0], stime[i0], ht1, ht2)
            elif drmethod==3:
                
                sspec = spec_OTF[i0,:]
                
                yvalid = np.nonzero((yfac[1,:].squeeze() > 1.0))[0]
                sRn = seq_hots / yfac[1,:].squeeze()
                try:
                    # Create parameters
                    params = lmfit.Parameters()
                    vals = np.zeros(n_shots) + 0.5
                    vals[0] = 0.0
                    for index, value in enumerate(vals):
                        params.add('a%i'%index, value=vals[index])
                    
                    # Minimize the objective function
                    mini = lmfit.Minimizer(scalefunc, params, fcn_args=(sspec[yvalid], sRn[:,yvalid]))
                    result = mini.minimize()
                except:
                    print('Problems with REF fitting for spectrum %i in %s ...'%(i0, dfile), flush=True)
                    data['ROW_FLAG'][i0] |= 1<<24   # flagged as failed fit
                    continue
                    

                sR = getsRf(result.params, sRn)
                spref[i0,:] = sR 
                hcoef[i0,:n_shots] = np.array(list(result.params.valuesdict().values()))

                ta[i0,:] = 2.*tsyseff[i0,:] * (sspec-sR)/sR
                trms[i0] = np.std(ta[i0,yvalid])

            elif drmethod==4:
                
                cflag = data['CHANNEL_FLAG'][i0,:]
                #  S-R / R implementation fitting R to S
                spec_otf1 = spec_OTF[i0,:]

                # Useable Y factors should be > 1.0
                quse = np.argwhere((y1 > 1.0) & (cflag == 0))
                nsamp = quse.shape[0]
                #print("Good channels",nsamp)
                if nsamp == 0:
                    #Try 2nd Y factor y2
                    print('switching to 2nd Y')
                    y1 = y2
                    quse = np.argwhere((y2 > 1.0) & (cflag == 0))
                    nsamp = quse.shape[0]

                if nsamp > 0:

                    testvar = []
                    mincoeffs = []
                    synRefs=[]
                    for hot1 in seq_hots:
                        sR1 = hot1/y1
                        xstart = [1.,0.0]
                        synRefs.append(sR1)
                        minresult = minimize(scalefunc,xstart,args=(spec_otf1[quse],sR1[quse]))
                        mincoeffs.append(minresult.x)
                        testvar.append(minresult.fun/nsamp)
                        #print(testvar, minresult.x)
                    testvar = np.array(testvar)
                    mincoeffs = np.array(mincoeffs)
                    #find the minimum variance
                    qmin = np.argmin(testvar)
                    #print(mincoeffs.shape)
                    #print(testvar[qmin],mincoeffs[qmin])
                    synRef = mincoeffs[qmin][0]*synRefs[qmin] + mincoeffs[qmin][1]
                    
                    Ta_otf = 2. * tsyseff[i0,:] * ( spec_otf1 - synRef ) / synRef
                    # Ta_otf = 2. * Tsys_nominal * ( spec_otf1 - synRef ) / synRef
                    ta[i0,:] = Ta_otf
                    trms[i0] = np.std(Ta_otf[quse])
                    # this line is wrong!
                    BW = hdr['CDELT1']*1e6 *2 
                    #Expected = np.median(tsyseff1[quse])/np.sqrt(BW*0.33)
                    Expected = 1100.0/np.sqrt(BW*0.33)
                    #print(Trms / Expected)
            elif drmethod==5:
                
                sspec = spec_OTF[i0,:]
                
                seq_hots = ghots[hgroup[i0],:]
                if hgroup[i0]+1 <= hgroup.max():
                    seq_hots = np.vstack([seq_hots, ghots[hgroup[i0]+1,:]])
                if hgroup[i0]-1 >= 1:
                    seq_hots = np.vstack([seq_hots, ghots[hgroup[i0]-1,:]])
                n_shots = seq_hots.shape[0]
                
                yvalid = np.nonzero((yfac[1,:].squeeze() > 1.0))[0]
                sRn = seq_hots / yfac[1,:].squeeze()
                try:
                    # Create parameters
                    params = lmfit.Parameters()
                    vals = np.zeros(n_shots) + 0.5
                    vals[0] = 0.0
                    for index, value in enumerate(vals):
                        params.add('a%i'%index, value=vals[index])
                    
                    # Minimize the objective function
                    mini = lmfit.Minimizer(scalefunc, params, fcn_args=(sspec[yvalid], sRn[:,yvalid]))
                    result = mini.minimize()
                except:
                    print('Problems with REF fitting for spectrum %i in %s ...'%(i0, dfile), flush=True)
                    data['ROW_FLAG'][i0] |= 1<<24   # flagged as failed fit
                    continue
                    

                sR = getsRf(result.params, sRn)
                spref[i0,:] = sR 
                hcoef[i0,:n_shots] = np.array(list(result.params.valuesdict().values()))

                ta[i0,:] = 2.*tsyseff[i0,:] * (sspec-sR)/sR
                trms[i0] = np.std(ta[i0,yvalid])



        # now we have to save the data in a FITS file
        # Note: changing the output data format to
        # primarily passing through the header and the data table
        # and only few variables and the header keys are updated or added
        
        tsyseff_avg[k] = np.nanmean(tsys[:,int(pxrange[0]):int(pxrange[1])+1])
        aTsyseff[osel,:] = tsyseff
        ahcorr[osel,:] = hcorr
        ascorr[osel,:] = scorr
        aspref[osel,:] = spref
        afrac[osel,0] = fraca
        afrac[osel,1] = fracb
        ahcoef[osel,:] = hcoef
        aTa[osel,:] = ta
        aTam[osel] = tam
        aTrms[osel] = trms
        asscl[osel] = sscl
        aspref2[osel,:] = spref2
        aTa2[osel,:] = ta2
        if drmethod==3:
            ahcorr_sm[osel,:] = hcorr_sm
            aTsyseff_sm[osel,:] = tsyseff_sm
            aspref_sm[osel,:] = spref_sm
            aspref2_sm[osel,:] = spref2_sm
    
        keys = data.dtype.names
        if 'spec' in keys:
            dkey = 'spec'
        elif 'DATA' in keys:
            dkey = 'DATA'
            
        data[dkey][osel,:] = ta.data
        #data['CHANNEL_FLAG'] [osel,:] = ta.mask

    # analize datavalid
    validdata = np.any(datavalid)
    #print('validdata: ', datavalid)
    
    if validdata == False:
        print('Not enough data available for saving. ')
        return 0

    # check if there is ringing in the uncalibrated spectra
    var = np.zeros(n_spec)
    for i in range(0,n_spec):
        if (data['scan_type'][i] == 'OTF'):
            if (np.any(np.isfinite(spec[i,pvrange[0]:pvrange[1]+1]))):
                tsp = spec[i,pvrange[0]:pvrange[1]+1]
                xx = np.arange(len(tsp))
                # remove the linear baseline
                # fit a line to the data
                p = np.polyfit(xx, tsp, 1)
                # remove the line
                tsp = tsp - np.polyval(p, xx)
                var[i] = np.abs((np.nanmax(tsp) - np.nanmin(tsp)) / np.nanmean(tsp))
                if (np.abs(var[i]) >= vcut):
                    if checkringflag:
                        data['ROW_FLAG'][i] |= (1 << 14) # set bit 14
            else:
                var[i] = 9999
                if checkringflag:
                    data['ROW_FLAG'][i] |= (1 << 13)   # set bit 13
        else:
            var[i] = 0.0
        # if debug:
        #     print('%4i  %7.3f  %2i  %6s  %i  %i'%(i, var[i], data['MIXER'][i], data['scan_type'][i], data['row_flag'][i], data['ROW_FLAG'][i]))

    
    tred = Time(datetime.datetime.now()).fits
    
    # updating header keywords
    hdr.set('DLEVEL', value = 0.9, after='PROCTIME')
    hdr['PROCTIME'] = tred
    
#    hdr.insert('VLSR', ('PROC_LEV', 0.9, 'pipeline processing level'), after=True)
#    hdr.add_comment('Level 0.9 Pipeline Processing', before='PROC_LEV')
    #hdr.add_comment('Level 0.9 Pipeline Processing', after='VLSR')
    hdr.set('L09PTIME', value=tred, comment=('L0.9 pipeline processing time'))
    hdr.set('goodpxst', value=pxrange[0], comment='pixel index where good pixels start')
    hdr.set('goodpxen', value=pxrange[1], comment='pixel index where good pixel stop')
    hdr.set('ringpxst', value=pvrange[0], comment='start pixel index for ringing test')
    hdr.set('ringpxen', value=pvrange[1], comment='end pixel index for ringing test')
    hdr.set('rwflfilt', value=rowflagfilter, comment='applied rowflag filter for useful spectra')
    hdr.set('RFSFLAG', value=rfsflag, comment='flag indicating which REFS are available')
    hdr.set('rhID1', value=rhIDs[0], comment='scan ID for REFHOT/REF before OTF')
    hdr.set('rhID2', value=rhIDs[1], comment='scan ID for REFHOT/REF after OTF')
    for l in range(n_umix):
        hdr.set('Tsysmean%i'%(l), value=atsys_mean[l], comment='mean Tsys for mixer %i over goodpx'%(l))
        hdr.set('MixerTsys%i'%(l), value=atsmix[l], comment='mixer for Tsys_mean%i'%(l))
    hdr.set('procmeth', value=drmethod, comment='data reduction processing method applied')
    hdr.set('', value='')
    hdr.set('', value='', after='VLSR')
    hdr.set('', value='          Level 0.9 Pipeline Processing', after='VLSR')
    hdr.set('', value='', after='VLSR')
    hdr.add_comment('L0.9 processing time: %s'%(tred))
    hdr.add_comment('L0.9 version: %s'%(__version__))
    hdr.add_comment('L0.9 last pipeline update: %s'%(__updated__))
    hdr.add_comment('L0.9 developer: %s'%(__author__))
    hdr.set('', value='', after='procmeth')
    
    os.makedirs(outDir, exist_ok=True)
    if debug:
        ofile = os.path.join(outDir+'_m%i'%(drmethod), os.path.split(dfile)[1].replace('.fits','%s_L09.fits'%(fadd)))
        os.makedirs(outDir+'_m%i'%(drmethod), exist_ok=True)
    else:
        ofile = os.path.join(outDir, os.path.split(dfile)[1].replace('.fits','_L09.fits'))
    fits.writeto(ofile, data=None, header=hdr, overwrite=True)
    fits.append(ofile, data=data, header=hdr1)
    
    # if in debug mode, add more extensions
    if debug & validdata:
        # data needed
        # Tsys1, Tsys2, fraca, fracb, Tsyseff, hcorr, spref
        # we have spectra/data for each otf spectrum: Tsyseff, hcorr, spref, fraca, fracb
        # single spectra: Tsys1, Tsys2, 
        # small groupf of spectra: hots
        col1 = Column(aTsyseff, name='Tsyseff', description='effective Tsys per OTF spectrum')
        col2 = Column(ahcorr, name='hcorr', description='effective hot spectrum')
        col3 = Column(aspref, name='spref', description='effective hot for ref spectrum')
        col4 = Column(aspref2, name='spref2', description='ref spectrum')
        col5 = Column(ascorr, name='scorr', description='effective hot spectrum correction')
        if drmethod==3:
            col10 = Column(ahcorr_sm, name='hcorr_sm', description='effective smoothed hot spectrum')
            col11 = Column(aspref_sm, name='spref_sm', description='effective smoothed hot for ref spectrum')
            col12 = Column(aspref2_sm, name='spref2_sm', description='smoothed ref spectrum')
            col13 = Column(aTsyseff_sm, name='Tsyseff_sm', description='smoothed effective Tsys per OTF spectrum')
        col6 = Column(afrac, name='frac', description='fraction of Tsys1/2')
        col7 = Column(spec_org, name='spec', description='original spectra')
        col8 = Column(aTa2, name='tant', description='spectrum without hot correction')
        col9 = Column(var, name='ringvar', description='value for determining if ringing in spectrum')
        col14 = Column(aTam, name='Ta_mean', description='mean of Ta over good pixel range')
        col15 = Column(ahcoef, name='hot coefs', description='coef for hots to create synth. REF')
        
        if drmethod==3:
            fT = Table([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15])
        else:
            fT = Table([col1, col2, col3, col4, col5, col6, col7, col8, col9, col14, col15])
        hdre1 = fits.Header()
        hdre1.set('vcut', value=vcut, comment='ringing treshold used with ringvar')
        fits.append(ofile, data=fT.as_array())
        
        #print('aghots: ', len(aghots))
        if len(aghots) > 0:
            col21 = Column(np.vstack(aghots), name='hots', description='averaged hot spectra')
            col22 = Column(np.hstack(aghtim), name='htime', description='utime of avg. hots')
            col23 = Column(np.hstack(aghmix), name='hmixer', description='mixer of avg. hots')
            col24 = Column(np.hstack(aghtint), name='htint', description='integration time of avg. hots')
            fgh = Table([col21, col22, col23, col24])
            fgh.write(ofile, append=True)
        else:
            fgh = Table()
            fgh.write(ofile, append=True)

        col31 = Column(np.array(atsys), name='Tsys', description='Tsys before/after OTF')
        col32 = Column(np.array(atsmix), name='tsmix', description='mixer of Tsys')
        col33 = Column(np.array(ayfac), name='yfactor', description='y-factor')
        fits.append(ofile, data=Table([col31, col32, col33]).as_array())

        col41 = Column(aTa, name='Ta', description='single mixer antenna temperature')
        fits.append(ofile, data=Table([col41]).as_array())
        
        # this is a crutch to properly name the extensions!
        with fits.open(ofile) as hdu:
            hdu[2].header['EXTNAME'] = 'DEBUG1'
            hdu[3].header['EXTNAME'] = 'HOTS'
            hdu[4].header['EXTNAME'] = 'Tsys'
            hdu[5].header['EXTNAME'] = 'Ta'
            hdu.writeto(ofile, overwrite=True)

        
    
    print('saved file: ', ofile)
    # logger.info('saved file: ', ofile)
        
        
    return dfile



    
    
