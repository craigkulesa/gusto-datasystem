#!/usr/bin/env python
"""
This is the GUSTO L09 Pipeline.
"""

__date__ = '9/19/2024'
__updated__ = '20241206'
__version__ = '0.3'
__author__ = 'V. Tolls, CfA | Harvard & Smithsonian'

from joblib import Parallel, delayed
from joblib import Memory
import glob
import numpy as np
import numpy.ma as ma
import time
import argparse
import textwrap
import pkg_resources
import parsl
import sys
import os
import logging
import matplotlib.pyplot as plt
import itertools
import multiprocessing
import astropy.wcs
from matplotlib.patches import Rectangle
from pathlib import Path
from pprint import pprint
from datetime import datetime
from astropy.table import QTable, Table, Column
from astropy import units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame, ICRS
from astropy.time import Time
from astropy.io import fits
from multiprocessing.pool import Pool

from .GL09PipelineSetupClass import GL09PipelineSetupClass
from .GL09PConfigClass import GL09PConfigClass, getValues, getRange
from .GL09PProcessControl import GL09PProcessControl
from .GL09PDataIO import loadL08Data
from .GL095Pipeline import GL095Pipeline
from .GL09PUtils import *
from .GL09PLogger import *


spectralLines = ['CII', 'NII', 'OI']
n_sL = len(spectralLines)

cfg_file0 = pkg_resources.resource_filename('gustoL09P', 'Data/GL09P_setup_pars.txt')
par_file0 = pkg_resources.resource_filename('gustoL09P', 'Data/GUSTO_BaselineData_draft3.xlsx')
#cii_file0 = pkg_resources.resource_filename('gustoL09P', 'Data/CIIconfig.txt')

tpipe = 'GUSTO L1 Pipeline'

runtime = time.strftime('%Y%m%d%H%M%S')

def runGL09P(verbose=False):
    r"""Function running the GUSTO Level 2 pipeline.
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
    parser = argparse.ArgumentParser(
        prog='gustoL1',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
            GUSTO Level 1 Pipeline CLI
            --------------------------------
                I have indented it
                exactly the way
                I want it
            '''))
    parser.add_argument('--configFile', '-c', nargs='?', type=str,
                        default=cfg_file0,
                        help='GUSTO L1 Pipeline configuration file. This file contains the data directories, etc.')
    parser.add_argument('--startLevel', '-s', nargs='?', type=str,
                        default='0.8',
                        help='GUSTO Data Level the pipeline should start processing (default=0.8); \npossible entries are 0.8, 0.9 and 0.95')
    parser.add_argument('--endLevel', '-e', nargs='?', type=str,
                        default='0.9',
                        help='GUSTO Data Level produced by pipeline (default=1.0); possible entries are 0.9 and 1.0')
    parser.add_argument('--scanRange', '-r', nargs=2, type=int,
                        default=[2000, 30000],
                        help='Range of scans to be processed by pipeline.')
    parser.add_argument('--loglevel', '-l', type=str,
                        default='INFO',
                        help='setting the log level of the {tpipe}')
    parser.add_argument('--verbose', '-v', type=bool,
                        default=False,
                        help='enables verbosity of the {tpipe}')
    parser.add_argument('--debug', '-d', type=bool,
                        default=False,
                        help='enables settings for debugging pipeline')
    args = parser.parse_args()
    if verbose:
        print('commandline arguments:\n', args, '\n')
        print('configFile: ', args.configFile)
        print()

    # this overrides the verbosity from above in case it is enabled
    if args.verbose:
        verbose = args.verbose

    # inspect any provided arguments

    # initialize the pipeline
    # this also sets all the directories
    if verbose: 
        print('\n%s: Initializing GUSTO L1 Pipeline'%(time.strftime("%c")))
        print()

    gL09P = GL09PipelineSetupClass()
    status = gL09P.initializePipeline(verbose=verbose, configFile=args.configFile)
    
    #########################################################
    # get the pipeline configuration
    if verbose:
        print('\n%s: Reading pipeline configuration file ...'%(time.strftime("%c")))
    cfi = gL09P.getConfigInfo(verbose=verbose)
    if verbose:
        print('\nProcessing settings:')
        #pprint(cfi)
        
    
    # initialize logging:
    logDir = cfi['gdirs']['logDir']
    numeric_level = getattr(logging, args.loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s\nValid options are: DEBUG, INFO, WARNING, ERROR, CRITICAL.' % loglevel)
    logfile = os.path.join(logDir,'gusto_pipeline_%s.log'%(time.strftime("%Y%m%d%H%M%S")))
    logger = init_logger(loglevel=numeric_level, logfile=None, loggername='GL09PLogger')
    
    logger.info('Started logging.')
    logger.info('Pipeline configuration file: %s'%(args.configFile))
    logger.info('Pipeline configuration:')
    logger.info(cfi)

    
    if args.debug:
        cfi['gprocs']['debug'] = args.debug
        

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
    print('levels: ', isl, iel, exlevels)
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
    
    if '1' in exlevels:
        print('Executing Level 1.0 pipeline: coordinate corrections')
        #res = GL10Pipeline(cfi, scanRange, verbose=verbose)
        print('Level 0.95 to 1.0 done.\n')
    

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
    
    
    #########################################################
    # determine the range of data files to be processed
    
    



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
    
    ignore = [8194, 8338, 9182, 9306, 9314, 9342, 10246, 10250, 10254, 10446, 10606, 
              10854, 11006, 11026, 11074, 11102, 11106, 11126, 11134, 24803, 24881, 26294, 26296]
    
    
    for line in lines:
        if verbose:
            print('Processing line: ', line)
        # identify the files for processing
        inDir = cfi['gdirs']['L08DataDir']
        outDir = cfi['gdirs']['L09DataDir']
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
            print('Selected data files: ', dfiles)
        
        
        # setup multiprocessing loop here to process each file in list
        with Pool(n_procs) as pool:
            # execute tasks in order
            for result in pool.imap(processL08, paramlist):
                print(f'Processed: {result}', flush=True)
        
    return n_ds


def processL08(params, verbose=False):
    """Function processing the Level 0.8 data. Input are uncalibrated 
    REF, HOT, and OTF spectra and output are calibrated OTF spectra


    Parameters
    ----------
    param1 : int

    Returns
    -------

    """
    
    #loadL08Data(dfile, verbose=True)
    line, inDir, outDir, dfile, drmethod, debug = params[0], params[1], params[2], params[3], params[4], params[5]
    
    # define some processing data first (maybe relocat to function later?)

    if 'ACS3' in dfile:
        pfr_ra = np.array(list([[-0.01,0.01],[0.25,0.38],[2.25,2.33]]), dtype=float)
        pp_ra = np.array([[20,83],[132,138],[200,205],[265,274],[323,340],[399,407],[498,509]], dtype=int)
        # all pixels before and after these values are masked as bad
        pixel_st = 80
        pixel_en = 600
        band = 2
        add = 'B2'
        Tsky = 45  # Kelvin
        rfreq = 1900.5369  # GHz
    else:
        pfr_ra = np.array(list([[-0.01,0.01],[0.33,0.37],[2.26,2.31],[3.92,3.96]]), dtype=float)
        pp_ra = np.array([[23,50],[65,71],[100,103],[163,170],[198,206],[300, 511]], dtype=int)
        # all pixels above this value are masked as bad
        # all pixels before and after these values are masked as bad
        pixel_st = 80
        pixel_en = 300
        band = 1
        add = 'B1'
        Tsky = 33.5  # Kelvin
        rfreq = 1461.131406

    datavalid = True

    #logger.info('loading: ', os.path.join(inDir,dfile), ' for line: ', line)
    print('Loading file: ', dfile)
    spec, data, hdr, hdr1 = loadL08Data(os.path.join(inDir,dfile), verbose=False)
    # we perform a global setting of the data mask to avoid math problems
    # like division by zero below for Tsys and the calibration
    # this defines the global good pixel (or channel) range of the spectra
    # the values in these pixels will be set to nan outside the good range and 
    # interpolated within the good range later in the pipeline processing after 
    # the baseline fit.
    spec.mask[:,pixel_en:] = 1
    spec.mask[:,:pixel_st] = 1
    
    rowFlag = data['ROW_FLAG']
    
    n_spec, n_pix = spec.shape
    
    # for now, process all mixers
    umixers = np.unique(data['MIXER'])
    
    tsyseff_avg = np.zeros(umixers.size)
    aTsyseff = np.zeros([n_spec,n_pix])
    ahcorr = np.zeros([n_spec,n_pix])
    aspref = np.zeros([n_spec,n_pix])
    aspref2 = np.zeros([n_spec,n_pix])
    afrac = np.zeros([n_spec,2])
    aTa = np.zeros([n_spec,n_pix])
    aTa2 = np.zeros([n_spec,n_pix])
    aghots = []   # hots
    aghtim = []   # unixtime of hots
    aghtint = []  # integration time of hots
    aghmix = []   # mixer of hots
    atsys = []
    atsmix = []

    for k, mix in enumerate(umixers):
        # first check crudely if we have enough data of various scan_types
        msel = np.argwhere(data['MIXER']==mix)
        otfID, rfsID, rhsID, hotID = getSpecScanTypes(mix, spec, data, hdr)
        check = (np.argwhere(data['scan_type']=='REF').size > 3) & \
                (np.argwhere(data['scan_type']=='HOT').size > 3) & \
                (np.argwhere(data['scan_type']=='REFHOT').size > 3) & \
                (np.argwhere(data['scan_type']=='OTF').size > 5) & \
                (otfID.size>0) & (rfsID.size>0) & (rhsID.size>0) & (hotID.size>0) & np.any(rowFlag[msel]==0)
        if not check:
            print('mix, dfile')
            print('specs: ', spec.shape)
            print('REFs: ', np.argwhere(data['scan_type']=='REF').size)
            print('HOTs: ', np.argwhere(data['scan_type']=='HOT').size)
            print('REFHOTs: ', np.argwhere(data['scan_type']=='REFHOT').size)
            print('OTFs: ', np.argwhere(data['scan_type']=='OTF').size)
            print('Not enough data available for processing. ROW_FLAGs are set appropriately. ')
            data['ROW_FLAG'][msel] = 4   # flagged as missing data
            datavalid = False
            return 0
        
        tsys, refs, rhots, rtime, htime, Thot, Tsky = getCalSpectra(mix, spec, data, hdr, verbose=True)
        # tsys is a masked array if valid or an int if no good
        if type(tsys)==type(0):
            print('No Tsys available! Stop processing mix of dfile ', mix, dfile, tsys)
            # logger.error('No Tsys available! Stop processing mix of dfile ', mix, dfile, tsys)
            # logger.info('Tsys: ', tsys)
            datavalid = False
            break
        #print('<Tsys>: ', np.nanmean(tsys))
        tsys.fill_value = 0.0
        #print('tsys shape: ', tsys.shape)
        #print(list(tsys))
        prange = [40, 350]
        #pxs = np.arange(n_pix)
        
        
        otfID, rfsID, rhsID, hotID = getSpecScanTypes(mix, spec, data, hdr, verbose=verbose)
        
        # osel = np.argwhere((otfID == data['scanID']) & (otfID.size>=1) & (rfsID.size>2) & (rhsID.size>2) & (hotID.size>2) & (mix == data['MIXER']) & (data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0))
        osel = np.argwhere((otfID == data['scanID']) & (rfsID.size>=1) & (rhsID.size>=1) & (hotID.size>=1) & 
                           (otfID.size>=1) & (mix == data['MIXER']) & (data['scan_type'] == 'OTF') & 
                           (data['ROW_FLAG'][mix]==0)).flatten()
        #print('otfID.size: ', otfID.size)
        if len(osel) > 0:
            # print('processing OTFs')
            # print('OTFs: ', otfID)
            # print('REFs: ', rfsID)
            # print('REFHOTs: ', rhsID)
            # print('HOTs: ', hotID)
            # print('mixer: ', mix)
            pass
        else:
            print('WARNING: No OTF spectra available.')
            # logger.warning('No OTF spectra available.')
            datavalid = False
            return 0
    
        spec_OTF = np.squeeze(spec[osel,:])
        stime = data['UNIXTIME'][osel]
        btime = (rtime[0] + htime[0]) / 2. # before OTFs
        atime = (rtime[1] + htime[1]) / 2. # after OTFs
        fracb = (stime - btime) / (atime - btime)
        fraca = (atime - stime) / (atime - btime)
        
        n_OTF, n_opix = spec_OTF.shape
        # antenna temperature is a masked array
        ta = ma.zeros([n_OTF, n_opix])
        ta.mask = spec.mask
        ta2 = ma.zeros([n_OTF, n_opix])
        ta2.mask = spec.mask
        tsyseff = np.zeros([n_OTF, n_opix])
        hcorr = np.zeros([n_OTF, n_opix])
        spref = np.zeros([n_OTF, n_opix])
        spref2 = np.zeros([n_OTF, n_opix])

        # this call returns the results for all spectra, but we need only everything for the OTF spectra
        # ahgroup is the assignment of hots to all spectra
        # aghots are the grouped hots
        # aghtim is the unixtime associated with the grouped hots
        # aglast is a flag indicating that there is a hot at the end of the OTF strip
        ahgroup, ghots, ghtim, ghtint, glast, htflag = getHotInfo(spec, data, mix, dfile=dfile, verbose=True)
        # reduce the assignment to the OTF spectra only
        hgroup = ahgroup[osel]
        gsz = ghots.shape
        n_ghots = gsz[0]
        if n_ghots > 0:
            aghots.append(ghots)
            aghtim.append(ghtim)
            aghtint.append(ghtint)
            aghmix.append(np.zeros(n_ghots)+mix)
        
        atsys.append(tsys)
        atsmix.append(mix)
        
    
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
                ht2 = ghtim[hgroup[i0]+1]
                hfrac = (stime[i0]-ht1)/(ht2-ht1)
                # determine the hots for the individual OTF spectra
                hcorr[i0,:] = ghots[hgroup[i0],:]*hfrac + (1-hfrac) * ghots[hgroup[i0]+1,:]
                # determine the hots-reduced REF spectra
                spref[i0,:] = fracb[i0] * refs[0,:] / ghots[0,:] + fraca[i0] * refs[1,:] / ghots[-1,:]
                spref2[i0,:] = fracb[i0] * refs[0,:] + fraca[i0] * refs[1,:]
                
                # put everything together. issue: divide by zero -> catch in masks
                ta[i0,:] = 2.*tsyseff[i0,:] * (spec_OTF[i0,:]/hcorr[i0,:] - spref[i0,:])/spref[i0,:]
                ta2[i0,:] = 2.*tsyseff[i0,:] * (spec_OTF[i0,:] - spref[i0,:])/spref2[i0,:]
                #print('%4i %i %7.2f %7.2f %7.2f %7.2f '%(i0, mix, np.nanmin(ta[i0,200:400]), np.nanmax(ta[i0,200:400]), ta2[i0,200:400].min(), ta2[i0,200:400].max()))
            
            if type(ta)==type(np.ndarray(0)):
                ta[i0,data['CHANNEL_FLAG'][i0,:]>0] = 0.0
            else:
                ta[i0,ta[i0,:].mask>0] = 0.0
                ta2[i0,ta2[i0,:].mask>0] = 0.0

    
        # now we have to save the data in a FITS file
        # Note: changing the output data format to
        # primarily passing through the header and the data table
        # and only few variables and the header keys are updated or added
        
        
        tsyseff_avg[k] = np.nanmean(tsys[:,200:400])
        aTsyseff[osel,:] = tsyseff
        ahcorr[osel,:] = hcorr
        aspref[osel,:] = spref
        afrac[osel,0] = fraca
        afrac[osel,1] = fracb
        aTa[osel,:] = ta
        aspref2[osel,:] = spref2
        aTa2[osel,:] = ta2
    
        keys = data.dtype.names
        if 'spec' in keys:
            dkey = 'spec'
        elif 'DATA' in keys:
            dkey = 'DATA'
            
        data[dkey][osel,:] = ta.data
        data['CHANNEL_FLAG'] [osel,:] = ta.mask
        print('ta: ', ta.data[:5,:10])
        print('data: ', (data[dkey][osel,:])[:30,:10])
        
        
    tred = Time(datetime.datetime.now()).fits
    
    # updating header keywords
    hdr['DLEVEL'] = 0.95
    hdr['PROCTIME'] = tred
#    hdr.insert('VLSR', ('PROC_LEV', 0.9, 'pipeline processing level'), after=True)
#    hdr.add_comment('Level 0.9 Pipeline Processing', before='PROC_LEV')
    hdr.add_comment('Level 0.9 Pipeline Processing', after='VLSR')
    hdr.set('pgpixst', value=pixel_st, comment='pixel index where good pixels start')
    hdr.set('pgpixen', value=pixel_en, comment='pixel index of upper good pixel range')
    hdr.add_comment('L0.9 processing time: '%(tred.split('T')))
    hdr.add_comment('L0.9 version: %s'%(__version__))
    hdr.add_comment('L0.9 last pipeline update: %s'%(__updated__))
    hdr.add_comment('L0.9 developer: %s'%(__author__))
    
    os.makedirs(outDir, exist_ok=True)
    ofile = os.path.join(outDir, os.path.split(dfile)[1].replace('.fits','_L09.fits'))
    fits.writeto(ofile, data=None, header=hdr, overwrite=True)
    fits.append(ofile, data=data, header=hdr1)
    
    # if in debug mode, add more extensions
    if debug & datavalid:
        # data needed
        # Tsys1, Tsys2, fraca, fracb, Tsyseff, hcorr, spref
        # we have spectra/data for each otf spectrum: Tsyseff, hcorr, spref, fraca, fracb
        # single spectra: Tsys1, Tsys2, 
        # small groupf of spectra: hots
        col1 = Column(aTsyseff, name='Tsyseff', description='effective Tsys per OTF spectrum')
        col2 = Column(ahcorr, name='hcorr', description='effective hot spectrum')
        col3 = Column(aspref, name='spref', description='effective hot for ref spectrum')
        col4 = Column(aspref2, name='spref2', description='ref spectrum')
        col5 = Column(afrac, name='frac', description='fraction of Tsys1/2')
        col6 = Column(spec, name='spec', description='original spectra')
        col7 = Column(aTa2, name='tant', description='spectrum without hot correction')
        fT = Table([col1, col2, col3, col4, col5, col6, col7])
        fits.append(ofile, data=fT.as_array())

        col21 = Column(np.vstack(aghots), name='hots', description='averaged hot spectra')
        col22 = Column(np.hstack(aghtim), name='htime', description='utime of avg. hots')
        col23 = Column(np.hstack(aghmix), name='hmixer', description='mixer of avg. hots')
        col24 = Column(np.hstack(aghtint), name='htint', description='integration time of avg. hots')
        fgh = Table([col21, col22, col23, col24])
        fgh.write(ofile, append=True)

        col31 = Column(np.array(atsys), name='Tsys', description='Tsys before/after OTF')
        col32 = Column(np.array(atsmix), name='tsmix', description='mixer of Tsys')
        fits.append(ofile, data=Table([col31,col32]).as_array())

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

    
