#!/usr/bin/env python
"""
This is the GUSTO L09 Pipeline.
"""

__date__ = '9/19/2024'
__updated__ = '20241031'
__version__ = '0.2'
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
        cfi['gprocs']['debug'] = arg.debug
        

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
    
    
    for line in lines:
        if verbose:
            print(line)
        # identify the files for processing
        inDir = cfi['gdirs']['L08DataDir']
        outDir = cfi['gdirs']['L09DataDir']
        if line=='CII':
            filter = '*.fits'
        else:
            filter = '*.fits'
        print('outDir: ', outDir)
        print('filter: ', os.path.join(inDir,filter))
        
        # sdirs = sorted(glob.glob(os.path.join(inDir,filter), root_dir=inDir))
        #print(glob.glob(os.path.join(inDir,filter)))
        sdirs = sorted(glob.glob(os.path.join(inDir,filter)))
        print('single result: ', sdirs[0], os.path.split(sdirs[0]))
        dsc = [int(os.path.split(sdir)[1].split('_')[1].split('.')[0]) for sdir in sdirs]
        
        sdirs.sort(key=lambda sdirs: dsc)
        
        dfiles = []
        for i,ds in enumerate(dsc):
            if (ds >= scanRange[0]) & (ds <= scanRange[1]):
                dfiles.append(sdirs[i])
                        
        n_ds = len(dfiles)
        if int(cfi['gprocs']['max_files']) > 0:
            n_ds = int(cfi['gprocs']['max_files'])
            dfiles = dfiles[:n_ds]
        
        paramlist = [[a, b, c, d] for a in [line] for b in [inDir] for c in [outDir] for d in dfiles]
        # paramlist = [[a, b, c, d, e] for a in [line] for b in [inDir] for c in [outDir] for d in dfiles for e in worker_configurer]
        
        if verbose:
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
    line, inDir, outDir, dfile = params[0], params[1], params[2], params[3]
    
    # define some processing data first (maybe relocat to function later?)

    if 'ACS3' in dfile:
        pfr_ra = np.array(list([[-0.01,0.01],[0.25,0.38],[2.25,2.33]]), dtype=float)
        pp_ra = np.array([[20,83],[132,138],[200,205],[265,274],[323,340],[399,407],[498,509]], dtype=int)
        # all pixels above this value are masked as bad
        pixel_cut = 600
        band = 2
        add = 'B2'
        Tsky = 45  # Kelvin
        rfreq = 1900.5369  # GHz
    else:
        pfr_ra = np.array(list([[-0.01,0.01],[0.33,0.37],[2.26,2.31],[3.92,3.96]]), dtype=float)
        pp_ra = np.array([[23,50],[65,71],[100,103],[163,170],[198,206],[300, 511]], dtype=int)
        # all pixels above this value are masked as bad
        pixel_cut = 300
        band = 1
        add = 'B1'
        Tsky = 33.5  # Kelvin
        rfreq = 1461.131406


    #logger.info('loading: ', os.path.join(inDir,dfile), ' for line: ', line)
    spec, data, hdr = loadL08Data(os.path.join(inDir,dfile), verbose=False)
    rowFlag = data['ROW_FLAG']
    
    # for now, process all mixers
    umixers = np.unique(data['MIXER'])
    for mix in umixers:
        # first check crudely if we have enough data of various scan_types
        otfID, rfsID, rhsID, hotID = getSpecScanTypes(mix, spec, data, hdr)
        check = (np.argwhere(data['scan_type']=='REF').size > 3) & \
                (np.argwhere(data['scan_type']=='HOT').size > 3) & \
                (np.argwhere(data['scan_type']=='REFHOT').size > 3) & \
                (np.argwhere(data['scan_type']=='OTF').size > 5) & \
                (otfID.size>0) & (rfsID.size>0) & (rhsID.size>0) & (hotID.size>0)
        if not check:
            print('mix, dfile')
            print('specs: ', spec.shape)
            print('REFs: ', np.argwhere(data['scan_type']=='REF').size)
            print('HOTs: ', np.argwhere(data['scan_type']=='HOT').size)
            print('REFHOTs: ', np.argwhere(data['scan_type']=='REFHOT').size)
            print('OTFs: ', np.argwhere(data['scan_type']=='OTF').size)
            print('Not enough data available for processing')
            return 0
        
        tsys, refs, rhots, rtime, htime, Thot, Tsky = getCalSpectra(mix, spec, data, hdr, verbose=True)
        # tsys is a masked array if valid or an int if no good
        if type(tsys)==type(0):
            print('No Tsys available! Stop processing mix of dfile ', mix, dfile, tsys)
            # logger.error('No Tsys available! Stop processing mix of dfile ', mix, dfile, tsys)
            # logger.info('Tsys: ', tsys)
            break
        #print('<Tsys>: ', np.nanmean(tsys))
        tsys.fill_value = 0.0
        #print('tsys shape: ', tsys.shape)
        #print(list(tsys))
        prange = [40, 350]
        #pxs = np.arange(n_pix)
        
        
        otfID, rfsID, rhsID, hotID = getSpecScanTypes(mix, spec, data, hdr, verbose=verbose)
        
        # osel = np.argwhere((otfID == data['scanID']) & (otfID.size>=1) & (rfsID.size>2) & (rhsID.size>2) & (hotID.size>2) & (mix == data['MIXER']) & (data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0))
        osel = np.argwhere((otfID == data['scanID']) & (rfsID.size>=1) & (rhsID.size>=1) & (hotID.size>=1) & (otfID.size>=1) & (mix == data['MIXER']) & (data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0))
        print('otfID.size: ', otfID.size)
        if len(osel) > 0:
            # print('processing OTFs')
            # print('OTFs: ', otfID)
            # print('REFs: ', rfsID)
            # print('REFHOTs: ', rhsID)
            # print('HOTs: ', hotID)
            pass
        else:
            print('WARNING: No OTF spectra available.')
            # logger.warning('No OTF spectra available.')
            return 0
    
        spec_OTF = np.squeeze(spec[osel,:])
        stime = data['UNIXTIME'][osel]
        btime = (rtime[0] + htime[0]) / 2. # before OTFs
        atime = (rtime[1] + htime[1]) / 2. # after OTFs
        fracb = (stime - btime) / (atime - btime)
        fraca = (atime - stime) / (atime - btime)
        
        n_OTF, n_pix = spec_OTF.shape
        # antenna temperature is a masked array
        ta = ma.zeros([n_OTF, n_pix])
        ta.mask = spec.mask
        tsyseff = np.zeros([n_OTF, n_pix])
    
        # create the calibrated spectra
        for i0 in range(n_OTF):
            tsyseff[i0,:] = fracb[i0] * tsys[0,:] + fraca[i0] * tsys[1,:]
            # we might have to replace flagged pixels in tsyseff to not cause a problem in the spectra
            # => skipped for now since pixel flags should be very similar to flagged pixels in spectra
            spref = fracb[i0] * refs[0,:] + fraca[i0] * refs[1,:]
            ta[i0,:] = 2.*tsyseff[i0,:] * (spec_OTF[i0,:] - spref)/spref
            
            if type(ta)==type(np.ndarray(0)):
                ta[i0,data['CHANNEL_FLAG'][i0,:]>0] = 0.0
            else:
                ta[i0,ta[i0,:].mask>0] = 0.0

    
        # now we have to save the data in a FITS file
        
        # Table for data extension:
        ra = data['RA'][osel]
        dec = data['DEC'][osel]
        ra = data['RA']
        dec = data['DEC']
        rflag = data['ROW_FLAG'][osel] 
        
        c0 = SkyCoord(ra[0], dec[0], unit=(u.deg, u.deg), frame='icrs')
        c1 = SkyCoord(ra[osel]*u.degree, dec[osel]*u.degree, frame='icrs')
        rc1 = c1.transform_to(SkyOffsetFrame(origin=c0)) 
        raoff = rc1.lon.deg
        decoff = rc1.lat.deg 
        
        tsyseff_avg = np.nanmean(tsyseff[:,200:400], axis=1)
        timeobs = Time(data['UNIXTIME'][osel], format='unix').fits
        dobs = [tt[0].split('T')[0] for tt in timeobs]
        tobs = [tt[0].split('T')[1] for tt in timeobs]
        tred = Time(datetime.datetime.now()).fits
        
        aa = np.zeros(n_OTF)
        ao = np.ones(n_OTF)
        
        
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
        
    return dfile

    
    
    
    
