"""
This is the GUSTO Level 0.9 Pipeline.
"""
import numpy as np
import numpy.ma as ma
import os
import subprocess
from tqdm import tqdm
from datetime import datetime
from astropy import units as u
from astropy.time import Time
from astropy.io import fits
from multiprocessing.pool import Pool
from PyAstronomy import pyasl

from .DataIO import *
from .Logger import *
from .Configuration import *
from .flagdefs import *

logger = logging.getLogger('pipelineLogger')

def despike_polyRes(x, data, cflags, start, stop, points=60, count=3, deg=2, dx=1, stdlim=4.0):
    mask = np.zeros(len(data), dtype=bool)
    iin, iout = pyasl.slidingPolyResOutlier(x[start:stop], data[start:stop], points=points, count=count, deg=deg, stdlim=stdlim, controlPlot=False, dx=dx, mode='both')
    shifted = [i+start for i in iout]
    mask[shifted] = True
    cflags[mask] |= ChanFlags.SPUR 
    newdata = np.ma.masked_array(data, mask)
    return newdata, cflags


def getSpecScanTypes(mixer, spec, data, hdr, rowflagfilter=0, verbose=False):
    """Function calculating the calibration spectrum for a single mixer.

    Parameters
    ----------
    mixer : int
        mixer number
    spec : masked float array
        array containing the latest masked spectra
    data : FITS rec array
        record array containing the GUSTO data fron the FITS file

    """
    mixers  = data['MIXER']
    scanID = data['scanID']
    scan_type = data['scan_type']
    rflags = checkRowflag(data['ROW_FLAG'], rowflagfilter=rowflagfilter)

    rfsel = (mixers == mixer) & (scan_type == 'REF') & (rflags)
    rfsID = np.unique(scanID[rfsel])
    rhsel = (mixers == mixer) & (scan_type == 'REFHOT') & (rflags)
    rhsID = np.unique(scanID[rhsel])

    # identify the HOT spectra
    htsel = (mixers == mixer) & (scan_type == 'HOT') & (rflags)
    hotID = np.unique(scanID[htsel])

    # identify the OTF spectra
    otsel = (mixers == mixer) & (scan_type == 'OTF') & (rflags)
    otfID = np.unique(scanID[otsel])
    return otfID, rfsID, rhsID, hotID


def getWeights(spec, hdr, selectS, selectR, cutoff=0.3):
    """Function getWeights imports a scanID and selections for S and R and computes
    the rms of (S-R)/R and returns weights accordingly. 

    Parameters
    ----------
    spec : spectra from a scanID
    hdr  : header from that scanID
    selectS : selection of indices corresponding to the S array
    selectR : index corresponding to the R array
    cutoff : optional cutoff in weights below which the weight is set to 0
    
    Returns
    -------
    an array of weights scaled by 1/rms, with a low cutoff below which the weight is set to 0

    Usage
    --------
    Typically used for the coadding of frames with multiple subscans like HOT, REF and REFHOT
    
    """
    rms = []
    chan = [512, 1024]
    fScale = [5000/511.0, 5000/1023.0]
    band = int(hdr['BAND'])
    xaxis = np.arange(0, chan[band-1]*fScale[band-1], fScale[band-1])

    for s in selectS:
        specS = spec[s,:]
        specR = spec[selectR,:]
        calspec = ((specS-specR)/specR).flatten()
        x_fit = xaxis[band*40:band*60]
        x_fit = np.append(x_fit, xaxis[band*250:band*300], axis=0)
        y_fit = calspec[band*40:band*60]
        y_fit = np.append(y_fit, calspec[band*250:band*300], axis=0)
        fit = np.polyfit(x_fit, y_fit, 1)
        baseline = np.poly1d(fit)
        calspec = calspec - baseline(xaxis)
        rms.append(np.std(calspec[band*75:band*95]) + np.std(calspec[band*250:band*300]))

    weights = np.mean(rms)/rms
    weights = np.where(weights < cutoff, 0.0, weights)
    return weights


def getCalSpectra(mixer, spec, data, hdr, rowflagfilter=0, Tsky=45., verbose=False):
    """Function calculating the calibration spectrum for a single mixer.

    Parameters
    ----------
    mixer : int
        mixer ID
    spec : masked float array
        array containing the latest masked spectra
    data : FITS rec array
        record array containing the GUSTO data fron the FITS file
    rowflagfilter : int
        filter value for rowflags
    Tsky : float
        Sky temperature at wavelength of observation.

    Returns
    -------
    returns averaged noise temperature Tsyss spectra, REFs spectra, RHOTs spectra, 
    and averaged REF and HOT times, rtime and htime, respectively for the spectra 
    before and after the OTF scans. The Tsys time would be the average of rtime 
    and htime for before and after the OTFs.

    """
    mixers  = data['MIXER']
    scanID = data['scanID']
    
    Thot   = float(hdr['THOT'])
    seqflag = int(hdr['SEQ_FLAG'])
    scan_type = data['scan_type']
    rflags = checkRowflag(data['ROW_FLAG'], rowflagfilter=rowflagfilter)
    ch_flag = data['CHANNEL_FLAG']   # spectral pixel (or channel) mask
    stime = data['UNIXTIME']

    otfID, rfsID, rhsID, hotID = getSpecScanTypes(mixer, spec, data, hdr, rowflagfilter=rowflagfilter, verbose=verbose)

    try:
        bef = rhsID[rhsID<otfID]
        aft = rhsID[rhsID>otfID]
    except:
        bef = rhsID[rhsID<otfID[0]]
        aft = rhsID[rhsID>otfID[0]]
    if (len(bef)>0) & (len(aft)>0):
        rhIDbef = bef[np.argmax(bef)]
        rhIDaft = aft[np.argmin(aft)]
        rhIDs = [rhIDbef, rhIDaft]
    elif (len(aft)>0):
        rhIDaft = aft[np.argmin(aft)]
        rhIDs = [rhIDaft]
        logger.debug('getCalSpectra: Only REFHOT after OTF available')
    elif (len(bef)>0):
        rhIDbef = bef[np.argmax(bef)]
        rhIDs = [rhIDbef]
        logger.debug('getCalSpectra: Only REFHOT before OTF available')
    else:
        logger.debug('getCalSpectra: Not enough REFHOT scans available (before/after OTF)')
        return -999, 0, 0, 0, 0, 0, [0,0], 0

    Tsyss = []
    REFs = []
    RHOTs = []
    rtimes = []
    htimes = []
    yfac = []
    for rhID in rhIDs:
        rsel = np.argwhere((rhID == scanID) & (mixer == mixers) & (scan_type == 'REF') & (rflags))
        hsel = np.argwhere((rhID == scanID) & (mixer == mixers) & (scan_type == 'REFHOT') & (rflags))
        htime = stime[hsel].mean()
        rtime = stime[rsel].mean()
        osel = np.argwhere((mixer == mixers) & (scan_type == 'OTF') & (rflags))
        closest = np.argmin(np.abs(stime[osel] - rtime))
        ohsel = np.argwhere((mixer == mixers) & (scan_type == 'HOT') & (rflags))
        hot_closest = np.argmin(np.abs(stime[ohsel] - htime))
        
        weights = getWeights(spec, hdr, rsel, osel[closest]) # compare REFs w/ closest OTF 
        spec_r = np.average(spec[rsel,:], axis=0, weights=weights)

        weights = getWeights(spec, hdr, hsel, ohsel[hot_closest]) # compare REFHOTs w/ closest OTF HOT
        spec_h = np.average(spec[hsel,:], axis=0, weights=weights)

        # determine yfactor
        # estimate Tsys for each Device, correct for backend gain slope of 1.3x
        y_factor  = spec_h/spec_r
        y_fixed = (y_factor-1.)/1.3 + 1.

        tsys = np.squeeze((Thot - Tsky*y_fixed[:])/(y_fixed[:] - 1.))

        Tsyss.append(tsys)
        yfac.append(y_factor)
        REFs.append(spec_r)
        RHOTs.append(spec_h)
        rtimes.append(rtime)
        htimes.append(htime)

    Tsyss = np.ma.array(Tsyss).squeeze()
    yfac = np.ma.array(yfac).squeeze()
    REFs = np.ma.array(REFs).squeeze()
    RHOTs = np.ma.array(RHOTs).squeeze()
    return Tsyss, REFs, RHOTs, rtimes, htimes, Thot, rhIDs, yfac


def getHotInfo(spec, data, hdr, mixer, dfile='', verbose=False, rowflagfilter=0):
    """Function analyzing and processing HOT spectra in an 
    OTF strip.
    There are caveats in the current data (tbc) including that there are duplicate
    REFs at the end of the sequence. This issue has been mitigated for now by
    ignoring these REFs. The current assumption, yet unverified, is that there is 
    always a REF at the beginning of the sequence. <- asumption is wrong! REF/REFHOT can be at end of sequence

    Parameters
    ----------
    spec : array
            array containing only the spectra for the OTF scan
    data : recarray
            recarray with the binary table information from the FITS file
    mixer : float
            mixer processed

    Returns
    -------
    The function returns 4 arrays: hgroup, ghots, ghtim 
    hgroup: the HOT group information for all spectra in data set
    ghots: the averaged HOTs for each HOT group
    ghtim: the average unixtime for each HOT group
    ghtint: the integration time for the averaged hots
    """
    n_spec, n_pix = spec.shape
    umixers = np.unique(data['MIXER'])
    mx = np.argwhere(umixers==mixer).flatten()
    
    unixtime = data['UNIXTIME']
    tint = data['INTTIME']
    ut0 = unixtime[0]
    
    hgroup = np.zeros(n_spec, dtype=int)
    hcnt = 0   # counter counting the total number of hots
    hgrp = 0   # counter for hot groups

    # determine the first hot or refhot scan based on time
    lasthot = np.min(unixtime[np.argwhere((data['MIXER']==umixers[mx])&((data['scan_type']=='HOT')|(data['scan_type']=='HOT')))])
    lasthot = np.min(np.argwhere(unixtime == lasthot))
    firstref = True if (np.argwhere((data['MIXER']==umixers[mx])&(data['scan_type']=='REF')).min() < 10) else False

    # added check for REFHOT duplicates at beginning and at end of sequence
    # there should only be one REFHOT at the beginning and one at
    # the end
    rflags = checkRowflag(data['ROW_FLAG'], rowflagfilter=rowflagfilter)
    rhscans = data['scanID'][(data['scan_type']=='REFHOT')&(rflags)&(data['MIXER']==mixer)]
    rfscans = data['scanID'][(data['scan_type']=='REF')&(rflags)&(data['MIXER']==mixer)]
    otscans = data['scanID'][(data['scan_type']=='OTF')&(rflags)&(data['MIXER']==mixer)]
    htscans = data['scanID'][(data['scan_type']=='HOT')&(rflags)&(data['MIXER']==mixer)]
    urhs = np.unique(rhscans)
    urfs = np.unique(rfscans)
    uots = np.unique(otscans)
    uhts = np.unique(htscans)

    # the first OTF scan ID should be the one determining the sequence
    if uots.size == 0 or uots.size > 3:
        # bad sequence with none or too many OTFs
        logger.debug("SEQUENCE HAS THE WRONG NUMBER OF OTFs")
        return 0,0,0,0
    
    uflag = np.zeros(n_spec)    
    for i in range(n_spec):
        if data['MIXER'][i] == mixer:
            
            if (data['scan_type'][i] == 'HOT')|((data['scan_type'][i] == 'REFHOT')&(data['scanID'][i] in urhs)):
                if hcnt==0:
                    lasthot = i
                else:
                    if unixtime[i] - unixtime[lasthot] < 4.0:
                        pass
                    else:
                        hgrp += 1
                    lasthot = i
                uflag[i] = 1
                hcnt += 1
            else:
                uflag[i] = 0
            hgroup[i] = hgrp
            if ((data['scan_type'][i] == 'REF')&(data['scanID'][i] in urhs)):
                uflag[i] = 1
            
    maxgrp = 0
    ghots = np.zeros((int(hgroup.max()+1), n_pix))
    ghtim = np.zeros(int(hgroup.max()+1))
    ghtint = np.zeros(int(hgroup.max()+1))
    
    # now we have to average the hots for each group flagged for use:
    # for i, mx in enumerate(umixers):
    for i, mx in enumerate([mixer]):
        maxgrp = int(hgroup[data['MIXER']==mx].max()+1)
        for j in range(maxgrp):
            sel = np.argwhere((data['MIXER']==mx) & (hgroup==j) & (rflags) & (uflag==1)).flatten()
            hsel = np.argwhere((data['MIXER']==mx) & (hgroup!=j) & (rflags) & (uflag==1)).flatten()
            htime = unixtime[sel].mean()
            closest = np.argmin(np.abs(unixtime[hsel] - htime))        

            weights = getWeights(spec, hdr, sel, hsel[closest])
            ghots[j,:] = np.average(spec[sel,:], axis=0, weights=weights)
            ghtim[j] = np.mean(unixtime[sel])
            ghtint[j] = np.sum(tint[sel])
    
    # loop through OTF spectra, identify closest HOT group in time, return index as hgroup
    for i in range(n_spec):
        if data['MIXER'][i] == mixer:
            hgroup[i] = np.argmin(np.abs(ghtim - unixtime[i]))

    return hgroup, ghots, ghtim, ghtint


def checkRowflag(rowflagvalue, rowflagfilter=0):
    r"""Function checking bitflags from spectra to a bitmask of allowed bits.
    See also the vectorized version checkRowflags() below.

    Parameters
    ----------
    rowflagvalue : numpy int array
                rowflag (bitmask) of spectra
    rowflagfilter : int
                allowed rowflags (bitmask)

    Returns
    -------
    boolean array of True/False values

    Examples
    --------
    import enum
    import numpy as np
    
    rowflagfilter = RowFlags.NO_HK | RowFlags.MISSING_INT | RowFlags.MIXER_UNPUMPED
    print('rowflagfilter input: ', rowflagfilter, 'contained flags: ', enum.show_flag_values(rowflagfilter))
    
    # enum.show_flag_values(rowflagvalue) works only on single values
    flag_values = [enum.show_flag_values(rfv) for rfv in rowflagvalue]
    print('rowflagvalue input: ', rowflagvalue, 'contained flags: ', flag_values)
    """
    if type(rowflagvalue) == type(0):
        rowflagvalue = np.array([rowflagvalue])
    negfilt = ~rowflagfilter#.__invert__()
    return np.array([False if int(negfilt).__and__(int(rfvalue)) else True for rfvalue in rowflagvalue], dtype=bool)

checkRowflags = np.vectorize(checkRowflag)


    #########################################################
    # setup the processing loop for each file
    # the loop should use multiprocessing (parameters like # of cores, should be in setup file)
    # - each file is loaded
    # - determine the data type: OTF or PS to be handled separately
    # - for each file process the refs/hots -> y-factor -> Tsys
    # - for each file process the OTFs
    # - re-evaluate the row flags for each OTF spectrum
    # save spectra to FITS file
        
def L09_Pipeline(args, scanRange, verbose=False):
    global logger
    prefix = ['NII_', 'CII_'] 
    """Function processing the Level 0.7 data. Input are uncalibrated 
    REF, HOT, and OTF spectra and output are calibrated OTF spectra
    """
    
    if args.debug ==True:
        logger.debug('Executing debug mode.')
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)
        n_procs = 1
    elif args.cpus.isnumeric():
        n_procs = int(args.cpus)
    else:
        n_procs = multiprocessing.cpu_count()

    logger.info('Number of cores used for processing: %i'%(n_procs))
    
    sum_files = 0    
    inDir = args.path + 'level0.7/'
    outDir = args.path + 'level0.9/'
    os.makedirs(outDir, exist_ok=True)
    if args.erase:
        clear_folder(outDir)

    commit_info = runGitLog('0.9', 'L09_calibrate.py')
    
    for band in args.band:
        logger.info(f'Processing band %s'%(band))
        dfiles = makeFileGlob(inDir, prefix[int(band)-1], 'fits', scanRange)
        sum_files += len(dfiles)        
        pxrange = [90, 700]
        rowflagfilter = 4294967295
        calmethod = args.calmethod
        despurmethod = args.despurmethod
        spurchannelfilter = args.spurchannelfilter
        polyorder = args.polyorder
        
        params = {'band': int(band), 'inDir': inDir, 'outDir': outDir, 'polyorder': polyorder, 
                  'calmethod': calmethod, 'despurmethod': despurmethod,
                  'spurchannelfilter': spurchannelfilter, 'debug': args.debug, 'verbose': verbose,
                  'pxrange': pxrange, 'rowflagfilter': rowflagfilter, 'commit_info': commit_info}
        paramlist = [[a, b] for a in dfiles for b in [params]]

        if verbose:
            logger.debug('Number of data files: %i'%(len(dfiles)))
            logger.debug('Selected data files: %s'%(dfiles))
        
        # setup multiprocessing loop here to process each file in list
        with Pool(n_procs) as pool:
            list(tqdm(pool.imap_unordered(processL07, paramlist), total=len(dfiles), colour='yellow', leave=False))
            
        print()
    return sum_files



def cal_scaledGainHOTs(sspec, band, cflags, hgroup, closest, ghots, tsys, yfac, polyorder):
    chan = [512, 1024]
    fScale = [5000/511.0, 5000/1023.0]
    oldmed=999999
    best = [0.5, 0.5]
    Ta = ma.zeros(sspec.shape)
    tsyseff = ma.zeros(sspec.shape)
                
    xaxis = np.arange(0, chan[band-1]*fScale[band-1], fScale[band-1])
    seq_hots = ghots[closest,:]
    
    if closest+1 < hgroup.max():
        seq_hots = np.vstack([seq_hots, ghots[closest+1,:]])
    else:
        seq_hots = np.vstack([seq_hots, ghots[closest-1,:]])
    n_shots = seq_hots.shape[0]
                
    if yfac.ndim > 1:
        sRntest = [x/y for x,y in zip(seq_hots, yfac.squeeze())]
        sRntest = np.array(sRntest)
        yvalid = np.nonzero((yfac[0,:].squeeze() > 1.0))[0]
    else:
        sRntest = seq_hots / yfac.squeeze()
        yvalid = np.nonzero((yfac.squeeze() > 1.0))
    test = sspec - sRntest
    med = np.ma.median(test[:,band*40:band*295], axis=1)
    scale = 1 + med/np.ma.median(sRntest[:, band*40:band*295], axis=1)
    
    for a in np.arange(0.0, 1.01, 0.25):  # now loop over a and b with c constant
        for b in np.arange(0.0, 1.01, 0.25):
            if yfac.ndim == 1:
                yfac_eff = yfac
                tsyseff = tsys
                sRn = seq_hots*scale[:,None] / yfac_eff
                Ta = 2*tsyseff * (sspec - (a*sRn[0,:] + (1.0-a)*sRn[1,:]))/(a*sRn[0,:] + (1.0-a)*sRn[1,:])
            else:
                yfac_eff = b*yfac[0,:] + (1-b)*yfac[1,:]
                tsyseff = b*tsys[0,:] + (1-b)*tsys[1,:]
                sRn1 = seq_hots*scale[0] / yfac_eff
                sRn2 = seq_hots*scale[1] / yfac_eff
                Ta = 2*tsyseff * (sspec - (a*sRn1[0,:] + (1.0-a)*sRn2[1,:]))/(a*sRn1[0,:] + (1.0-a)*sRn2[1,:])
            x_fit = xaxis[band*40:band*60]
            x_fit = np.append(x_fit, xaxis[band*250:band*300], axis=0)
            y_fit = Ta[band*40:band*60]
            y_fit = np.append(y_fit, Ta[band*250:band*300], axis=0)
            fit = np.polyfit(x_fit, y_fit, 1)
            baseline = np.poly1d(fit)
            Ta = Ta - baseline(xaxis)
            med = np.std(Ta[band*75:band*95]) + np.std(Ta[band*250:band*300])
            if med < oldmed:
                oldmed = med
                best[0] = a
                best[1] = b
            if yfac.ndim == 1:
                break
    #print(best, oldmed)
    if yfac.ndim == 1:
        yfac_eff = yfac
        tsyseff = tsys
    else:
        yfac_eff = best[1]*yfac[0,:] + (1-best[1])*yfac[1,:]
        tsyseff = best[1]*tsys[0,:] + (1-best[1])*tsys[1,:]
    sRn = seq_hots*scale[:,None] / yfac_eff
    Ta = 2*tsyseff * (sspec - (best[0]*sRn[0,:] + (1.0-best[0])*sRn[1,:]))/(best[0]*sRn[0,:] + (1.0-best[0])*sRn[1,:])
    x_fit = xaxis[band*40:band*60]
    x_fit = np.append(x_fit, xaxis[band*75:band*95], axis=0)
    x_fit = np.append(x_fit, xaxis[band*250:band*300], axis=0)
    y_fit = Ta[band*40:band*60]
    y_fit = np.append(y_fit, Ta[band*75:band*95], axis=0)
    y_fit = np.append(y_fit, Ta[band*250:band*300], axis=0)
    fit = np.polyfit(x_fit, y_fit, polyorder)
    baseline = np.poly1d(fit)
    Ta = Ta - baseline(xaxis)
    
    Ta, cflags = despike_polyRes(xaxis, Ta, cflags, band*40, band*75, points=20*band, count=1, deg=1, stdlim=3)
    Ta, cflags = despike_polyRes(xaxis, Ta, cflags, band*80, band*105, points=20*band, count=1, deg=1, stdlim=3)
    Ta, cflags = despike_polyRes(xaxis, Ta, cflags, band*150, band*180, points=20*band, count=1, deg=1, stdlim=3)
    Ta, cflags = despike_polyRes(xaxis, Ta, cflags, band*210, band*245, points=20*band, count=1, deg=1, stdlim=3)
    Tsys_median = 2.0*np.ma.median(tsyseff[band*40:band*240])
    rms = 0.33*(np.std(Ta[band*40:band*60]) + np.std(Ta[band*75:band*95]) + np.std(Ta[band*250:band*300]))
    return Ta, cflags, Tsys_median, rms


def processL07(paramlist):
    """Function processing the Level 0.7 data. Input are uncalibrated 
    REF, HOT, and OTF spectra and output are calibrated OTF spectra
    """
    TSKY = [33, 45]
    dfile = paramlist[0]
    params = paramlist[1]
    band, inDir, outDir, polyorder, calmethod, debug, verbose, rowflagfilter, commit_info = \
        params['band'], params['inDir'], params['outDir'], params['polyorder'], params['calmethod'], \
        params['debug'], params['verbose'], params['rowflagfilter'], params['commit_info']
    #pxrange = (int(params['pxrange'][0]), int(params['pxrange'][1]))      # good pixel ranges

    logger.debug(f'Loading file: {dfile}')
    spec, data, hdr, hdr1 = loadSDFITS(os.path.join(inDir,dfile), verbose=False)
    band = hdr['band']
    Tsky = TSKY[band-1]
    if params['spurchannelfilter']:
        indices = np.where(data['CHANNEL_FLAG'] & ChanFlags.VARIABLE_SPUR)
        spec.mask[indices] = True
    else:
        spec.mask = np.zeros(spec.shape, dtype=bool)    
            
    n_spec, n_pix = spec.shape
    umixers = np.unique(data['MIXER'])     # for now, process all mixers
    n_umix = umixers.size
    amixer = np.zeros(n_umix, dtype=int)
    datavalid = np.ones(n_umix, dtype=bool)

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

        if not check:
            logger.debug(f'{mix} {dfile}')
            logger.debug(f'check: {check}')
            logger.debug(f'specs: {spec.shape}')
            logger.debug(f"REFs: {np.argwhere(data['scan_type']=='REF').size} {(np.argwhere(data['scan_type']=='REF').size > 3)}")
            logger.debug(f"HOTs: {np.argwhere(data['scan_type']=='HOT').size} {(np.argwhere(data['scan_type']=='HOT').size > 3)}")
            logger.debug(f"REFHOTs: {np.argwhere(data['scan_type']=='REFHOT').size} {(np.argwhere(data['scan_type']=='REFHOT').size > 3)}")
            logger.debug(f"OTFs: {np.argwhere(data['scan_type']=='OTF').size} {(np.argwhere(data['scan_type']=='OTF').size > 5)}")
            logger.debug(f'other: {(otfID.size>0)} {(rfsID.size>0)} {(rhsID.size>0)} {(hotID.size>0)}')
            logger.debug(f'IDs: {otfID} {rfsID} {rhsID} {hotID}')
            logger.debug(f"rowflagfilter: {rowflagfilter} {np.any(checkRowflag(data['ROW_FLAG'][msel], rowflagfilter=rowflagfilter))}")
            logger.debug('GUSTO Pipeline: Not enough data available for processing. ROW_FLAGs are set appropriately. ')
            data['ROW_FLAG'][msel] |= RowFlags.BAD_PHASE   # flagged as missing data
            datavalid[k] = False
            return 0
        
        tsys, refs, rhots, rtime, htime, Thot, rhIDs, yfac = getCalSpectra(mix, spec, data, hdr, rowflagfilter=rowflagfilter, Tsky=Tsky, verbose=True)
        # tsys is a masked array if valid or an int if no good
        if type(tsys)==type(0):
            logger.debug('No Tsys available! Stopped processing of mixer %i in dfile %s'%(mix, dfile))
            datavalid[k] = False
            continue
                
        osel = np.argwhere((np.isin(data['scanID'], otfID)) & (rfsID.size>=1) & (rhsID.size>=1) & 
                           (otfID.size>=1) & (mix == data['MIXER']) & (data['scan_type'] == 'OTF') & 
                           (rflag)).flatten()
        if len(osel) > 2:
            pass
        else:
            logger.debug('WARNING: No OTF spectra available for mixer: %i'%(mix))
            datavalid[k] = False
            data['ROW_FLAG'][msel] |= RowFlags.BAD_DATA   # flagged as missing data
            continue

        spec_OTF = np.squeeze(spec[osel,:])
        cflags_OTF =  np.squeeze(data['CHANNEL_FLAG'][osel,:])
        stime = data['UNIXTIME'][osel]
        Tsys_OTF = data['Tsys'][osel]
        rms_OTF = data['rms'][osel]
        n_OTF, n_opix = spec_OTF.shape
        # antenna temperature is a masked array
        ta = ma.zeros([n_OTF, n_opix])
        ahgroup, ghots, ghtim, ghtint = getHotInfo(spec, data, hdr, mix, dfile=dfile, rowflagfilter=rowflagfilter, verbose=True)
        if type(ahgroup)==type(0):
            logger.debug('Encountered problem with HOT groups. Flaging rows.')
            data['ROW_FLAG'][msel] |= RowFlags.BAD_DATA   # flagged as missing data
            continue

        # reduce the assignment to the OTF spectra only
        hgroup = ahgroup[osel]
        # create the calibrated spectra
        for i0 in range(n_OTF):
            # fixme: make this conditional.  if calmethod == 'cal_scaledGainHOTs'
            ta[i0,:], cflags_OTF[i0], Tsys_OTF[i0], rms_OTF[i0] = cal_scaledGainHOTs(spec_OTF[i0,:], band, cflags_OTF[i0], hgroup, hgroup[i0], ghots, tsys, yfac, int(polyorder))
            
        # now we have to save the data in a FITS file
        data['DATA'][osel,:] = ta.data        
        data['CHANNEL_FLAG'] [osel,:] = cflags_OTF
        data['Tsys'][osel] = Tsys_OTF
        data['rms'][osel] = rms_OTF
        
    if np.any(datavalid) == False:
        logger.debug('Not enough data available for saving. ')
        return 0

    tred = Time(datetime.now()).fits
    
    # updating header keywords
    hdr.set('DLEVEL', value = 0.9)
    hdr.set('rwflfilt', value=rowflagfilter, comment='applied rowflag filter for useful spectra')
    hdr.set('rhID1', value=rhIDs[0], comment='scan ID for first REFHOT/REF')
    if len(rhIDs) > 1:
        hdr.set('rhID2', value=rhIDs[1], comment='scan ID for second REFHOT/REF')
    hdr.set('polyordr', value=int(params['polyorder']), comment='order of baseline polynomial fit')
    hdr.set('spurfltr', value=params['spurchannelfilter'], comment='was a spur channel filter pre-applied')
    hdr.set('despur', value=params['despurmethod'], comment='despur processing method applied')
    hdr.set('calmethd', value=calmethod, comment='calibration processing method applied')
    hdr.set('', value='')
    hdr.set('', value='', after='VLSR')
    hdr.set('', value='          Level 0.9 Pipeline Processing', after='VLSR')
    hdr.set('', value='', after='VLSR')
    hdr.add_comment(commit_info)
    hdr.set('', value='', after='calmethd')
    hdr.add_history('Level 0.9 processed at %s'%(tred))

    os.makedirs(outDir, exist_ok=True)
    ofile = os.path.join(outDir, os.path.split(dfile)[1].replace('_L07.fits','_L09.fits'))
    fits.writeto(ofile, data=None, header=hdr, overwrite=True)
    fits.append(ofile, data=data, header=hdr1)
    logger.debug(f'Processed and saved file: {ofile}')
    
    return dfile

    
