#!/usr/bin/env python3
"""
This is the GUSTO L09 Pipeline.
"""

__updated__ = '20251003'

import glob
import numpy as np
import numpy.ma as ma
import time
import sys
import os
import subprocess
import logging
from datetime import datetime
from astropy import units as u
from astropy.time import Time
from astropy.io import fits
from multiprocessing.pool import Pool
from PyAstronomy import pyasl

from .DataIO import loadL08Data
from .Utils import *
from .Logger import *
from .Configuration import *
from .flagdefs import *

commit_info = ''

def runGitLog():
    try:
        result = subprocess.run(['git', 'log', '-1', '--format=%cd', '--date=format-local:%Y-%m-%d %H:%M:%S %Z', '--pretty=format:Level 0.9 commit %h by %an %ad', '--', 'L09_calibrate.py'], capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        return f"Error: {e}"

    
def despike_polyRes(x, data, cflags, start, stop, points=60, count=3, deg=2, dx=1, stdlim=4.0):
    mask = np.zeros(len(data), dtype=bool)
    iin, iout = pyasl.slidingPolyResOutlier(x[start:stop], data[start:stop], points=points, count=count, deg=deg, stdlim=stdlim, controlPlot=False, dx=dx, mode='both')
    shifted = [i+start for i in iout]
    mask[shifted] = True
    cflags[mask] |= ChanFlags.SPUR 
    newdata = np.ma.masked_array(data, mask)
    return newdata, cflags


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
    """Function processing the Level 0.8 data. Input are uncalibrated 
    REF, HOT, and OTF spectra and output are calibrated OTF spectra
    """
    
    if args.debug ==True:
        print('\nExecuting debug mode.\n')
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)
        n_procs = 1
    elif args.cpus.isnumeric():
        n_procs = int(args.cpus)
    else:
        n_procs = multiprocessing.cpu_count()

    print('Number of cores used for processing: %i'%(n_procs))
    
    # get lines for processing
    band = args.band
    
    # these are for the Feb 25 data sets
    ignore = [10086, 13638, 17751, 27083, 28089, 4564, 7165, 7167]
    
    # load ranges for 2nd pixel masking
    
    inDir = args.path + 'level0.8/'
    outDir = args.path + 'level0.9/'
    os.makedirs(outDir, exist_ok=True)

    commit_info = runGitLog()
    
    for bandNum in band:
        if verbose:
            print('\nProcessing band:', int(bandNum))
        # identify the files for processing
        if int(bandNum) == 1:
            filter='NII*.fits'
        else:
            filter='CII*.fits'
        print('outDir: ', outDir)
        print('filter: ', os.path.join(inDir,filter))
        
        sdirs = sorted(glob.glob(os.path.join(inDir,filter)))
        dsc = [int(os.path.split(sdir)[1].split('_')[2].split('.')[0]) for sdir in sdirs]
        
        sdirs.sort(key=lambda sdirs: dsc)
        
        dfiles = []
        for i,ds in enumerate(dsc):
            if (ds >= scanRange[0]) & (ds <= scanRange[1]) & (ds not in ignore):
                dfiles.append(sdirs[i])
                        
        n_ds = len(dfiles)
        pxrange = [90, 700]
        rowflagfilter = 4294967295
        calmethod = args.calmethod
        despurmethod = args.despurmethod
        spurchannelfilter = args.spurchannelfilter
        
        params = {'band': int(bandNum), 'inDir': inDir, 'outDir': outDir,
                  'calmethod': calmethod, 'despurmethod': despurmethod,
                  'spurchannelfilter': spurchannelfilter, 'debug': args.debug, 'verbose': verbose,
                  'pxrange': pxrange, 'rowflagfilter': rowflagfilter}
        paramlist = [[a, b] for a in dfiles for b in [params]]

        if verbose:
            print('Number of data files: ', n_ds, len(sdirs))
            print('Selected data files: ', dfiles)
        
        # setup multiprocessing loop here to process each file in list
        with Pool(n_procs) as pool:
            for result in pool.imap(processL08, paramlist):
                print(f'Processed: {result}', flush=True)
        
    return n_ds



def cal_scaledGainHOTs(sspec, band, cflags, hgroup, closest, ghots, tsys, yfac):
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
    fit = np.polyfit(x_fit, y_fit, 1)
    baseline = np.poly1d(fit)
    Ta = Ta - baseline(xaxis)
    
    Ta, cflags = despike_polyRes(xaxis, Ta, cflags, band*40, band*75, points=20*band, count=1, deg=1, stdlim=3)
    Ta, cflags = despike_polyRes(xaxis, Ta, cflags, band*80, band*105, points=20*band, count=1, deg=1, stdlim=3)
    Ta, cflags = despike_polyRes(xaxis, Ta, cflags, band*150, band*180, points=20*band, count=1, deg=1, stdlim=3)
    Ta, cflags = despike_polyRes(xaxis, Ta, cflags, band*210, band*245, points=20*band, count=1, deg=1, stdlim=3)
    Tsys_median = 2.0*np.ma.median(tsyseff[band*40:band*240])
    rms = 0.33*(np.std(Ta[band*40:band*60]) + np.std(Ta[band*75:band*95]) + np.std(Ta[band*250:band*300]))
    return Ta, cflags, Tsys_median, rms


def processL08(paramlist):
    """Function processing the Level 0.8 data. Input are uncalibrated 
    REF, HOT, and OTF spectra and output are calibrated OTF spectra
    """
    TSKY = [33, 45]
    dfile = paramlist[0]
    params = paramlist[1]
    band, inDir, outDir, calmethod, debug = params['band'], params['inDir'], params['outDir'], params['calmethod'], params['debug']
    verbose = params['verbose']
    rowflagfilter = params['rowflagfilter']
    #pxrange = (int(params['pxrange'][0]), int(params['pxrange'][1]))      # good pixel ranges

    spec, data, hdr, hdr1 = loadL08Data(os.path.join(inDir,dfile), verbose=False)
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
            print('GUSTO Pipeline: Not enough data available for processing. ROW_FLAGs are set appropriately. ')
            data['ROW_FLAG'][msel] |= RowFlags.BAD_PHASE   # flagged as missing data
            datavalid[k] = False
            return 0
        
        tsys, refs, rhots, rtime, htime, Thot, rhIDs, yfac = getCalSpectra(mix, spec, data, hdr, rowflagfilter=rowflagfilter, Tsky=Tsky, verbose=True)
        # tsys is a masked array if valid or an int if no good
        if type(tsys)==type(0):
            print('No Tsys available! Stopped processing of mixer %i in dfile %s'%(mix, dfile), tsys)
            datavalid[k] = False
            continue
                
        otfID, rfsID, rhsID, hotID = getSpecScanTypes(mix, spec, data, hdr, rowflagfilter=rowflagfilter, verbose=verbose)
        osel = np.argwhere((otfID == data['scanID']) & (rfsID.size>=1) & (rhsID.size>=1) & 
                           (otfID.size>=1) & (mix == data['MIXER']) & (data['scan_type'] == 'OTF') & 
                           (rflag)).flatten()
        if len(osel) > 2:
            pass
        else:
            print('WARNING: No OTF spectra available for mixer: ', mix)
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
            print('Encountered problem with HOT groups. Flaging rows.')
            data['ROW_FLAG'][msel] |= RowFlags.BAD_DATA   # flagged as missing data
            continue

        # reduce the assignment to the OTF spectra only
        hgroup = ahgroup[osel]
        # create the calibrated spectra
        for i0 in range(n_OTF):
            # fixme: make this conditional.  if calmethod == 'cal_scaledGainHOTs'
            ta[i0,:], cflags_OTF[i0], Tsys_OTF[i0], rms_OTF[i0] = cal_scaledGainHOTs(spec_OTF[i0,:], band, cflags_OTF[i0], hgroup, hgroup[i0], ghots, tsys, yfac)
            
        # now we have to save the data in a FITS file
        data['DATA'][osel,:] = ta.data        
        data['CHANNEL_FLAG'] [osel,:] = cflags_OTF
        data['Tsys'][osel] = Tsys_OTF
        data['rms'][osel] = rms_OTF
        
    if np.any(datavalid) == False:
        print('Not enough data available for saving. ')
        return 0

    tred = Time(datetime.datetime.now()).fits
    
    # updating header keywords
    hdr.set('DLEVEL', value = 0.9, after='PROCTIME')
    hdr['PROCTIME'] = tred
    
    hdr.set('L09_TIME', value=tred, comment=('L0.9 pipeline processing time'))
    hdr.set('rwflfilt', value=rowflagfilter, comment='applied rowflag filter for useful spectra')
    hdr.set('rhID1', value=rhIDs[0], comment='scan ID for first REFHOT/REF')
    if len(rhIDs) > 1:
        hdr.set('rhID2', value=rhIDs[1], comment='scan ID for second REFHOT/REF')
    hdr.set('spurfltr', value=params['spurchannelfilter'], comment='was a spur channel filter pre-applied')
    hdr.set('despur', value=params['despurmethod'], comment='despur processing method applied')
    hdr.set('calmethd', value=calmethod, comment='calibration processing method applied')
    hdr.set('', value='')
    hdr.set('', value='', after='VLSR')
    hdr.set('', value='          Level 0.9 Pipeline Processing', after='VLSR')
    hdr.set('', value='', after='VLSR')
    hdr.add_comment('L0.9 processing time: %s'%(tred))
    hdr.add_comment('L0.9 code update date: %s'%(__updated__))
    hdr.add_comment(commit_info)
    hdr.set('', value='', after='calmethd')
    
    os.makedirs(outDir, exist_ok=True)
    ofile = os.path.join(outDir, os.path.split(dfile)[1].replace('.fits','_L09.fits'))
    fits.writeto(ofile, data=None, header=hdr, overwrite=True)
    fits.append(ofile, data=data, header=hdr1)
    
    return dfile

    
