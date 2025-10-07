"""
GUSTO Pipeline utility functions
"""

import numpy as np
from astropy.io import ascii, fits
from astropy.utils.exceptions import AstropyWarning
import sys
import os
import inspect
import datetime
import warnings
from enum import Enum, Flag, auto
from astropy.io import fits
from .flagdefs import *


__version__ = 0.12
__date__ = '20250927'
__updated__ = '20250928'
fheader, tail = os.path.split(inspect.stack()[0][1])
__pyfile__ = tail
__path__ = fheader


warnings.filterwarnings('ignore', category=Warning,
                        message=' FITSFixedWarning: ', append=True)

def lprint(*args, **kwargs):    
    print(*args, **kwargs)


def getRange(icpar, dtype='float', endpoint=True):
    """Function reading an input str and convert ranges to 
    an array of floats. 
    Uses getValues() and returns a numpy float array if <3 parameters are
    provided. 
        
    Parameters
    ----------
    icpar: str
        parameter str like '[start, stop, step]'
    """
    pars = getValues(icpar)    
    if pars.size < 3:
        return pars
    else:
        # create a grid from the parameters
        if endpoint:
            return np.arange(pars[0], pars[1]+pars[2], pars[2])
        else:
            return np.arange(pars[0], pars[1], pars[2])


def getValues(icpar, dtype='float'):
    """Function converting str data arrays to numpy float array.
        
    Parameters
    ----------
    icpar: str
        parameter str like '[start, stop, step]' or '3.45'
    """
    icpars = icpar.replace('[','').replace(']','').replace(' ','').split(',')
    ipars = np.array((icpars), dtype=dtype)
    return ipars


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
    if len(otfID) > 1:
        otfID = otfID[1:]
        print('reducing number of OTF IDs to 1')
    return otfID, rfsID, rhsID, hotID


def getWeights(spec, hdr, selectS, selectR, cutoff=0.33):
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
    if (len(otfID)>1):
        print('getCalSpectra mixer %i: Too many OTF scan IDs for processing: '%mixer, otfID)
        return -999, 0, 0, 0, 0, 0, 0, [0,0], 0

    bef = rhsID[rhsID<otfID]
    aft = rhsID[rhsID>otfID]
    if (len(bef)>0) & (len(aft)>0):
        rhIDbef = bef[np.argmax(bef)]
        rhIDaft = aft[np.argmin(aft)]
        rhIDs = [rhIDbef, rhIDaft]
    elif (len(aft)>0):
        rhIDaft = aft[np.argmin(aft)]
        rhIDs = [rhIDaft]
        print('getCalSpectra: Only REFHOT after OTF available')
    elif (len(bef)>0):
        rhIDbef = bef[np.argmax(bef)]
        rhIDs = [rhIDbef]
        print('getCalSpectra: Only REFHOT before OTF available')
    else:
        print('getCalSpectra: Not enough REFHOT scans available (before/after OTF): ', bef, aft)
        return -999, 0, 0, 0, 0, 0, 0, [0,0], 0

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
        osel = np.argwhere((mixer == mixers) & (scan_type == 'OTF') & (rflags) & (abs(stime - rtime) < 900))
        closest = np.argmin(np.abs(stime[osel] - rtime))
        ohsel = np.argwhere((mixer == mixers) & (scan_type == 'HOT') & (rflags) & (abs(stime - htime) < 900))
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
    The function returns 4 arrays: hgroup, ghots, ghtim, and glast
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
    if uots.size == 0 or uots.size > 1:
        # bad sequence with none or too many OTFs
        print("SEQUENCE HAS THE WRONG NUMBER OF OTFs")
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
    
    rowflagvalue = 255
    rowflagvalue = [80, 255, 4096, 4]
    
    # enum.show_flag_values(rowflagvalue) works only on single values
    flag_values = [enum.show_flag_values(rfv) for rfv in rowflagvalue]
    print('rowflagvalue input: ', rowflagvalue, 'contained flags: ', flag_values)
    print()
    
    print('checkRowflag: ', checkRowflag(rowflagvalue))


    """
    if type(rowflagvalue) == type(0):
        rowflagvalue = np.array([rowflagvalue])
    negfilt = ~rowflagfilter#.__invert__()
    return np.array([False if int(negfilt).__and__(int(rfvalue)) else True for rfvalue in rowflagvalue], dtype=bool)


checkRowflags = np.vectorize(checkRowflag)


def string_to_enum_combination(istring):
    """Converts a string of color names to an enum combination.
    
    """
    if type(istring)==type('m'):
        inames = istring.split()
    else:
        inames = str(istring).split()
    print(inames)
    cflags = RowFlags(0)  # Initialize with no value
    for name in inames:
        if '|' in name:
            continue
        else:
            try:
                flag = RowFlags[name.replace("'","").replace('"','').split('.')[1].upper()]
                cflags |= flag
            except KeyError:
                raise ValueError(f"Invalid flag: {name}")
    return cflags


def anaFlagString(aa):
    r"""Function .

    Parameters
    ----------
    aa : int or string
            rowflagfilter integer or string representation of allowed flags

    Examples
    --------
    aa = 'RowFlags.NO_HK | RowFlags.MISSING_INT | RowFlags.MIXER_UNPUMPED'
    enum_combination = anaFlagString(aa)
    
    print(repr(enum_combination), type(enum_combination))
    #<RowFlags.NO_HK|MISSING_INT|MIXER_UNPUMPED: 84> <flag 'RowFlags'>
    
    """
    if aa.isnumeric():
        return RowFlags(int(aa))
    else:
        aa = aa.replace("|","")
        return string_to_enum_combination(aa)

