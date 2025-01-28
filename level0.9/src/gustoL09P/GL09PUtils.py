"""
GUSTO Pipeline utility functions

V. Tolls, SAO

created: 9/19/2024
"""

import numpy as np
from astropy import units as u
import astropy.constants as const
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS, WCSSUB_LONGITUDE, WCSSUB_LATITUDE, WCSSUB_SPECTRAL
from astropy.wcs import validate as WCS_validate
from astropy.io import ascii, fits
from astropy.utils.exceptions import AstropyWarning
from scipy.interpolate import krogh_interpolate, barycentric_interpolate
from pprint import pprint
import parsl
from parsl.app.app import python_app, bash_app
from parsl.configs.local_threads import config
from parsl.providers import LocalProvider
from parsl.channels import LocalChannel
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
import sys
import os
import re
import inspect
import datetime
import warnings
import pkg_resources
from astropy.io import fits
from tqdm import tqdm
import logging
log = logging.getLogger(__name__)


__version__ = 0.11
__date__ = '20240919'
__updated__ = '20240919'
fheader, tail = os.path.split(inspect.stack()[0][1])
__pyfile__ = tail
__path__ = fheader


# warnings.filterwarnings('ignore', category=Warning, message=' FITSFixedWarning: The WCS transformation has more axes (3) than the image it is associated with (2) [astropy.wcs.wcs]', append=True)
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
    """Function converting str data to numpy float array.
        
    Parameters
    ----------
    icpar: str
        parameter str like '[start, stop, step]' or '3.45'
    """
    icpars = icpar.replace('[','').replace(']','').replace(' ','').split(',')
    ipars = np.array((icpars), dtype=dtype)
    
    return ipars


def simpleDespikeSpectra(spec0, data, hdr, pwidth=10, verbose=False, interpflag=False):
    """Function determining spikes in data array and save the
    spike information in the channel mask / spectra mask.


    Parameters
    ----------
    param1 : int

    Returns
    -------

    """
    if verbose:
        print('Started despiking function.')
    n_spec, n_pix = spec0.shape
    spec = spec0.copy()

    # there is the channel mask already in the table for flagging pixels
    ch_mask = data['CHANNEL_FLAG']
    row_flag = data['ROW_FLAG']
    
    psz = pfr_ra.shape
    if verbose:
        print(psz)

    # this flagging is according to pixel position
    # then replace flux with interpolated value
    for i in range(pp_ra.shape[0]):
        ch_mask[:, pp_ra[i,0]:pp_ra[i,1]] = 1
        if pp_ra[i,0]<pp_ra[i,1]:
            apix = np.arange(pp_ra[i,0]-pwidth, pp_ra[i,1]+pwidth+1, 1)
            if interpflag:
                pargs = apix[(apix<pp_ra[i,0])|(apix>pp_ra[i,1])]
                margs = apix[(apix>=pp_ra[i,0])&(apix<=pp_ra[i,1])]
                for k in trange(0,n_spec):
                    spec[k,pp_ra[i,0]:pp_ra[i,1]+1] = np.interp(margs, pargs, spec[i,pargs])
                    # if np.any(np.isfinite(spec[k,:]))==False:
                    #     print(i, k, spec[k,pp_ra[i,0]-pwidth, pp_ra[i,1]+pwidth+1])

    # for k in trange(0,n_spec):
    #     print(i, k, np.mean(spec[k,:]), np.mean(spec[k,:]))
        
    # try to determine "good" spectra
    # args = range(300, 400)
    # mask = np.ones(n_spec)
    
    # for i in trange(0,n_spec):
    #     sp = spec[i,:]
    #     if np.any(np.isfinite(spec[i,args])):
    #         if np.nanmean(spec[i,args]) > 100.0:
    #             var[i] = (np.nanmax(spec[i,args] - np.nanmin(spec[i,args]))) / np.abs(np.nanmean(spec[i,args]))
    #             if (var[i] < vcut) & (var[i] > 0.0):
    #                 mask[i] = 0    
    #     else:
    #         # there is no valid spectrum
    #         print(i, np.mean(spec[i,args]), spec[i,args])
    #         row_flag[i] = 2

    if verbose:
        print('Done despiking.')
    return spec


def getSpecScanTypes(mixer, spec, data, hdr, verbose=False):
    """Function calculating the calibration spectrum for a single mixer.


    Parameters
    ----------
    mixer : int
        mixer number
    spec : masked float array
        array containing the latest masked spectra
    data : FITS rec array
        record array containing the GUSTO data fron the FITS file

    Returns
    -------

    """

    mixers  = data['MIXER']
    scanID = data['scanID']
    # ra     = data['RA']
    # dec    = data['DEC']
    THOT   = data['THOT']
    scan_type = data['scan_type']
    row_flag = data['ROW_FLAG']      # spectra mask (one entry per spectrum)
    ch_flag = data['CHANNEL_FLAG']   # spectral pixel (or channel) mask

    rfsel = (mixer == mixers) & (scan_type == 'REF') & (row_flag==0)
    rfsID = np.unique(scanID[rfsel])
    rhsel = (mixer == mixers) & (scan_type == 'REFHOT') & (row_flag==0)
    rhsID = np.unique(scanID[rhsel])
    rhsel = (mixer == mixers) & (scan_type == 'REFHOT') & (row_flag==0)
    rhsID = np.unique(scanID[rhsel])

    # identify the HOT spectra
    htsel = (mixer == mixers) & (scan_type == 'HOT') & (row_flag==0)
    hotID = np.unique(scanID[htsel])

    # identify the OTF spectra
    otsel = (mixer == mixers) & (scan_type == 'OTF') & (row_flag==0)
    otfID = np.unique(scanID[otsel])

    # if verbose:
    #     print('Mixer: ', mixer)
    #     print('REF scan IDs: ', rfsID)
    #     print('REFHOT scan IDs: ', rhsID)
    #     print('HOT scan IDs: ', hotID)
    #     print('OTF scan IDs: ', otfID)

    return otfID, rfsID, rhsID, hotID


def getCalSpectra(mixer, spec, data, hdr, Tsky=45., verbose=False):
    """Function calculating the calibration spectrum for a single mixer.


    Parameters
    ----------
    mixer : int
        mixer number
    spec : masked float array
        array containing the latest masked spectra
    data : FITS rec array
        record array containing the GUSTO data fron the FITS file
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
    THOT   = data['THOT']
    scan_type = data['scan_type']
    row_flag = data['ROW_FLAG']      # spectra mask (one entry per spectrum)
    ch_flag = data['CHANNEL_FLAG']   # spectral pixel (or channel) mask
    stime = data['UNIXTIME']

    otfID, rfsID, rhsID, hotID = getSpecScanTypes(mixer, spec, data, hdr, verbose=verbose)
    if (len(otfID)<1) | (len(rfsID)<1) | (len(rhsID)<1) | (len(hotID)<1):
        print('Not enough scans types for processing (otf/refs/refhots/hots): ', otfID, rfsID, rhsID, hotID)
        return -999, 0, 0, 0, 0, 0, 0, [0,0]
    
    # determine the REFHOTs that bracket the OTFs
    bef = rhsID[rhsID<otfID]
    aft = rhsID[rhsID>otfID]
    if (len(bef)>0) & (len(aft)>0):
        rhIDbef = bef[np.argmax(bef)]
        rhIDaft = aft[np.argmin(aft)]
    else:
        print('Not enough ref scans available (REFs before/after OTF): ', bef, aft)
        return -999, 0, 0, 0, 0, 0, 0, [0,0]

    rhIDs = [rhIDbef, rhIDaft]

    Tsyss = []
    REFs = []
    RHOTs = []
    ttimes = []
    rtimes = []
    htimes = []
    Thot = []
    for rhID in rhIDs:
        # determine yfactor
        rsel = (rhID == scanID) & (mixer == mixers) & (scan_type == 'REF') & (row_flag==0)
        hsel = np.argwhere((rhID == scanID) & (mixer == mixers) & (scan_type == 'REFHOT') & (row_flag==0))
        # if we have 3 (or more, more than3 is less likely) take only the last 2 since the first one is causing offsets!
        if hsel.size>2:
            data['ROW_FLAG'][hsel[0]] = 2<<15
            hsel = hsel[1:]
        #print(rselbef)
        
        spec_h =  spec[hsel,:].sum(axis=0)/len(spec[hsel,:])
        spec_r =  spec[rsel,:].sum(axis=0)/len(spec[rsel,:])
        htime = stime[hsel].mean()
        rtime = stime[rsel].mean()
        #print(list(spec[,:]))
        
        THOT_avg = THOT[hsel].sum(axis=0)/len(THOT[hsel])
        
        # estimate Tsys for each Device
        y_factor  = spec_h/spec_r

        tsys = np.squeeze((THOT_avg - Tsky*y_factor[:])/(y_factor[:] - 1.))

        Thot.append(THOT_avg)
        Tsyss.append(tsys)
        REFs.append(spec_r)
        RHOTs.append(spec_h)
        rtimes.append(rtime)
        htimes.append(htime)

    Tsyss = np.ma.array(Tsyss)
    REFs = np.ma.array(REFs)
    RHOTs = np.ma.array(RHOTs)
    return Tsyss, REFs, RHOTs, rtimes, htimes, Thot, Tsky, rhIDs


def getHotInfo(spec, data, mixer, dfile='', verbose=False):
    """Function analyzing and processing HOT spectra in an 
    OTF strip.
    There are caveats in the current data (tbc) including that there are duplicate
    REFs at the end of the sequence. This issue has been mitigated for now by
    ignoring these REFs. The current assumption, yet unverified, is that there is 
    always a REF at the beginning of the sequence.


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
    glast: flag if there is a HOT after the last OTF spectrum
    htflag: flag indicating the availability of REFs and REFHOTs
        0: clean sequence with REFs at start and end
        1:
        2: double REFs at start
        3: double REFs at end
        4: REF at start
        5: REF at end
        -1: no or too man OTFs
        -2: too many OTF scan IDs
        -3: 2 REFs at start or end; reverting to method 1
        -4: 

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
    
    # determine the first hot scan
    lasthot = np.argmin(unixtime[(data['MIXER']==umixers[mx])&(data['scan_type']=='HOT')])
    firstref = True if (np.argwhere((data['MIXER']==umixers[mx])&(data['scan_type']=='REF')).min() < 10) else False

    # added check for REFHOT duplicates at beginning and at end of sequence
    # there should only be one REFHOT at the beginning and one at
    # the end
    rhscans = data['scanID'][(data['scan_type']=='REFHOT')&(data['ROW_FLAG']==0)&(data['MIXER']==mixer)]
    rfscans = data['scanID'][(data['scan_type']=='REF')&(data['ROW_FLAG']==0)&(data['MIXER']==mixer)]
    otscans = data['scanID'][(data['scan_type']=='OTF')&(data['ROW_FLAG']==0)&(data['MIXER']==mixer)]
    urhs = np.unique(rhscans)
    urfs = np.unique(rfscans)
    uots = np.unique(otscans)
    htflag = 0

    # print('urhs: ', urhs)
    # print('urfs: ', urfs)
    # print('uots: ', uots)
    # print()
    
    # the first OTF scan ID should be the one determining the sequence
    if uots.size == 0:
        # bad sequence with no OTFs
        htflag = -1
        return 0,0,0,0,0, -1
    elif uots.size == 1:
        htflag = 0
        tscan = uots[0]
    elif uots.size > 1:
        # bad sequence with too many OTFs
        htflag = -2
        return 0,0,0,0,0, -2

    # print('htflag: ', htflag)
    
    # case 2 REFHOTs
    if urfs.size==1:
        # print('size: 1\n')
        if (urfs[0]<tscan):
            # REF at beginning
            htflag = 4
            good_rhs = urhs
        elif (urfs[0]>tscan):
            htflag = 5
            good_rhs = urhs
    elif urfs.size==2:
        # print('size: 2\n')
        # check if on rfscan is smaller and one larger than tscan
        if (urfs[0]<tscan) & (urfs[1]>tscan):
            # clean sequence REF OTF REF
            htflag = 0
            good_rhs = urhs
        else:
            # start double REFs
            # end double REFs
            # either case should revert to cal method 1!
            good_rhs = urhs
            htflag = -3
            return 0,0,0,0,0, -1
    elif urfs.size==3:
        # print('size: 3\n')
        # print('tscan: ', tscan)
        # print()
        # print('tscan < urfs[2]: ', tscan < urfs[2])
        # print('tscan > urfs[1]: ', tscan > urfs[1])
        # print()
        # now we duplicate REFs
        if (tscan < urfs[1]):
            # double at end
            htflag = 3
            good_rhs = urhs[:2]
        elif tscan > urfs[1]:
            # double at beginning
            htflag = 2
            good_rhs = urhs[1:3]
        elif tscan < urfs[2]:
            # double at beginning, 
            htflag = 2
            good_rhs = urhs[1:3]
    elif urfs.size>3:
        # print('size: >3\n')
        if tscan > urfs[1]:
            # double at beginning
            htflag = 2
            good_rhs = urhs[1:3]
        elif tscan < urfs[1]:
            # double at end
            htflag = 3
            good_rhs = urhs[:2]
        else:
            # totally wrong sequence
            htflag = -4
            return 0,0,0,0,0, -1
            
    # print()
    # print('htflag: ', htflag)
    # print('good_rhs: ', good_rhs)
    uflag = np.zeros(n_spec)
    
    for i in range(n_spec):
        if data['MIXER'][i] == mixer:
            if (data['scan_type'][i] == 'HOT')|((data['scan_type'][i] == 'REFHOT')&(data['scanID'][i] in good_rhs)):
                # print(hcnt, data['scanID'][i], good_rhs, (data['scanID'][i] in good_rhs), (data['scan_type'][i] == 'HOT'))
                if hcnt==0:
                    #hgroup[i] = hgrp
                    # if unixtime[i] - unixtime[lasthot] < 4.0:
                    #     hgroup[i] = hgrp
                    # else:
                    #     hgrp += 1
                    #     hgroup[i] = hgrp
                    lasthot = i
                else:
                    if unixtime[i] - unixtime[lasthot] < 4.0:
                        pass
                        #hgroup[i] = hgrp
                    else:
                        hgrp += 1
                        #hgroup[i] = hgrp
                    lasthot = i
                uflag[i] = 1
                hcnt += 1
            else:
                uflag[i] = 0
            hgroup[i] = hgrp
            if ((data['scan_type'][i] == 'REF')&(data['scanID'][i] in good_rhs)):
                uflag[i] = 1
            
                    
            # if verbose:
            #     if i == 0:
            #         #     #  0   8897    0    0    0     REF  5  0          0.000     1.833
            #         print('#sp scanID hcnt hgrp hgrp huse scantyp Mx rf           time   inttime')
            #     print('%3i %6i  %3i %4i %4i %4i  %6s  %i  %i  %13.3f  %8.3f'%(i, data['scanID'][i], hcnt, hgroup[i], hgrp, uflag[i], data['scan_type'][i], data['MIXER'][i], data['ROW_FLAG'][i], unixtime[i]-ut0, data['INTTIME'][i]))
    
    
    
    # maxgrp = np.zeros(umixers.size, dtype=int)
    # ghots = np.zeros((int(umixers.size), int(hgroup.max()+1), n_pix))
    # ghtim = np.zeros((int(umixers.size), int(hgroup.max()+1)))
    # ghtint = np.zeros((int(umixers.size), int(hgroup.max()+1)))
    # glast = np.zeros(int(umixers.size), dtype=bool)
    #
    # # now we have to average the hots for each group:
    # # for i, mx in enumerate(umixers):
    # for i, mx in enumerate([mixer]):
    #     maxgrp[i] = int(hgroup[data['MIXER']==mx].max()+1)
    #     if verbose:
    #         print('mixer/# groups: ', mx, maxgrp[i])
    #
    #     for j in range(maxgrp[i]):
    #         sel = (data['MIXER']==mx) & (hgroup==j) & (data['ROW_FLAG']==0)
    #         ghots[i,j,:] = np.mean(spec[sel,:], axis=0)
    #         ghtim[i,j] = np.mean(unixtime[sel])
    #         ghtint[i,j] = np.sum(tint[sel])
    #
    #     # final check if the last OTF is followed by a HOT/REFHOT
    #     sel = np.argwhere((data['MIXER']==mx) & (data['scan_type']=='OTF') & (data['ROW_FLAG']==0))
    #     if np.max(data['UNIXTIME'][sel])<ghtim[i,-1]:
    #         glast[i] = True
    maxgrp = 0
    ghots = np.zeros((int(hgroup.max()+1), n_pix))
    ghtim = np.zeros(int(hgroup.max()+1))
    ghtint = np.zeros(int(hgroup.max()+1))
    glast = False
    
    # now we have to average the hots for each group flagged for use:
    # for i, mx in enumerate(umixers):
    for i, mx in enumerate([mixer]):
        maxgrp = int(hgroup[data['MIXER']==mx].max()+1)
        # if verbose:
        #     print('mixer/# groups: ', mx, maxgrp)
    
        for j in range(maxgrp):
            sel = (data['MIXER']==mx) & (hgroup==j) & (data['ROW_FLAG']==0) & (uflag==1)
            ghots[j,:] = np.mean(spec[sel,:], axis=0)
            ghtim[j] = np.mean(unixtime[sel])
            ghtint[j] = np.sum(tint[sel])
    
        # final check if the last OTF is followed by a HOT/REFHOT
        sel = np.argwhere((data['MIXER']==mx) & (data['scan_type']=='OTF') & (data['ROW_FLAG']==0))
        if np.max(data['UNIXTIME'][sel])<ghtim[-1]:
            glast = True
        
    # if verbose:
    #     print('mx: ', mixer, '  good_rhs: ', good_rhs, '  firstref: ', firstref, '  glast: ', glast, '  #hgrp: ', hgroup.max()+1, maxgrp, dfile,'\n')

    return hgroup, ghots, ghtim, ghtint, glast, htflag


# old version of getHotInfo - not working correctly
def getHotInfo_old(spec, data, mixer, verbose=False):
    """Function analyzing and processing HOT spectra in an 
    OTF strip.


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
    glast: flag if there is a HOT after the last OTF spectrum

    """
    n_spec, n_pix = spec.shape
    
    # mixer index for testing: 0, 1, or 2
    #mx = mixer
    umixers = np.unique(data['MIXER'])
    mx = np.argwhere(umixers==mixer).flatten()
    
    unixtime = data['UNIXTIME']
    tint = data['INTTIME']
    ut0 = unixtime[0]
    
    hgroup = np.zeros(n_spec, dtype=int)
    hcnt = 0   # counter counting the total number of hots
    hgrp = 0   # counter for hot groups
    
    # determine the first hot scan
    lasthot = np.argmin(unixtime[(data['MIXER']==umixers[mx])&(data['scan_type']=='HOT')])
    
    for i in range(n_spec):
        for mx in umixers:
            if data['MIXER'][i] == mx:
                if (data['scan_type'][i] == 'HOT')|(data['scan_type'][i] == 'REFHOT'):
                    if hcnt==0:
                        #hgroup[i] = hgrp
                        # if unixtime[i] - unixtime[lasthot] < 4.0:
                        #     hgroup[i] = hgrp
                        # else:
                        #     hgrp += 1
                        #     hgroup[i] = hgrp
                        lasthot = i
                    else:
                        if unixtime[i] - unixtime[lasthot] < 4.0:
                            pass
                            #hgroup[i] = hgrp
                        else:
                            hgrp += 1
                            #hgroup[i] = hgrp
                        lasthot = i
                    hcnt += 1
                hgroup[i] = hgrp
                        
                if verbose:
                    # if i == 0:
                    #     #     #  0   8897    0    0    0     REF  5  0          0.000     1.833
                    #     print('#sp scanID hcnt hgrp hgrp scantyp Mx rf           time   inttime')
                    # print('%3i %6i  %3i %4i %4i  %6s  %i  %i  %13.3f  %8.3f'%(i, data['scanID'][i], hcnt, hgroup[i], hgrp, data['scan_type'][i], data['MIXER'][i], data['ROW_FLAG'][i], unixtime[i]-ut0, data['INTTIME'][i]))
                    pass
    
    if verbose:
        print('\nNumber of hot groups per mixer: ', hgroup.max()+1,'\n')
    
    # maxgrp = np.zeros(umixers.size, dtype=int)
    # ghots = np.zeros((int(umixers.size), int(hgroup.max()+1), n_pix))
    # ghtim = np.zeros((int(umixers.size), int(hgroup.max()+1)))
    # ghtint = np.zeros((int(umixers.size), int(hgroup.max()+1)))
    # glast = np.zeros(int(umixers.size), dtype=bool)
    #
    # # now we have to average the hots for each group:
    # # for i, mx in enumerate(umixers):
    # for i, mx in enumerate([mixer]):
    #     maxgrp[i] = int(hgroup[data['MIXER']==mx].max()+1)
    #     if verbose:
    #         print('mixer/# groups: ', mx, maxgrp[i])
    #
    #     for j in range(maxgrp[i]):
    #         sel = (data['MIXER']==mx) & (hgroup==j) & (data['ROW_FLAG']==0)
    #         ghots[i,j,:] = np.mean(spec[sel,:], axis=0)
    #         ghtim[i,j] = np.mean(unixtime[sel])
    #         ghtint[i,j] = np.sum(tint[sel])
    #
    #     # final check if the last OTF is followed by a HOT/REFHOT
    #     sel = np.argwhere((data['MIXER']==mx) & (data['scan_type']=='OTF') & (data['ROW_FLAG']==0))
    #     if np.max(data['UNIXTIME'][sel])<ghtim[i,-1]:
    #         glast[i] = True
    maxgrp = 0
    ghots = np.zeros((int(hgroup.max()+1), n_pix))
    ghtim = np.zeros(int(hgroup.max()+1))
    ghtint = np.zeros(int(hgroup.max()+1))
    glast = False
    
    # now we have to average the hots for each group:
    # for i, mx in enumerate(umixers):
    for i, mx in enumerate([mixer]):
        maxgrp = int(hgroup[data['MIXER']==mx].max()+1)
        if verbose:
            print('mixer/# groups: ', mx, maxgrp)
    
        for j in range(maxgrp):
            sel = (data['MIXER']==mx) & (hgroup==j) & (data['ROW_FLAG']==0)
            ghots[j,:] = np.mean(spec[sel,:], axis=0)
            ghtim[j] = np.mean(unixtime[sel])
            ghtint[j] = np.sum(tint[sel])
    
        # final check if the last OTF is followed by a HOT/REFHOT
        sel = np.argwhere((data['MIXER']==mx) & (data['scan_type']=='OTF') & (data['ROW_FLAG']==0))
        if np.max(data['UNIXTIME'][sel])<ghtim[-1]:
            glast = True
        
    if verbose:
        print('\nglast: ', glast)

    return hgroup, ghots, ghtim, ghtint, glast


    