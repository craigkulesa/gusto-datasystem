#!/usr/bin/env python3.11 

import matplotlib.pyplot as plt
import matplotlib as mpl

import numpy as np

import pandas as pd
import os
import sys
import glob

from astropy.io import fits
from astropy.table import Table
from astropy.table import vstack
from astropy import units as u
from astropy import constants as const
from scipy import interpolate
from scipy import signal
#from sklearn.neighbors import LocalOutlierFactor
#from scipy import ndimage as ndi
#from skimage.util import random_noise

#import skimage as ski
from multiprocessing import Pool
from functools import partial

from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic

#from pybaselines import Baseline, utils
#from PyAstronomy import pyasl

from flagdefs import *

import math

from BayesicFitting import GaussModel
from BayesicFitting import PolynomialModel
from BayesicFitting import ConstantModel
from BayesicFitting import NestedSampler
from BayesicFitting import formatter as fmt
#from BayesicFitting import plotFit
#from BayesicFitting import SineModel
#from BayesicFitting import ExpModel
#from BayesicFitting import ExponentialPrior
from BayesicFitting import LevenbergMarquardtFitter
from BayesicFitting import RobustShell
from BayesicFitting import Tools




GUSTO_DIR='/data/scratch/GUSTO/gusto-datasystem/'

L9DATA = f'{GUSTO_DIR}/Data/level0.9/'
L8DATA = f'{GUSTO_DIR}/Data/level0.8/'
L7DATA = f'{GUSTO_DIR}/Data/level0.7/'
L10DATA= f'{GUSTO_DIR}/Data/level1/'
#obj = 'G348'
obj = sys.argv[1]
b = int(sys.argv[2])

if b == 2:
    line = 'CII'
else:
    line = 'NII'

print('OBJECT', obj, 'LINE',line, b*1)
#L10DATA= f'{GUSTO_DIR}/Data/level1/G337/'


def idspike(x,data,cflags,start,stop,points=10,count=3,deg=2,dx=1,stdlim=0.05):
    mask = np.zeros(len(data),dtype=bool)
    pm = PolynomialModel(deg)
    model = pm
    lmf = LevenbergMarquardtFitter(x[start:stop],model)
    ftr = RobustShell( lmf )
    par = ftr.fit( data[start:stop],verbose=0)
    rwgt = ftr.wieghts
    qmask = rwgt < stdlim
    cflags[qmask] |= ChanFlags.SPUR
    newdata = np.ma.masked_array(data,qmask)
    return newdata, cflags

def bayesspike(x,data,cflags,deg=1,spurlim = 0.5):

    pm = PolynomialModel( deg )
    model = pm
    lmf = LevenbergMarquardtFitter( x, model )
    ftr = RobustShell( lmf )

    pars = ftr.fit( data, verbose=0 )
    #print( "params  ", fmt( pars ) )
    #print( "stdevs  ", fmt( ftr.stdevs ) )
    #print( "scale   ", fmt( ftr.scale ) )
    #lo = [ 0.9, 0.9, -0.1, -0.1]
    #hi = [ 1.1, 1.1, +0.1, +0.1]
    newchf = cflags.copy()

    rwgt = ftr.weights  # weights of points ignored are ~0
    qmask = (rwgt < spurlim) & (cflags == 0)
    #make an exclusive OR with existing flags
    newchf[qmask] |= ChanFlags.SPUR_CANDIDATE
    #update input data
    newdata = np.ma.masked_array(data,qmask)
    
    return newdata, rwgt, newchf, model, qmask

def TestSPUR(x,y,rms):
    mdl = GaussModel( )
    #print( mdl )
    lo = [rms,x.min(),1]
    hi = [30,x.max(),4]
    mdl.setLimits( lo, hi )
    ns = NestedSampler( x, mdl, y )
    ns.distribution.setLimits( [0.1, 10] )
    evid = ns.sample( plot=False )
    sl = ns.samples
    par = sl.parameters
    print( "Parameters :", fmt( ns.parameters ) )
    #print( "StDevs     :", fmt( ns.stdevs ) )
    #print( "Scale      :", fmt( ns.scale ) )
    print( "Evidence   :", fmt( evid ) )
    # Test no line
    mdl0 = ConstantModel()
    ns0 = NestedSampler( x, mdl0, y )
    ns0.distribution.setLimits( [0.1,10] )
    evid0 = ns0.sample( plot=False )
    sl0 = ns0.samples
    spurfound = evid > evid0
    return evid, par, mdl, evid0




def ticks_IF(x,pos):
    #print(x)
    delIF = IF_freq[1]-IF_freq[0]
    zeroIF = IF_freq[0]
    #print(delIF,zeroIF)
    return f'{x*delIF - zeroIF:4.0f}'

def mark_spurs(rwgt,chanflag,ax,threshold=0.2):
    newblob = []

    qch = np.where(chanflag)
    qarg = np.argwhere(chanflag)
    for q in qarg:
        #print(q.shape)
        x = q[1]
        y = q[0]
        r = 1
        
        if x < 300: # & (rwgt[x,y] < threshold):
            cpatch = plt.Circle((x,y),r,color='r',linewidth=1 ,fill=False)
            newblob.append([x,y,rwgt[y,x]])
            ax.add_patch(cpatch)
    
    return( newblob ,ax )

#def find_SPUR(data1,hdr1,hdr0,vlsr):
def find_SPUR(file):

    print(file)
    hdu = fits.open(file)
    hdr0 = hdu[0].header
    hdr1 = hdu[1].header
    data1 = hdu[1].data
    # prepare velocity/IF or Channel axis
    npix  = hdr0['NPIX']
    chanarr = np.arange(npix)
    VLSR = hdr0['VLSR']
    IF_vlsr0  = hdr0['IF0']
    line_freq = hdr0['linefreq']

    def to_vlsr(x):
        c = const.c.value/1e3
        x = (IF_vlsr0 - x )/line_freq * c + VLSR # Vlsr in km/s 
        return(x) 
    
    def to_IF(x): 
        c = const.c.value/1e3
        x = IF_vlsr0 - ( x - VLSR) * line_freq / c
        near_zero = np.isclose(x, 0)
        x[near_zero] = 0 
        return(x) 

    if hdr0['DLEVEL'] < 1:
        IF_pix = hdr0['CRPIX1']
        IF_val = hdr0['CRVAL1']
        IF_del = hdr0['CDELT1']
        IF_freq = (chanarr-IF_pix)*IF_del+IF_val
    
        vlsr = (IF_vlsr0 - IF_freq)/line_freq*const.c.value/1.e3 +VLSR # Vlsr in km/s
        freq = IF_freq + hdu[0].header['SYNTFREQ'] * hdu[0].header['SYNTMULT']
        cdelt = vlsr[1]-vlsr[0]
    else:
        crpix1 = hdr0['CRPIX1']
        crval = hdr0['CRVAL1']
        cdelt = hdr0['CDELT1']
        vlsr = (chanarr - crpix1)*cdelt + crval
        IF_freq = to_IF(vlsr)
    hdu.close()
    V0 = vlsr[0]

    def channel_to_vlsr(x):
        x = x * cdelt + V0 
        return (x)
    
    def vlsr_to_channel(x):
        x = (x - V0) / cdelt
        return (x)



    #if line == 'CII':
    #    fig, ((ax1,ax1a),(ax2,ax2a)) = plt.subplots(2,2,figsize=(16,9),layout= 'compressed')
    #    ax1.set_ylim(-15,40)
    #    ax2.set_ylim(-15,40)
    #    #secax = ax1.secondary_xaxis('top', functions = (to_IF, to_vlsr))
    #    secax = ax1.secondary_xaxis('top', functions = (vlsr_to_channel, channel_to_vlsr))
    #    
    #    ax2.set_xlabel('$V_{LSR} [km/s]$')
    #    #secax.set_xlabel('IF (MHz)')
    #    secax.set_xlabel('Channel')
    #    ax1.set_ylabel('T (K)')
    #    ax2.set_ylabel('T (K)')
    #    xxs = [ax1,ax2]
    #    xxim = [ax1a,ax2a]
    #else:
    fig, ((ax1,ax1a),(ax2,ax2a),(ax3,ax3a)) = plt.subplots(3,2,figsize=(16,9),layout= 'compressed')
    ax1.set_ylim(-10,30)
    ax2.set_ylim(-10,30)
    ax3.set_ylim(-10,30)
    #secax = ax1.secondary_xaxis('top', functions = (to_IF, to_vlsr))
    secax = ax1.secondary_xaxis('top', functions = (vlsr_to_channel, channel_to_vlsr))
    
    ax3.set_xlabel('$V_{LSR} [km/s]$')
    #secax.set_xlabel('IF (MHz)')
    secax.set_xlabel('Channel')
    ax1.set_ylabel('T (K)')
    ax2.set_ylabel('T (K)')
    ax3.set_ylabel('T (K)')
    xxs = [ax1,ax2,ax3]
    xxim = [ax1a,ax2a,ax3a]



    newspurs = []
    mixer  = data1['MIXER']
    mixers = np.unique(mixer)
    #print(f'Found mixers: {mixers}')
    scanIDs = data1['scanID']
    rms = data1['rms']
    scanID = np.unique(scanIDs)[0]

    seqID = hdr0['SEQ_ID']
    OBJECT = hdr0['OBJECT']
    seqflag = hdr0['SEQ_FLAG']

    ra     = data1['RA']
    dec    = data1['DEC']
    
    THOT   = hdr0['THOT']
    rflag  = data1['ROW_FLAG']
    chflag = data1['CHANNEL_FLAG']
    rowflag= data1['ROW_FLAG']
    scan_type = data1['scan_type']
    # sky temperature Callen Welton hnu/2k
    line_freq = hdr0['LINEFREQ']
    Tsky = (const.h*line_freq*1e6*u.Hz /(2*const.k_B)).value



    # index based on channels
    qv = (chanarr > 40*b) & (chanarr < 300*b)
    #qv = (vlsr > -450) & (vlsr < 40)
    qvarg = np.argwhere(qv).flatten()

    # Iterate over mixers here    
    # qotf0 = (data1['scan_type'] == 'OTF') & (data1['mixer'] == mixers[m])
    for ix, mx in enumerate(mixers):
        ax = xxs[ix]
        rfl = RowFlags.RINGING_BIT0 | RowFlags.RINGING_BIT1 | RowFlags.MIXER_UNPUMPED 

        qotf0 = (data1['scan_type'] == 'OTF')&(data1['mixer'] == mx) & ((rowflag & rfl)==0)
        if qotf0.max():
            c_radec = SkyCoord(ra = ra[qotf0]*u.degree, dec = dec[qotf0]*u.degree, frame = 'icrs')
            c_gal = c_radec.transform_to(Galactic)
        
            
            #l=c_gal.l.value
            #b=c_gal.b.value
            medlat = np.median(c_gal.l.value)
    
            unixtime = data1[qotf0]['unixtime']
            ttime=unixtime-unixtime[0]
            norm= plt.Normalize(ttime[0],ttime[-1])
            #norm=plt.Normalize(b.min(),b.max())
            #norm=plt.Normalize(-1,1)
            otfspec = data1[qotf0]['data']
            mixspec = data1[qotf0]['mixer']
            badchan = data1[qotf0]['CHANNEL_FLAG']
            cf_update = badchan.copy()
            tsys = data1[qotf0].Tsys
            rms = data1[qotf0].rms
            mspec = np.ma.array(otfspec,mask = badchan)
             
            spurlim = 0.5
            rwgt = np.zeros(mspec.shape)
            newchf = np.zeros(mspec.shape,'>i2')

            
            # step through all otf spectra an create a weights map
            for ii, mspec1 in enumerate(mspec):
                mask = np.zeros(len(mspec1),dtype=bool)
                ttime1 = ttime[ii]
                rms1 = rms[ii]
                cflags1 = badchan[ii]

                x = vlsr[qv]
                y = mspec1[qv]
                cflg1 = cflags1[qv]

                # create x   and y vectors of only usable points 
                qkeep = np.invert(y.mask)
                qkarg = np.argwhere(qkeep)
                xx = x[qkeep]
                yy = y.data[qkeep]
                #Don't consider data with Existing channel flags
                cf = cflg1[qkeep]

                #evidspur, spurpar, mdl, evid0 = TestSPUR(xx,yy,5*rms1)
                yyout, rwgt1, cfl1, mdl, qspur = bayesspike(xx, yy, cf, deg = 3, spurlim = spurlim)
                #put results back into original time order
                rwgt[ii,qvarg[qkeep]] = rwgt1
                newchf[ii,qvarg[qkeep]] = cfl1

                ax.plot(vlsr[qv],mspec1[qv],alpha=0.3, color=cmap(norm(ttime1)) )
                #IF new spurs are found, mark them on the plot
                if qspur.max():
                    ax.scatter( xx[qspur], yy[qspur], alpha=1, marker='+', color='r' )
    
            
            ax.annotate(f'Mixer {mx}' ,(vlsr[qvarg[qkeep]].min(),25))
            ax.annotate(f'Lat {medlat:0.3f}',(vlsr[qvarg[qkeep]].min(),20))
            # Now process weights image
            axa = xxim[ix]
            mxspurs = mark_spurs(rwgt,newchf,axa,threshold = spurlim)
            axa.imshow(1-rwgt) 
            #
            # Create the output table for the mixer
            #
            syntfreq= hdr0["SYNTFREQ"]
            syntmult = hdr0['SYNTMULT']
            try:
                syntemp = hdr0['B1_SYNTH']
            except:
                # no temperature present in header
                syntemp = -1000
            
            arr = np.array(mxspurs[0])
            if arr.any():
                t = Table(arr, names = ('Channel','T_idx','weight'), dtype=('i4','i4','f8')) #,dtype=[('x','i4'),('y','i4'),('w','f8')])
                test = [ut for ut in t['T_idx']]
                t.add_column(unixtime[test],name='UNIXT')
                t.add_column(seqID, name='SEQ') #,dtype='i4')
                t.add_column(seqflag, name='SEQ FLAG')
                t.add_column(scanID, name='SCAN')
                t.add_column(mx,name='MIXER')
                t.add_column(IF_vlsr0,name='IF0')
                t.add_column(VLSR,name='VLSR')
                t.add_column(syntfreq,name='SYNTFREQ')
                t.add_column(syntmult,name='SYNTMULT')
                t.add_column(syntemp,name='B1_SYNTH')
            
                if len(newspurs) == 0:
                    newspurs = t
                else:
                    newspurs = vstack([newspurs,t])
    
                   
    ptitle = f'SPUR_CANDIDATES: {line} {OBJECT} Seq:{seqID:05}  Scan:{scanID:06} ' 
    
    plt.suptitle(ptitle)
    fnameroot = f'{line}_{seqID:05}_{scanID:6}_SPURs'
    if len(newspurs) > 0:
        print(f'Saving found SPURS in {fnameroot}')
        #newspurs.write(f'/home/russ/GUSTO/scripts/{line}/{fnameroot}.fits',overwrite=True)
        #plt.savefig(f'/home/russ/GUSTO/scripts/{line}/{fnameroot}.png')
        newspurs.write(f'{L10DATA}/{obj}/{line}/{fnameroot}.fits',overwrite=True)
        plt.savefig(f'{L10DATA}/{obj}/{line}/{fnameroot}.png')
    else:
        print(f'No spurs found for {fnameroot}')
    #plt.show()
    plt.close()

files = glob.glob(f'{L10DATA}/{obj}/{line}*fits')
files.sort()


cmap = mpl.colormaps.get_cmap('BuGn')
cmap =mpl.colormaps.get_cmap('hot')
cmap = mpl.colormaps.get_cmap('viridis')
print(f'processing data for {obj}')

#for file in files:
if __name__ == "__main__":



    #find_SPUR(data1, hdr1, hdr0, vlsr) 
    #find_SPUR(file)
    cpus = 10
    
    if cpus > 1:
        pool = Pool(processes = cpus)
    else:
        pool = Pool()

    pool.map(partial(find_SPUR),files)
    

