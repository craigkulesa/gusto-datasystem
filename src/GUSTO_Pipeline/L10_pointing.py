#!/usr/bin/env python
"""
This is the executable for the GUSTO L2 Pipeline.
"""
import glob
import numpy as np
import time
import sys
import os
from multiprocessing.pool import Pool
from importlib.resources import files
from astropy import constants as const
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time

from .DataIO import loadL08Data
from .Utils import *
from .Logger import *

offsetsfile0 = files('GUSTO_Pipeline') / 'calib/offsets.txt'

def L10_Pipeline(args, scanRange, verbose=False):
    """Function processing the Level 0.95 data injecting 
    coordinate corrections. 

    1: read in the data from file
    2: apply coordinate offset correction to each spectrum
    3: save back to files with only setting a flag about the baseline correction,
    but all other parameters are passed through

    Parameters
    ----------
    args : list 
            list of configuration parameters provided by pipeline framework
    scanRange : int array
            array with firat and last scan number

    """
    if args.debug ==True:
        print('\nExecuting debug mode.\n')
        logger = logging.getLogger()
        logger.setLevel(10)
        n_procs = 1
    elif args.cpus.isnumeric():
        n_procs = int(args.cpus)
    else:
        n_procs = multiprocessing.cpu_count()
        
    print('Number of cores used for processing: %i\n'%(n_procs))
    ignore = [0]
    
    # read offsets file
    offsetsfile = offsetsfile0
    if not os.path.exists(offsetsfile):
        print('offsetsfile: ', offsetsfile)
        print('No file with coordinate offsets available. Terminating pipeline run!')
        sys.exit(1)
        
    offs = np.genfromtxt(offsetsfile, delimiter='\t', skip_header=2, 
                     dtype=[('mxpix', 'U4'), ('az', '<f8'), ('el', '<f8'), ('comment', 'U16')])

    
    for band in args.band:
        if verbose:
            print('\nProcessing line: ', line)
        # identify the files for processing
        inDir = args.path + 'level0.9/'
        outDir = args.path + 'level1/'
        if int(band) == 2:
            filter = 'CII*.fits'
        else:
            filter = 'NII*.fits'
        
        sdirs = sorted(glob.glob(os.path.join(inDir,filter)))
        dsc = [int(os.path.split(sdir)[1].split('_')[2].split('.')[0]) for sdir in sdirs]
        
        sdirs.sort(key=lambda sdirs: dsc)
        
        dfiles = []
        for i,ds in enumerate(dsc):
            if (ds >= scanRange[0]) & (ds <= scanRange[1]) & (ds not in ignore):
                dfiles.append(sdirs[i])
                        
        n_ds = len(dfiles)
                    
        paramlist = [[a, b, c, d, e, f] for a in [band] for b in [inDir] for c in [outDir] for d in dfiles for e in [offs] for f in [args.debug]]
        if verbose:
            print('Number of data files: ', n_ds, len(sdirs))
            #print('Selected data files: ', dfiles)
        
        # setup multiprocessing loop here to process each file in list
        with Pool(n_procs) as pool:
            # execute tasks in order
            for result in pool.imap(processL10, paramlist):
                print(f'Processed: {result}', flush=True)
        
    return n_ds



def processL10(params, verbose=True):
    """Function applying the actual baseline fit


    Parameters
    ----------
    params : list
            list of parameters used in process

    """
    
    line, inDir, outDir, dfile, offsets, debug = params[0], params[1], params[2], params[3], params[4], params[5]
    
    # define some processing data first (maybe relocat to function later?)

    #logger.info('loading: ', os.path.join(inDir,dfile), ' for line: ', line)
    spec, data, hdr, hdr1 = loadL08Data(os.path.join(inDir,dfile), verbose=False)
    
    umixers = np.unique(data['MIXER'])
    band = hdr['BAND']
    # insert the coordinate corrections
    # Note: the coordinate correction is not yet final
    # and will be (iteratively) improved    
    mxoffs = getMixerOffsets(band, umixers, verbose=verbose)
    
    for i, mix in enumerate(umixers):
        azoff = mxoffs['az'][i]
        eloff = mxoffs['el'][i]
        
        msel = np.argwhere((data['scan_type'] == 'OTF') & (data['MIXER']==mix)).flatten()
        if msel.size > 0:
            ras = data['RA'][msel]
            decs = data['DEC'][msel]
            utime = data['UNIXTIME'][msel]
            
            glat = hdr['GON_LAT']
            glon = hdr['GON_LON']
            galt = hdr['GON_ALT']
            
            cc = SkyCoord(ra=ras*u.deg, dec=decs*u.deg, frame='icrs')
            balloon = EarthLocation(lat=glat, lon=glon, height=galt*u.m)
            otime = Time(utime, format='unix')
            
            aa = AltAz(location=EarthLocation(lat=glat, lon=glon, height=galt*u.m), obstime=otime)
            altaz = cc.transform_to(aa)
            
            # apply beam offsets
            naz = altaz.az.deg*u.deg + azoff*u.deg
            nalt = altaz.alt.deg*u.deg + eloff*u.deg
            
            ncc = SkyCoord(AltAz(az=naz, alt=nalt, obstime=otime, location=balloon))
            nradec = ncc.transform_to('icrs')
            
            data['RA'][msel] = nradec.ra.deg
            data['DEC'][msel] = nradec.dec.deg
        
    # now we have to save the data in a FITS file
    
    # insert delta_ra and delta_dec in table for compatibility with WCS
    # insert WCS information in header
    # all will be relative to the first ra/dec pair
#    cc0 = SkyCoord(data['ra'][osel[0]]*u.deg, data['dec'][osel[0]]*u.deg, frame='icrs')
#    cc = SkyCoord(data['ra']*u.deg, data['dec']*u.deg, frame='icrs')    
#    dra, ddec = cc0.spherical_offsets_to(cc)
    # add the columns to the data
#    cols = data.columns
#    new_cols = fits.ColDefs([
#        fits.Column(name='dra', format='D', array=dra.deg),
#        fits.Column(name='ddec', format='D', array=ddec.deg)])
#    nhdu = fits.BinTableHDU.from_columns(cols + new_cols)
#    data = nhdu.data

    
    # change the spectral axis from IF frequency to velocity
    IF_freq = (np.arange(hdr['NPIX'])-hdr['CRPIX1'])*hdr['CDELT1']+hdr['CRVAL1']
    vlsr    = (hdr['IF0'] - IF_freq)/hdr['LINEFREQ']*const.c.value/1.e3 + hdr['VLSR'] # Vlsr in km/s
    hdr.set('CUNIT1', value='km/s', comment='Spectral unit: velocity')
    hdr.set('IFPIX0', value=hdr['CRPIX1'], comment='')
    hdr.set('IFDELT0', value=hdr['CDELT1'], comment='')
    hdr.set('IFVAL0', value=hdr['CRVAL1'], comment='')
    hdr.set('CRPIX1', value=0.000, comment=(''))
    hdr.set('CRVAL1', value=vlsr[0], comment=(''))
    hdr.set('CDELT1', value=np.diff(vlsr).mean(), comment=(''))

# add 
#    hdr.set('CDELT2', value=0.000000001, comment=(''), after='CDELT1')
#    hdr.set('CRVAL2', value=(data['ra'][osel[0]]*u.deg).value, comment=(''), after='CDELT1')
#    hdr.set('CRPIX2', value=0, comment=(''), after='CDELT1')
#    hdr.set('CUNIT2', value='deg', comment=(''), after='CDELT1')
#    hdr.set('CTYPE2', value='RA---GLS', comment=(''), after='CDELT1')
#    hdr.set('CDELT3', value=0.000000001, comment=(''), after='CDELT2')
#    hdr.set('CRVAL3', value=(data['dec'][osel[0]]*u.deg).value, comment=('CDELT2'))
#    hdr.set('CRPIX3', value=0, comment=(''), after='CDELT2')
#    hdr.set('CUNIT3', value='deg', comment=(''), after='CDELT2')
#    hdr.set('CTYPE3', value='DEC--GLS', comment=(''), after='CDELT2')
        
    hdr['DLEVEL'] = 1.0

    tred = Time(datetime.datetime.now()).fits
    hdr['PROCTIME'] = tred
    
    hdr.set('', value='', after='CALMETHD')
    hdr.set('', value='          Level 1.0 Pipeline Processing', after='CALMETHD')
    hdr.set('', value='', after='CALMETHD')
    hdr.set('L10PTIME', value=tred, comment=('L1.0 pipeline processing time'))
    
    # select only the processed good data
    osel = np.argwhere(data['scan_type'] == 'OTF').flatten()
    odata = data[osel]

    os.makedirs(outDir, exist_ok=True)
    ofile = os.path.join(outDir, os.path.split(dfile)[1].replace('_L09.fits','_L10.fits'))
    fits.writeto(ofile, data=None, header=hdr, overwrite=True)
    fits.append(ofile, data=odata, header=hdr1)

    return dfile


def getMixerOffsets(band, mixers, type='THEORY', offsetfile=None, verbose=False):
    """Function retrieving the GUSTO on-sky mixer offsets from file.

    usage:
    ------
    
    aa = getMixerOffsets(1, [3, 5, 8], verbose=True)
    print(aa['az'])     # prints: [0.06079  0.062584 0.093356]
    print(aa[0]['az'])  # prints: 0.06079
    print(aa['az'][0])  # prints: 0.06079
    
    Parameters
    ----------
    band : str
        GUSTO band: 'B1' or 'B2'
    mixers : int array
        array of mixer indices
    type : str
        type of offset determination: 'Theory' or 'AS_MEASURED'
    offsetfile : str
        path to file with mixer offsets
    verbose : bool
        flag for output to STDOUT
        
    Returns
    -------
    function returns rec array with mixer offsets
    """
    if offsetfile is None:
        offsetfile = offsetsfile0

    data = np.genfromtxt(offsetfile, delimiter='\t', skip_header=2, 
                         dtype=[('mxpix', 'U4'), ('az', '<f8'), ('el', '<f8'), ('type', 'U16')])
    
    cmixers = ['B%iM%i'%(band, i) for i in mixers]    
    ar = [np.argwhere((cmixer == data['mxpix'])&((type==data['type'])|(data['type']=='FIDUCIAL'))).flatten() for cmixer in cmixers]

    return data[ar].flatten()
