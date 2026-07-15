#from Find_stw import find_fringes

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic
from astropy import units as u
from astropy import wcs
import numpy as np
import numpy.ma as ma
import scipy as sp
import glob
import math
from astropy import constants as const
from grid_otf_optimized import grid_otf
import sys
import warnings
#from progressbar import ProgressBar

import math
import numpy
import sys
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
import scipy
import time
import os,errno
#from gustoL09P.GL09PDataIO import loadL08Data
from GUSTO_Pipeline.DataIO import loadSDFITS 
from multiprocessing import Process, Queue
from flagdefs import *

import datetime
import argparse
import configargparse

def silentremove(filename):
# remove files without raising error
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred


def call_fits(dir):
	file = glob.glob(dir+'*.fits')
	nfile = len(file)
	hdu_tot = []
	
	print('Reading calibrated spectra', flush = True)
    #bar0 = ProgressBar()
	for i in range(0,nfile):
		hdu = fits.open(file[i])
		hdu_tot.append(hdu)

	return hdu_tot

def make_gusto_array(directory, linename, mx, vel_vector, coordType):
    # level 1 calibrated spectra in directory
    # Line to make cube of in line_str (either NII or CII)
    # velocity vector to interpolate Level 1 data onto
    input_files = glob.glob(f'{directory}/{linename}*fits')
    nfile = len(input_files)
    print(f'nchan: {vel_vector.size}')
    print(f'{nfile} files found for cube generation in {directory}.')

    # Accumulate per-file results in lists; single vstack/concatenate at end
    # avoids O(N^2) copying from repeated ma.vstack / np.append in a loop.
    all_specs = []
    all_chf   = []
    all_x     = []
    all_y     = []
    all_wgt   = []

    for ifile in input_files:

        spec, data, hdr, hdr1 = loadSDFITS(ifile, verbose=False)
        rowFlag  = data['ROW_FLAG']
        n_spec, n_pix = spec.shape
        chanflag = data['CHANNEL_FLAG']

        # velocity axis — already in V_lsr
        npix     = hdr['NPIX']
        VLSR_pix = hdr['CRPIX1']
        VLSR_val = hdr['CRVAL1']
        VLSR_del = hdr['CDELT1']
        vlsr = (np.arange(npix) - VLSR_pix) * VLSR_del + VLSR_val
        qsort = vlsr.argsort()
        vlsr  = vlsr[qsort]          # sorted in-place equivalent, keeps qsort valid

        line_freq = hdr['LINEFREQ']

        if (linename == "CII") & (line_freq < 1900500):
            continue
        if (linename == "NII") & (line_freq > 1900500):
            continue

        rfl = RowFlags.MIXER_UNPUMPED | RowFlags.MIXER_MISPUMPED
        rfl = RowFlags.RINGING_BIT0  | RowFlags.RINGING_BIT1

        if (linename == "CII"):
            if mx == 2:
                osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & rfl)==0) & (data['MIXER']==2 )).flatten()
            elif mx == 5:
                osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & rfl)==0) & (data['MIXER']==5 )).flatten()
            elif mx == 8:
                osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & rfl)==0) & (data['MIXER']==8 )).flatten()
            else:
                osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & rfl)==0) & ((data['MIXER']==5) | (data['MIXER']==8))).flatten()
        if (linename == "NII"):
            if mx == 2:
                osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & rfl)==0) & (data['MIXER']==2) ).flatten()
            elif mx == 3:
                osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & rfl)==0) & (data['MIXER']==3) ).flatten()
            elif mx == 6:
                osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & rfl)==0) & (data['MIXER']==6) ).flatten()
            else:
                osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & rfl)==0) & ((data['MIXER']==2) | (data['MIXER']==3) | (data['MIXER']==6))).flatten()

        if len(osel) <= 0:
            print('WARNING: No OTF spectra available in ', ifile)
            continue

        spec_OTF = np.squeeze(spec[osel, :])
        chan_OTF = np.squeeze(chanflag[osel, :])
        data_OTF = np.squeeze(data[osel])
        mxrms    = data_OTF['rms']
        n_OTF, n_otfpix = spec_OTF.shape

        # Coordinate transform
        c_ra_dec = SkyCoord(ra=data_OTF['RA']*u.degree, dec=data_OTF['DEC']*u.degree, frame='icrs')
        if coordType[0][0] == 'G':
            c_l_b = c_ra_dec.transform_to(Galactic)
            leg_y = c_l_b.b
            leg_x = c_l_b.l.wrap_at(180*u.deg)
        else:
            leg_y = c_ra_dec.dec
            leg_x = c_ra_dec.ra

        # ------------------------------------------------------------------
        # Vectorized interpolation — replaces the per-spectrum Python loop.
        #
        # Pre-compute interpolation indices and weights once for vel_vector
        # vs vlsr, then apply to all spectra simultaneously.
        #
        # Also fixes a bug in the original: wgt1 was never reset between
        # spectra, so bad-channel flags from spectrum N leaked into N+1.
        # ------------------------------------------------------------------

        # Sort all spectra and channel flags along the velocity axis at once
        spec_sorted = spec_OTF[:, qsort]   # (n_OTF, n_otfpix)
        chan_sorted  = chan_OTF[:, qsort]   # (n_OTF, n_otfpix)

        # Interpolation index and fractional weight (same for every spectrum)
        idx = np.searchsorted(vlsr, vel_vector)
        idx = np.clip(idx, 1, len(vlsr) - 1)
        lo  = idx - 1
        hi  = idx
        dv  = vlsr[hi] - vlsr[lo]
        # Guard against zero-width intervals (shouldn't occur with real data)
        dv  = np.where(dv == 0, 1.0, dv)
        t   = (vel_vector - vlsr[lo]) / dv  # (nchan_out,)

        # Interpolate all spectra at once: (n_OTF, nchan_out)
        arr_lines = spec_sorted[:, lo] + t * (spec_sorted[:, hi] - spec_sorted[:, lo])
        arr_chans = chan_sorted[:, lo]  + t * (chan_sorted[:, hi] - chan_sorted[:, lo])

        # Build masked array: mask where channel flag > 0
        leg_spec = ma.MaskedArray(arr_lines, mask=(arr_chans > 0))
        leg_chf  = arr_chans

        all_specs.append(leg_spec)
        all_chf.append(leg_chf)
        all_x.append(leg_x.degree)
        all_y.append(leg_y.degree)
        all_wgt.append(1.0 / mxrms**2)

    # Single stack at the end — avoids O(N^2) copies from repeated vstack
    arr_line  = ma.vstack(all_specs)
    arr_chf   = np.vstack(all_chf)
    xpos      = np.concatenate(all_x)
    ypos      = np.concatenate(all_y)
    legweight = np.concatenate(all_wgt)

    print(np.min(xpos), np.max(xpos), np.min(ypos), np.max(ypos), np.median(xpos), np.median(ypos))

    # Filter out positions more than 1.5 degrees from median
    xmed  = np.median(xpos)
    ymed  = np.median(ypos)
    dliml  = 21.5
    dlimb = 2.0
    qkeep = np.argwhere((np.abs(xpos - xmed) < dliml) & (np.abs(ypos - ymed) < dlimb))
    arr_linekeep = np.squeeze(arr_line[qkeep, :])
    xkeep        = np.squeeze(xpos[qkeep])
    ykeep        = np.squeeze(ypos[qkeep])
    wgtkeep      = np.squeeze(legweight[qkeep])
    arr_chfkeep  = np.squeeze(arr_chf[qkeep, :])
    nchan        = vel_vector.shape[0]
    print(arr_line.shape, legweight.shape, xpos.shape, ypos.shape)
    return arr_linekeep, xkeep, ykeep, wgtkeep, nchan, line_freq, arr_chfkeep
            


	
	
def get_restfreq(hdu):
	data_hdu = hdu[1].data
	restfreq = data_hdu.field('RESTFREQ')
	return restfreq

def get_vel_freq(hdu):
	header = hdu[1].header
	data_hdu = hdu[1].data
	n_pixl = data_hdu.field('MAXIS1')[0]
	restfreq = data_hdu.field('RESTFREQ')
	n_line = len(restfreq)
	vv = np.zeros(n_pixl*n_line).reshape(n_line, n_pixl)
	freq = np.zeros(n_pixl*n_line).reshape(n_line, n_pixl)
	for j0 in range(0,n_line):
		vv[j0,:] = (np.float(hdu[1].header['CRVAL1']) + (1 + np.arange(n_pixl) - data_hdu.field('CRPIX1')[j0]) * data_hdu.field('CDELT1')[j0]) 
		vv[j0,:] = vv[j0,:]*1.e2
		freq[j0,:] = restfreq[j0]* (1.- vv[j0,:]/const.c.cgs.value)
	return vv, freq
		
def make_header(xref, yref, xsize, ysize, pix_scale, xref_pix, yref_pix, coordType, radesys, equinox, frest, faxis, beam_fwhm, veldef, specsys, proj="SFL"):

    hdr = fits.Header()

    # BASIC stuff, the WCS code needs this
    hdr['SIMPLE'] = True
    #hdr['NAXIS'] = 4
    hdr['NAXIS1'] = xsize
    hdr['NAXIS2'] = ysize
    hdr['NAXIS3'] = len(faxis)
    #hdr['NAXIS4'] = 1

    ctypeDashes = '----'

    xctype = coordType[0] + ctypeDashes[len(coordType[0]):]
    yctype = coordType[1] + ctypeDashes[len(coordType[1]):]

    # MAKE THE POSITION AXES
    hdr['CTYPE1'] = xctype + '-' + proj
    hdr['CRVAL1'] = xref
    hdr['CRPIX1'] = xref_pix
    hdr['CDELT1'] = -1.0*pix_scale

    hdr['CTYPE2'] = yctype + '-' + proj
    hdr['CRVAL2'] = yref
    hdr['CRPIX2'] = yref_pix
    hdr['CDELT2'] = pix_scale

    # MAKE THE VELOCITY AXIS (ALONG THE THIRD DIMENSION)
    # the frame is now indicated via SPECSYS.  Check on any other
    # needed WCS keywords for use here.
    hdr['CTYPE3'] = 'VEL-LSR'
    hdr['CUNIT3'] = 'km/s'
    hdr['CRVAL3'] = faxis[0]
    hdr['CRPIX3'] = 1.0
    hdr['CDELT3'] = faxis[1]-faxis[0]

    # STOKES axis - always I
    #hdr['CTYPE4'] = 'STOKES'
    #hdr['CRVAL4'] = 1.0
    #hdr['CRPIX4'] = 1.0
    #hdr['CDELT4'] = 1.0

    hdr['SPECSYS'] = specsys

    # AIPS velocity type
    hdr['VELREF'] = 0
    if specsys == "LSRK":
        hdr['VELREF'] = 1
    elif specsys == "HELIOCEN":
        hdr['VELREF'] = 2
    elif specsys == "TOPOCENT":
        hdr['VELREF'] = 3
    # no others are defined in the original AIPS memo, should search for updates
    # for now, leave everything else at 0
    if veldef == "RADI":
        # radio definition adds 256
        hdr['VELREF'] = hdr['VELREF'] + 256
    # AIPS memo doesn't say what to do for relativistic velocity definition

    # Set the ALT* axis keywords if possible
    #if hdr['CDELT3'] != 0. and frest > 0.:
    #    # zero velocity
    #    hdr['ALTRVAL'] = 0.0
    #    # is at channel here the frequency axis equals the rest frequency
    #    hdr['ALTRPIX'] = hdr['CRPIX3'] + (frest - hdr['CRVAL3'])/hdr['CDELT3']
        
    hdr['RESTFRQ'] = frest

    # ADD THE RADESYS and EQUINOX when appropriate
    if radesys is not None and len(radesys) > 0:
        hdr['RADESYS'] = radesys
    if equinox is not None and equinox > 0.0:
        hdr['EQUINOX'] = equinox

    return hdr
	
def create_wcsheader(xpos,ypos,restfreq,vv_in,coordType,pix_scale,beam_fwhm):
    # image size
    xRange = np.max(xpos)-np.min(xpos)
    yRange = np.max(ypos)-np.min(ypos)
    xsize = int(math.ceil(xRange*1.1/pix_scale))+20
    ysize = int(math.ceil(yRange*1.1/pix_scale))+20
    # set image center
    refXsky = np.min(xpos)+0.5*xRange
    refYsky = np.min(ypos)+0.5*yRange
    refXpix = math.ceil(xsize*0.5)
    refYpix = math.ceil(ysize*0.5)
    #xcoord = 'GLON'
    #ycoord = 'GLAT'
    specSysDict = {'OBS':'TOPOCENT','GEO':'GEOCENTR','BAR':'BARYCENT','HEL':'HELIOCEN','GAL':'GALACTOC','LSD':'LSRD','LSR':'LSRK','LGR':'LOCALGRP','COB':'CMBDIPOL'}
    #coordType = [xcoord,ycoord]
    
    radesys = ''
    equinox = 0.
    veldef = 'RADI'
    #specsys = specSysDict[header['VELFRAME']]
    specsys = 'LSR'  #or LSRK
    # create header for the spectral cube
    hdr = make_header(refXsky, refYsky, xsize, ysize, pix_scale, refXpix, refYpix, coordType, radesys, equinox, restfreq, vv_in, beam_fwhm, veldef, specsys)
    # create wcs object from STO2 header (non-trivial header)
    w = wcs.WCS(hdr,relax=True)
    return hdr, w, xsize, ysize

#
# begine of main program
#	
def main(args=None,verbose=True):
    if args==None:
        # Create the input parser
        my_parser = argparse.ArgumentParser(prog='GUSTOgridder',
                                            usage='%(prog)s source band',
                                            description='Regrid level1 data in source directory of specified band')
        my_parser.version = "Version 0.0.2 (22 Jan 2026) "
        my_parser.add_argument('-v', action='version')

        my_parser.add_argument('-s',
                               metavar='--source',
                               required=True,
                               type=str,
                               help = 'Name of source directory in level1. galactic coords maps for source G???, RADEC otherwise')
        my_parser.add_argument('-b',
                               metavar='--band',
                               type=str,
                               required=True,
                               help='NII or CII')
        my_parser.add_argument('-k', 
                               metavar='--kernel',
                               type=str,
                               required=False,
                               help='gridding kernel: gaussbessel, gauss, or nearest.  Default is gaussbessel ',
                               default='gaussbessel')
        my_parser.add_argument('-x', 
                               metavar='--mixer',
                               required=False,
                               help='NII mixer: 2, 3, 6 or CII mixer: 5, 8 or 0 for all mixers in band',
                               default=0)
        my_parser.add_argument('-o', 
                               metavar='--ofile',
                               required=False,
                               help='Output cube name.  If none provided a name based on target and regridding parameters is created')
        my_parser.add_argument('-P', 
                               metavar='--pixBeam',
                               required=False,
                               help='pixels per beam, default is 3',
                               default='3')
        my_parser.add_argument('-Beam', 
                               metavar='--BeamFWHM',
                               required=False,
                               help='beam FWHM in decimal arcmin, default use data header',
                               default='header')
        my_parser.add_argument('-dv', 
                               metavar='--vel_spacing',
                               required=False,
                               help='Velocity spacing in output cube in km/s',
                               default='header')
        my_parser.add_argument('-l', 
                               metavar='--VLSRrange',
                               nargs=2,
                               required=False,
                               help='minimum maximum velocity channel, default -200 200 km/s',
                               default=[-200, 200])
        my_parser.add_argument('-wf', 
                               metavar='--wcsfile',
                               required=False,
                               help='Input fits cube to match WCS if not present the WCS will be made based on input L1 scans')


        args = my_parser.parse_args()

    source = args.s
    line_str = args.b
    kern = args.k
    ofile=args.o
    vinput = args.l[0]
    mx = int(args.x)
    print(float(vinput))
    vmin = float(args.l[0])
    vmax = float(args.l[1])
    wcsfile = args.wf
    

    print(args)
    
    #dir = '/Users/umit/Desktop/STO2_etacar5_data-redution/Pipeline_HOTneeded/Gum31_4591-4733/'
    datadir = '/data/scratch/GUSTO/gusto-datasystem/Data/'
    dir_level1 = f'{datadir}/level1/{source}/'
    dir_write = f'{datadir}/level2/{source}/'
    
    dvNII = 2.0076146439883598  # band 1 native resolution
    dvCII = 0.7709722465531635  # band 2 native resolution
    # velocity steps
    if args.dv == 'header':
        if args.b == 'NII':
            vel_spacing = dvNII
        else:
            vel_spacing = dvCII
    else:
        vel_spacing = float(args.dv)

    vv_in = np.arange(vmin,vmax,vel_spacing)
    ktypes = ['B', 'G', 'N']
    match kern:
        case 'gaussbessel':
            KT = ktypes[0]
        case 'gauss': 
            KT = ktypes[1]
        case 'nearest':
            KT = ktypes[2]
        case _:
            print('Not a valid kernel: use gaussbessel (default), gauss or nearest')
            return

    if source[0] == 'G':
        xcoord='GLON'
        ycoord='GLAT'
        print('Data cube in Galactic Cooordinates')
    else:
        xcoord='RA'
        ycoord='DEC'
        print('Data cube in RA/DEC Cooordinates')
    coordType = [xcoord,ycoord]
    # read all calibrated fits data, at all positions
    print(f'Input dir {dir_level1} Line {line_str} Velocity array {vv_in.shape}')
    arr_line0, xpos0, ypos0, weight, nchan0, restfreq, arr_chf = make_gusto_array(dir_level1,line_str,mx,vv_in,coordType)
    #os.system('ls')
            
    restfreq *= 1e6 # convert to Hz
    # dish size of STO2 in cm
    dish_diam = 90.
    # wavelength of lines in cm
    wavelength = const.c.cgs.value/restfreq
    #
    # beam size in array 
    if args.Beam == 'header':
        beam_fwhm = 1.2 * wavelength/dish_diam * np.rad2deg(1.)
    else:
        beam_fwhm = float(args.Beam)/60.0

    # create spectra array to put in regridder
    #arr_line0, xpos0, ypos0, nchan0 = make_line_array(hdu)
    #
    # pixel size
    pixPerBeam = float(args.P) 
    #
    # mask nan channels and arrange variables for header and re-gridding
    
    
    
    
    
    nchan_in = len(vv_in)
    # 
    arr_line_in = arr_line0
    print('Input array is masked',ma.is_masked(arr_line_in))
    
    #xpos_in = np.append(xpos0[0,:],xpos1[1,:])
    #ypos_in = np.append(ypos0[0,:],ypos1[1,:])
    xpos_in = xpos0
    ypos_in = ypos0
    beam_fwhm_in = beam_fwhm
    #print(beam_fwhm_in)
    pix_scale = int(3600.0*beam_fwhm_in/pixPerBeam)/3600.0
    #
    print(f'Beam: {beam_fwhm*60:.2} arcmin.  Pixel scale: {pix_scale*60:.2}')
    # Use existing WCS or create header, wcs, and image size from header and given parameters
    if wcsfile != None:
        hdu_wcs = fits.open(wcsfile)
        hdr = hdu_wcs[0].header
        wcsObj = wcs.WCS(hdr,relax = True)
        xsize = hdr['NAXIS1']
        ysize = hdr['NAXIS2']
    else:
        hdr, wcsObj, xsize, ysize = create_wcsheader(xpos_in,ypos_in,restfreq,vv_in,coordType,pix_scale,beam_fwhm_in)
    #
    # create spectral map 
    cube, weight, beam_size = grid_otf(arr_line_in, xpos_in, ypos_in, wcsObj, nchan_in, xsize, ysize, pix_scale, beam_fwhm_in, weight=weight ,kern = kern)
    #cube, weight, beam_size = grid_otf(arr_line_in, xpos_in, ypos_in, wcsObj, nchan_in, xsize, ysize, pix_scale, beam_fwhm_in, kern = kern)
    #
    qzero = weight == 0.0
    bzero = -9999
    cube[qzero] = bzero
    
    #
    hdr['CTYPE3'] = 'VRAD'
    hdr['CUNIT3'] = 'km/s'
    hdr['CRVAL3'] = vv_in[0]
    hdr['CRPIX3'] = 1.0
    hdr['CDELT3'] = (vv_in[1]-vv_in[0])
    hdr['BZERO']  = bzero
    hdr['BSCALE'] = 1.0
    hdr['OBJECT'] = source 
    hdr['LINE'] = (line_str, 'Observed line') 
    hdr['KERNEL'] = (KT, 'Regridder, B:GAUSSBESSEL, G:GAUSS, N:NEAREST ') 
    hdr['BEAMFWHM'] = (beam_fwhm*60,'beam in arcminutes')
    hdr['PIX_BEAM'] = (pixPerBeam,'Pixels per beam')
    hdr['MIXER'] = (mx,'mixer ID, 0 for all usable mixers') 
    hdr['MINVEL'] = args.l[0]
    hdr['MAXVEL'] = args.l[1]
    #wcsfile = args.f
    

    #
    silentremove(dir_write+f'cube_{line_str}.fits')
    hdu_cube_out = fits.PrimaryHDU(cube, header=hdr)

    hdulist = fits.HDUList(hdu_cube_out)

    hduw = fits.ImageHDU(data = weight,name = 'WEIGHT')
    
    hdulist.append(hduw)
    if ofile == None:
        outcube = dir_write+f'{source}_{line_str}_{mx}_at_{beam_fwhm*60:0.2}_{KT}.fits'
    else:
        outcube = dir_write + ofile

    #hdu_cube_out.writeto(outcube ,overwrite = True)
    hdulist.writeto(outcube, overwrite = True)
    #
    #
    #
    #
if __name__ == "__main__":
    main()
