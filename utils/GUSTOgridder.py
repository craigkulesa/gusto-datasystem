#from Find_stw import find_fringes

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic
from astropy import units as u
from astropy import wcs
import numpy as np
import scipy as sp
import glob
import math
from astropy import constants as const
from grid_otf import grid_otf
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

def make_gusto_array(directory, linename, vel_vector, coordType):
    # level 1 calibrated spectra in directory 
    # Line to make cube of in line_str (either NII or CII)
    # velocity vector to interpolate Leve 1 data onto
    input_files = glob.glob(f'{directory}/{linename}*fits')
    #ifile = os.path.join(directory,input_filename)
    #print(ifile)        
    nfile= len(input_files)
    print(f'nchan: {vel_vector.size}')
    arr_line=np.array([],[])
    xpos = np.array([])
    ypos = np.array([])
    legweight = np.array([])

    print(f'{nfile} files found for cube generation in {directory}.')

    for ifile in input_files:
        
        spec, data, hdr, hdr1 = loadSDFITS(ifile, verbose=False)
        rowFlag = data['ROW_FLAG']
        n_spec, n_pix = spec.shape
        chanflag = data['CHANNEL_FLAG']
        
        
        # new level 1 files are already in V_lsr, no longer need this code
        ## compute velocity                                                           
        npix    = hdr['NPIX']            
        #IF_pix  = hdr['CRPIX1']
        #IF_val  = hdr['CRVAL1']
        #IF_del  = hdr['CDELT1']
        #IF_freq = (np.arange(npix)-IF_pix)*IF_del+IF_val
        #VLSR    = hdr['VLSR']        
        #
        #IF_vlsr0= hdr['IF0']
        #vlsr    = (IF_vlsr0 - IF_freq)/line_freq*constants.c.value/1.e3 + VLSR # Vlsr in km/s 
        
        # new level 1 V_lsr import
        VLSR_pix  = hdr['CRPIX1']
        VLSR_val  = hdr['CRVAL1']
        VLSR_del  = hdr['CDELT1']
        vlsr = (np.arange(npix)-(VLSR_pix))*VLSR_del+VLSR_val
        #print(vlsr.min(),vlsr.max())
        qsort = vlsr.argsort()
        vlsr.sort()
    
        #print(np.all(vlsr[1:] > vlsr[:-1]))



        
        line_freq = hdr['LINEFREQ']
        
        if (linename == "CII") & (line_freq < 1900500):
            continue
        if (linename == "NII") & (line_freq > 1900500):
            continue
        
        rfl = RowFlags.MIXER_UNPUMPED | RowFlags.MIXER_MISPUMPED
        rfl = RowFlags.RINGING_BIT0 | RowFlags.RINGING_BIT1

        #osel = np.argwhere((data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0)).flatten()
        #osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & 0x60)==0) & ((data['MIXER']==5) | (data['MIXER']==8)) & (data['rms']<5)).flatten()
        if (linename == "CII"):
            #osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & 0x60)==0) & ((data['MIXER']==8))).flatten()
            osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & rfl)==0) & ((data['MIXER']==5) | (data['MIXER']==8))).flatten()
        if (linename == "NII"):
            osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & rfl)==0) & ((data['MIXER']==2) | (data['MIXER']==3) | (data['MIXER']==6))).flatten()
        if len(osel) <= 0:
            print('WARNING: No OTF spectra available in ',input_filename)
            # logger.warning('No OTF spectra available.')                           
        else:
            spec_OTF = np.squeeze(spec[osel,:])
            chan_OTF = np.squeeze(chanflag[osel,:])
            data_OTF = np.squeeze(data[osel])
            mxrms    = data_OTF['rms']
            n_OTF, n_otfpix = spec_OTF.shape
            #x=np.arange(n_otfpix)
    
            # Instantiate for ra,dec->l,b transform                                      
            c_ra_dec = SkyCoord(ra=data_OTF['RA']*u.degree, dec=data_OTF['DEC']*u.degree, frame='icrs')
    
            #basecorr = np.zeros(spec_OTF.shape)
            #rmsotf = np.zeros(n_OTF)
            #rf = np.zeros(n_OTF)
            
            # empty arrays to fill                                                        
            
            

            if coordType[0][0] == 'G':
                c_l_b = c_ra_dec.transform_to(Galactic)    # transform to l,b
                leg_y = c_l_b.b
                leg_x = c_l_b.l
            else:
                leg_y=c_ra_dec.dec
                leg_x=c_ra_dec.ra
            #print(len(leg_l),len(leg_b))
            leg_spec=[]
            leg_chf=[]
            wgt1 = np.ones(n_otfpix)
            for i0,spec_OTF1 in enumerate(spec_OTF):
                chan_OTF1 = chan_OTF[i0]
                chan1 = chan_OTF1[qsort]
                qchan = np.argwhere(chan1 == 0)
                spec1 = spec_OTF1[qsort]
                wgt1[qchan] = 0.00
                arr_line1 = np.interp(vel_vector,vlsr, spec1)
                arr_wgt1 = np.interp(vel_vector,vlsr, wgt1)
                leg_spec.append(arr_line1)
                leg_chf.append(arr_wgt1)
            leg_spec=np.array(leg_spec)
            leg_chf =np.array(leg_chf)

            #stack all usable spectra into new array
            if arr_line.shape[0] == 0:
                arr_line = leg_spec
                arr_chf = leg_chf
            else:
                arr_line = np.vstack((arr_line,leg_spec))
                arr_chf = np.vstack((arr_chf,leg_chf))
                
            #
            xpos = np.append(xpos,leg_x.degree)
            ypos = np.append(ypos,leg_y.degree)
            legweight = np.append(legweight,1.0/mxrms**2)

            
    arr_line = np.array(arr_line)
    xpos = np.array(xpos)
    ypos = np.array(ypos)
    nchan = vel_vector.shape[0]
    print(arr_line.shape,legweight.shape,xpos.shape,ypos.shape)
    return arr_line, xpos, ypos, legweight, nchan, line_freq, arr_chf
            


	
	
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
        my_parser.add_argument('-l', 
                               metavar='--VLSRrange',
                               nargs=2,
                               required=False,
                               help='minimum maximum velocity channel, default -200 200 km/s',
                               default=[-200, 200])
        args = my_parser.parse_args()

    source = args.s
    line_str = args.b
    kern = args.k
    vinput = args.l[0]
    print(float(vinput))
    vmin = float(args.l[0])
    vmax = float(args.l[1])
    

    print(args)
    
    #dir = '/Users/umit/Desktop/STO2_etacar5_data-redution/Pipeline_HOTneeded/Gum31_4591-4733/'
    datadir = '/data/scratch/GUSTO/gusto-datasystem/Data/'
    dir_level1 = f'{datadir}/level1/{source}/'
    dir_write = f'{datadir}/level2/{source}/'
    
    vel_spacing = 2.0 # km/s
    
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
    arr_line0, xpos0, ypos0, weight, nchan0, restfreq, arr_chf = make_gusto_array(dir_level1,line_str,vv_in,coordType)
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
    
    #xpos_in = np.append(xpos0[0,:],xpos1[1,:])
    #ypos_in = np.append(ypos0[0,:],ypos1[1,:])
    xpos_in = xpos0
    ypos_in = ypos0
    beam_fwhm_in = beam_fwhm
    #print(beam_fwhm_in)
    pix_scale = int(3600.0*beam_fwhm_in/pixPerBeam)/3600.0
    #
    print(f'Beam: {beam_fwhm*60:.2} arcmin.  Pixel scale: {pix_scale*60:.2}')
    # create header, wcs, and image size from header and given parameters
    hdr, wcsObj, xsize, ysize = create_wcsheader(xpos_in,ypos_in,restfreq,vv_in,coordType,pix_scale,beam_fwhm_in)
    #
    # create spectral map 
    cube, weight, beam_size = grid_otf(arr_line_in, xpos_in, ypos_in, wcsObj, nchan_in, xsize, ysize, pix_scale, beam_fwhm_in, weight=weight ,kern = kern)
    #cube, weight, beam_size = grid_otf(arr_line_in, xpos_in, ypos_in, wcsObj, nchan_in, xsize, ysize, pix_scale, beam_fwhm_in, kern = kern)
    #
    #
    hdr['CTYPE3'] = 'VRAD'
    hdr['CUNIT3'] = 'm/s'
    hdr['CRVAL3'] = vv_in[0]
    hdr['CRPIX3'] = 1.0
    hdr['CDELT3'] = (vv_in[1]-vv_in[0])
    #
    silentremove(dir_write+f'cube_{line_str}.fits')
    hdu_cube_out = fits.PrimaryHDU(cube.data, header = hdr)
    hdu_cube_out.writeto(dir_write+f'{source}_{line_str}_{pixPerBeam}pix_in_{beam_fwhm*60:0.2}_{KT}.fits',overwrite=True)
    #
    #
    #
    #
if __name__ == "__main__":
    main()
