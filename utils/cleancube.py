import os
import sys
import glob
import numpy as np
import numpy.ma as ma
from astropy.io import fits
from astropy import units
from astropy.coordinates import SkyCoord, Galactic
from astropy.wcs import WCS
from astropy.wcs.wcsapi import SlicedLowLevelWCS
from astropy.utils import data
from spectral_cube import SpectralCube
import math


gustodir = '/data/scratch/GUSTO/gusto-datasystem/'
l2dir = f'{gustodir}/Data/level2/'


#target = 'GP'
#line = 'CII'
target = sys.argv[1]
line = sys.argv[2]

file = glob.glob(f'{l2dir}/{target}/{target}_{line}_6.fits')

hdu = fits.open(file[0])
cube0 = hdu[0].data
hdr0  = hdu[0].header
weight = hdu[1].data
hdrw   = hdu[1].header

#cube = np.zeros_like(cube0)
#w = np.zeros_like(weight)

qgood = np.where((weight > 5.0) & (weight < 1000) & ( cube0 > -10) & (cube0 < 20) & np.isfinite(cube0) & np.isfinite(weight))
cube = ma.masked_invalid(cube0)
w = ma.masked_invalid(weight)
ma.masked_greater(cube,20,False)
ma.masked_less(cube,-10,False)
cube = ma.masked_where((w < 3.0),cube,True) 
cube = ma.masked_where((w > 1000),cube,True) 
print(cube.shape)
vels=cube.shape[0]
lats=cube.shape[1]
lons=cube.shape[2]

#clean the data along legs at each velocity (simple median subtract)
for v in np.arange(vels):
    velslice = cube[v,:,:]
    bgd = ma.median(velslice,axis=0)
    cube[v,:,:] -= bgd

#outcube = np.zeros_like(cube0)

medw = ma.median(w)
print(cube.min(),cube.max(),w.max(),w.min(),medw)

hdulist = fits.PrimaryHDU(ma.filled(cube,-9999),hdr0)
outfile = f'{target}_{line}_clean.fits'
hdulist.writeto(outfile,overwrite=True)

cleaned = SpectralCube.read(outfile)
int_image = cleaned.moment0()
int_out = int_image.write(f'{target}_{line}_integrated.fits',overwrite = True)


