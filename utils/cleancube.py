import os
import sys
import glob
import numpy as np
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

cube = np.zeros_like(cube0)
w = np.zeros_like(weight)

qgood = np.where((weight > 5.0) & (weight < 1000) & ( cube0 > -4) & (cube0 < 10) & np.isfinite(cube0) & np.isfinite(weight))
cube[qgood] = cube0[qgood].copy()
w[qgood] = weight[qgood].copy()


medw = np.median(w[qgood])

print(cube.min(),cube.max(),w.max(),w.min(),medw)

hdulist = fits.PrimaryHDU(cube,hdr0)
outfile = f'{target}_{line}_clean.fits'
hdulist.writeto(outfile,overwrite=True)

cleaned = SpectralCube.read(outfile)
int_image = cleaned.moment0()
int_out = int_image.write(f'{target}_{line}_integrated.fits',overwrite = True)


