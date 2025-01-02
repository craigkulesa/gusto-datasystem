#
#
import os
import sys
import glob
import math
import random
import numpy as np
import argparse
import configparser
from pathlib import Path

from scipy import interpolate
from scipy.interpolate.interpnd import _ndim_coords_from_arrays
from scipy.spatial import cKDTree

from astropy import units as u
from astropy import constants
from astropy import constants

from astropy.io import fits
from astropy.io.fits import Header
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic

from pybaselines import Baseline
import matplotlib.pyplot as plt

# 10/29/2024 Makes L2 maps
# Little script to look at Volker's L0.9 data 
# NOT FULLY FEATURED YET, porting over from May 2024 codebase
# Does a disciplined least squares baseline fit to flatten L1 data
# Makes an integrated intensity map.
# ra,dec and l,b coordinates mish-mashed, working it...

# 11/19/24 make an animated gif of calibrated spectra
# Output TARGET-%05d.png files
# make palette:   ffmpeg -pattern_type glob -i 'NGC6334-*.png' -vf palettegen palette.png
# combine files:  ffmpeg -r 3 -pattern_type glob -i 'NGC6334-*.png' -i palette.png -lavfi paletteuse test.gif
# view:           mpv test.gif

# 12/12/2024
# TODO: Use as-is with the 4 Dec level09 commit CLASS format, or accept new 0.9/0.95 FITS format
# 4 Dec FITS format             new format
# =================             ==========
# single file per scanID        singlle file per scanID
# single mixer per file         all mixers from band in a single file (2-4)
# s-r/r with hots only          s-r/r with optional at 0.95 of baseline corr
# primary header zero           primary header is preserved from level 0.8
#   -- deltas primarily impact how velocity scale information is computed


################################################################################
def doStuff(self, args):
   my_file = Path(self)
   try:
       my_file.resolve(strict=True)
   except FileNotFoundError:
       print("file not found ", end="")
       return None
   else:
       hdu    = fits.open(self)

   ##
   ## Get Header Data
   ## 
   header = hdu[1].header
   data   = hdu[1].data
   spec   = data['DATA']
   ROW_FLAG  = data['ROW_FLAG']

   # compute velocity
   npix    = len(spec[0]) # 'NPIX' doesn't exist in 0.9 header
   IF_pix  = header['CRPIX1']
   IF_val  = header['CRVAL1']
   IF_del  = 4.887586
   IF_freq = (np.arange(npix)-IF_pix)*IF_del+IF_val
   VLSR    = header['VELO-LSR'] # 'VLSR' doesn't exist in 0.9 header

   IF_vlsr0= header['IF0']
   line_freq = header['LINEFREQ']
   vlsr    = (IF_vlsr0 - IF_freq)/line_freq*constants.c.value/1.e3 + VLSR # Vlsr in km/s

   # Onyl use some rows
   #mixer  = data['mixer']
   #row_mask = np.logical_and(mixer == 5, ROW_FLAG & mask == 0)

   n_OTF  = len(data)
   x = np.arange(1024)


   # Reference Coordinate System
   crval2 = header['CRVAL2']   # reference RA pixel
   crval3 = header['CRVAL3']   # reference DEC pixel

   # Delta Coordinate
   cdelt2 = hdu[1].data['CDELT2']   # array of RA offsets
   cdelt3 = hdu[1].data['CDELT3']   # array of DEC offsets


   # Instantiate for ra,dec->l,b transform
   c_ra_dec = SkyCoord(ra=(crval2+cdelt2)*u.degree, dec=(crval3+cdelt3)*u.degree, frame='icrs')

   # Setup baseline fitting and integrated intensity
   # baseline from -80 km/s to +16 km/s
   xlow   = np.argmin(np.abs(vlsr - 16))
   xhigh  = np.argmin(np.abs(vlsr - -80))
   # integrated intensity from -15.5 km/s to +4.5 km/s
   iilow  = np.argmin(np.abs(vlsr - 4.5)) - xlow
   iihigh = np.argmin(np.abs(vlsr - -15.5)) - xlow

   base = np.zeros([n_OTF,xhigh-xlow])
   baseline_fitter = Baseline(x_data=x[xlow:xhigh])
   y_flat = np.zeros(xhigh-xlow)

   # empty arrays to fill
   ii=np.array([])
   l=np.array([])
   b=np.array([])

   if args.plot:
       plt.clf()
   # Main routine run on every Level 1 fits file
   for i in range(n_OTF):
       mask = (1<<26)|(1<<27)   # Mask off ringing rows
       if ((int(ROW_FLAG[i]) & mask) == 0):
          c_l_b = c_ra_dec[i].transform_to(Galactic)    # transform to l,b

          # Baseline fitting
          spec_to_draw = spec[i,xlow:xhigh] - np.median(spec[i,xlow:xhigh])
          spec_new = spec_to_draw - base[i,0:(xhigh-xlow)]
          base2 = baseline_fitter.aspls(spec_new, 1e5)
          y_flat = spec_new - base2[0]

          # Fill integrated intensity and l,b
          ii = np.append(ii, sum(y_flat[iilow:iihigh]))
          l  = np.append(l, c_l_b.l.value[0])
          b  = np.append(b, c_l_b.b.value[0])
          if args.plot:
              plt.plot(vlsr[xlow:xhigh],y_flat)
       else:
          pass

   if args.plot:
       plt.ylim(-10, 80)
       plt.vlines(4.5, -5, 40)
       plt.vlines(-15.5, -5, 40)

       plt.savefig(f'NGC6334-{self[-16:-11]}.png') 

   #plt.show(block=False)
   #plt.pause(1)

   print(l.mean())

   # Return the (l,b) position and integrated intensity
   data = (l, b, ii)

   return data


def regrid(l, b, T, beam):
   # Interpolate data onto rectangular grid
   THRESHOLD = 0.05 # Any gridded data point Threshold (deg) distance from real data point is NaN

   # Calculate the range of glon and glat values
   l_min, l_max = np.min(l), np.max(l)
   b_min, b_max = np.min(b), np.max(b)

   # Calculate number of grid points
   N_l = int(np.ceil((l_max - l_min) / beam))
   N_b = int(np.ceil((b_max - b_min) / beam))

   # Create meshgrid
   l_grid, b_grid = np.meshgrid(np.linspace(l_min, l_max, N_l),np.linspace(b_min, b_max, N_b))

   # Initialize array
   image = interpolate.griddata((l, b), T, (l_grid, b_grid), method='linear')

   # KD-TREE and RESHAPE are used to NaN pixels greated than THRESHOLD distance from nearest real data points
   # Construct kd-tree, functionality copied from scipy.interpolate
   tree = cKDTree(np.c_[l, b])
   grid_points = np.c_[l_grid.ravel(), b_grid.ravel()]
   distances, _ = tree.query(grid_points, k=1)

   # Reshape distances back to grid shape
   distances = distances.reshape(l_grid.shape)

   # Copy original result and mask missing values with NaNs
   image[distances > THRESHOLD] = np.nan

   return image


def get_filenames(directory, scan_file):
    # Read filenames from text file generated by findSCANS.py
    with open(scan_file, 'r') as file:
        filenames = [line.strip() for line in file]

    fullpath = [os.path.join(directory, filename) for filename in filenames]

    return fullpath


################################################################################



# ConfigParser Object
config = configparser.ConfigParser()

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--files", help="\tFile containing scans", default="scans.txt")
parser.add_argument("--plot", help="\tMake pngs", action=argparse.BooleanOptionalAction, default=False)
args = parser.parse_args()

# Read config file for Data Paths
config.read('../../common/config.ini')
path = config.get('Paths', 'L09_path')

search_files = get_filenames(path, args.files)

# Initialize empty lists to accumulate data
glon_list = []
glat_list = []
Ta_list = []

for file in search_files:
    # get glon, glat, and calibrated spectra from each OTF file
    print("trying OTF file: ", file, end="   ")

    result = doStuff(file, args)

    if result is None:
        print("nothing returned")
        pass
    else:
        (glon, glat, Ta) = result

        for i in range(0,len(glon)):
            # why am i doing this test? document things better!
            if glat[i]<2 and glat[i]>-2:
                glon_list.append(glon[i])
                glat_list.append(glat[i])
                Ta_list.append(Ta[i])

# Convert lists to numpy arrays
glon = np.array(glon_list)
glat = np.array(glat_list)
Ta  = np.array(Ta_list)

# Regrid parameters
beam = 0.004 # beam width (deg)

# open a new blank FITS file
hdr = fits.Header()
hdr['NAXIS']   = 2
hdr['DATAMIN'] = min(Ta)
hdr['DATAMAX'] = max(Ta)
hdr['BUNIT']   = 'K (Ta*)     '

hdr['CTYPE1']  = 'GLON        '
hdr['CRVAL1']  = min(glon)
hdr['CDELT1']  = beam               # 1 arcmin beam
hdr['CRPIX1']  = 0                  # reference pixel array index
hdr['CROTA1']  = 0
hdr['CUNIT1']  = 'deg         '

hdr['CTYPE2']  = 'GLAT        '
hdr['CRVAL2']  = min(glat)
hdr['CDELT2']  = beam               # 1 arcmin beam
hdr['CRPIX2']  = 0                  # reference pixel array index
hdr['CROTA2']  = 0
hdr['CUNIT2']  = 'deg         '

hdr['OBJECT']  = 'NGC6334     '
hdr['GLON']    = min(glon)            # Fiducial is arbitrarily (glat,glon) min
hdr['GLAT']    = min(glat)

# Do the regridding
image = regrid(glon, glat, Ta, beam)

# Write the data cube and header to a FITS file
hdu = fits.PrimaryHDU(data=image, header=hdr)
hdu.writeto('my_data_image.fits', overwrite=True)



