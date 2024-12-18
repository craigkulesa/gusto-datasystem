import gc
import os
import sys
import glob
import numpy as np
from astropy.io import fits

def doStuff(fits_file):
# Open file
    with fits.open(fits_file, memmap=True) as hdul:

        # Read in hdu[2] DATA_TABLE
        data = hdul['DATA_TABLE'].data

        # Mask of all REF type rows
        mask = data['scan_type']=='OTF'
        # All REF type scanID
        listOTF = hdul['DATA_TABLE'].data['scanID'][mask]
        if (len(list(set(listOTF))) == 1):
            otf = list(set(listOTF))[0]
        elif (len(list(set(listOTF))) == 0):
            print(fits_file, "error no otf!")
            new_file = f"{fits_file}.bad"
            os.rename(fits_file, new_file)
            return
        elif (len(list(set(listOTF))) > 1):
            print(fits_file, "error more than one otf")
            new_file = f"{fits_file}.bad"
            os.rename(fits_file, new_file)
            return

        # Mask of all REF type rows
        mask = data['scan_type']=='REF'
        # All REF type scanID
        listREF = hdul['DATA_TABLE'].data['scanID'][mask]
        if not listREF.size:
            print(fits_file, "error no REFs!")
            new_file = f"{fits_file}.bad"
            os.rename(fits_file, new_file)
            return

        # Mask of all HOT type rows
        mask = data['scan_type']=='HOT'
        # All REF type scanID
        listHOT = hdul['DATA_TABLE'].data['scanID'][mask]
        if not listHOT.size:
            print(fits_file, "error no HOTs!")
            new_file = f"{fits_file}.bad"
            os.rename(fits_file, new_file)
            return

        earliestREF = list(sorted(set(listREF)))[0]
        latestREF = list(sorted(set(listREF)))[len(list(set(listREF)))-1]

        if( earliestREF < otf < latestREF ):
            print(earliestREF, otf, latestREF, " OK", len(set(listREF)))
        else:
            print(earliestREF, otf, latestREF, " ERROR")
            new_file = f"{fits_file}.bad"
            os.rename(fits_file, new_file)

directory = "./"
partial = sys.argv[1]
fits_files = sorted(glob.glob(f"{directory}/{partial}.fits"))

for fits_file in fits_files:
    doStuff(fits_file)
    gc.collect()


