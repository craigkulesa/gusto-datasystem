#!/usr/bin/env python3.11
import os
import glob
import shutil
from astropy.io import fits


def doCheckTHOT(scan):
    with fits.open(scan) as hdu:
        header = hdu[0].header
        try:
            object = header['THOT']
        except Exception as e:
            print("missing THOT:", scan)
            return
    return None


if __name__ == "__main__":
    path = '/home/obs/data/GUSTO/level0.8'
    files = sorted(glob.glob(f"{path}/*.fits"))

    for file in files:
        doCheckTHOT(file)
