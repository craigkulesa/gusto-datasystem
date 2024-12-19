#
# Search RA,DEC and return scanIDs
#
import os
import glob
import time
import numpy as np
import datetime
from influxdb import InfluxDBClient

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5

import matplotlib.pyplot as plt
from PyAstronomy import pyasl

#################################################################################
global numi
global numj
numi = 0
numj = 0

username = ''
password = ''
database = 'gustoDBlp'
retention_policy = 'autogen'
bucket = f'{database}/{retention_policy}'

# Use influxdb for V1.0
# https://influxdb-python.readthedocs.io/en/latest/index.html
client = InfluxDBClient('localhost', 8086, '', '', 'gustoDBlp')

# half-beam size for searching l,b grid
size = 1*u.degree

# image coordinates
coords = "17h20m50.9s -36d06m54.0s" # NGC6334
c = SkyCoord(coords, unit=(u.hourangle, u.deg), frame='icrs')
lat = c.galactic.l.value
lon = c.galactic.b.value

#galactic image coordinates
gc = SkyCoord(l=lat*u.degree, b=lon*u.degree, frame='galactic')

# image size
l_img = 6*u.degree
b_img = 2*u.degree

# create l index
start_l  = gc.l.deg - l_img.value/2
end_l    = gc.l.deg + l_img.value/2
N_l      = int(l_img.value/(2*size.value))
l_indx   = np.linspace(start_l, end_l, N_l)

# create b index
start_b  = gc.b.deg - b_img.value/2
end_b    = gc.b.deg + b_img.value/2
N_b      = int(b_img.value/(2*size.value))
b_indx   = np.linspace(start_b, end_b, N_b)

# list of scanIDs in requested image
my_list=[]

for l in l_indx:
    numi += 1
    numj = 0
    for b in b_indx:
        numj += 1
        # find all points in udpPointing where we pointed at ra, dec
        gci = SkyCoord(l=l*u.degree, b=b*u.degree, frame='galactic')
        c = gci.transform_to(FK5)
        myquery = f'SELECT * FROM "udpPointing" WHERE RA< {(c.ra.value +size.to(u.deg).value)} AND \
                                                      RA> {(c.ra.value -size.to(u.deg).value)} AND \
                                                      DEC<{(c.dec.value+size.to(u.deg).value)} AND \
                                                      DEC>{(c.dec.value-size.to(u.deg).value)}' 
        points = client.query(myquery).get_points()

        # For loop over all of these pointings
        # POINTS contains a (time, scanID) for each pointing at (ra,dec)
        #tmp_list = []   # temporary list for (l,b,0/1) output
        #for point in points:
        #    tmp_list.append(point.get("scanID"))

        #my_list.append(tmp_list)
        for point in points:
            my_list.append(point.get("scanID"))

        # output of l,b and yes/no for udpPointing
        # used for plotting where we may have data
        #if tmp_list:
        #    print(l, b, 1)
        #elif not tmp_list:
        #    print(l, b, 0)


prefix = "ACS3_"
suffix = "_5_L09.fits"
scanID_list = sorted(list(set(my_list)))
filenames = ["{:s}{:05d}{:s}".format(prefix,int(num),suffix) for num in scanID_list]

# output sorted list of L1 fits filenames
print("\n".join(filenames))


