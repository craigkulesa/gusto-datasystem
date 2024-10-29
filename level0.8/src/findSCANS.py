#
# Search RA,DEC and return scanIDs
#
import os
import glob
import time
import numpy as np
import datetime
from influxdb import InfluxDBClient
#from influxdb_client import InfluxDBClient
from tqdm import tqdm

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
size = 0.5*u.degree

# image coordinate center
gc = SkyCoord(l=351.0*u.degree, b=0.0*u.degree, frame='galactic')
c = gc.transform_to(FK5)
# image size 
l_img = 6.0*u.degree
b_img = 2.0*u.degree


start_l  = gc.l.deg - l_img.value/2
end_l    = gc.l.deg + l_img.value/2
N_l      = int(l_img.value/(2*size.value))
l_indx   = np.linspace(start_l, end_l, N_l)

start_b  = gc.b.deg - b_img.value/2
end_b    = gc.b.deg + b_img.value/2
N_b      = int(b_img.value/(2*size.value))
b_indx   = np.linspace(start_b, end_b, N_b)
print(N_l)
print(N_b)
print(l_indx)
print(b_indx)

my_list=[]

for l in tqdm(l_indx, desc='Main Loop', leave=True):
    numi += 1
    numj = 0
    for b in tqdm(b_indx, desc=f"loop for l = {l}", leave=False):
        numj += 1
        # find all points in udpPointing where we pointed at ra, dec
        gci = SkyCoord(l=l*u.degree, b=b*u.degree, frame='galactic')
        c = gc.transform_to(FK5)
        myquery = f'SELECT * FROM "udpPointing" WHERE RA<{(c.ra.value +size.to(u.deg).value)}  AND \
                                                      RA>{(c.ra.value -size.to(u.deg).value)}  AND \
                                                     DEC<{(c.dec.value+size.to(u.deg).value)}  AND \
                                                     DEC>{(c.dec.value-size.to(u.deg).value)}' 
        points = client.query(myquery).get_points()

        # For loop over all of these pointings
        # POINTS contains a (time, scanID) for each pointing at (ra,dec)
        for point in points:
            my_list.append(point.get("scanID"))


prefix = "ACS5_"
suffix = "_6_L09.fits"
scanID_list = sorted(list(set(my_list)))
filenames = [f"{prefix}{num}{suffix}" for num in scanID_list]
print("\n".join(filenames))


