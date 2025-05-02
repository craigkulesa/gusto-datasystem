import os
import shutil
import glob
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from astropy.io import fits
from pybaselines import Baseline, utils
from astropy import constants as const
from astropy.table import Table
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from pprint import pprint
#from tqdm.notebook import tqdm
from tqdm import tqdm

plt.rcParams['font.size'] = 16


# identify the files for processing
opath = '/rdat/Projects/GUSTO/Data'
inDir = '/rdat/Projects/GUSTO/Data/level0.9'   #cfi['L08DataDir']
filter = 'ACS*L09.fits'
mtfile = os.path.join(opath,'master_table_20250109.csv')

scanRange = [19600, 19700]

sdirs = sorted(glob.glob(filter, root_dir=inDir))
print(sdirs[0])

dsc = []
dsc = [int(sdir.split('_')[1]) for sdir in sdirs]

sdirs.sort(key=lambda sdirs: dsc)

dseries = []
dfiles = []
for i,ds in enumerate(dsc):
    # if (ds >= scanRange[0]) & (ds <= scanRange[1]):
    dseries.append(ds)
    dfiles.append(sdirs[i])
# print(dseries)
# print(dfiles)


# dtype = [('obs_ID', '<i8'), ('scan_ID', '<i8'), ('obs_seq', '<i8'), ('obs_time', '<f8'), 
#          ('unixtime', '<i8'), ('target_name', '<U16'), ('target_ID', '<i8'), ('lon', '<f8'), 
#          ('lat', '<f8'), ('cframe', '<U16'), ('cal_ID1', '<i8'), ('cal_ID2', '<i8'), 
#          ('freq_cal_ID', '<i8'), ('RxBand', '<i8'), ('activeMixers', '<i8'), ('rowMixer', '<i8'), 
#          ('mixMode', '<U3'), 
#          ('sideBand', '<U3'), ('obsMode', '<U3'), ('obsType', '<U3'), ('processingLev', '<i8'), 
#          ('frequency', '<f8'), ('LOfreq', '<f8'), ('minVelocity', '<f8'), ('maxVelocity', '<f8'), 
#          ('velocityRes', '<f8'), ('exposTime', '<f8'), ('rms', '<f8'), ('Tsys', '<i8'), 
#          ('Ta_min', '<f8'), ('Ta_max', '<f8'), ('quality', '<i8'), ('specFlag', '<i8'), 
#          ('file_name', '<U256'), ('delivDate', '<i8'), ('procVers', '<U16')]

with open(mtfile, 'w', encoding='utf-8') as f:
    ostr = 'scanID,line,band,sctype,umxs,ra,dec,l,b,LINEFREQ,SYNTFREQ,SYNTMULT,VLSR,IF0,ELEVATON,dfile\n'
    f.write(ostr)
    # f.write("Hello, World!\n")

    # dfile = dfiles[0]
    for dfile in tqdm(dfiles):

        with fits.open(os.path.join(inDir,dfile)) as hdu:
            #hdu.info()
            data0  = hdu[0].data
            hdr0   = hdu[0].header
            data1 = hdu[1].data
            hdr1   = hdu[1].header
        
            umixers = np.unique(data1['MIXER'])
            umxs = ' '.join(map(str, np.unique(data1['MIXER'])))
            ra = data1['RA'][data1['scan_type']=='OTF']
            dec = data1['DEC'][data1['scan_type']=='OTF']
            cc0 = SkyCoord(np.median(ra)*u.deg, np.median(dec)*u.deg, frame='icrs')
            l0 = cc0.galactic.l.deg
            b0 = cc0.galactic.b.deg
            scanIDs = np.unique(data1['scanID'])
            scanID = int(dfile.split('_')[1])
            scantypes = np.unique(data1['scan_type'])
            if 'PS' in scantypes:
                sctype = 'PS'
            elif 'OTF' in scantypes:
                sctype = 'OTF'
            else:
                sctype = '?'
            band = hdr0['BAND']
            line = hdr0['LINE']
            
            # ostr = '%5i, %i, %3s, %3s, %5s, %12.6f, %12.6f'%(scanID, band, line, sctype, umxs, ra.mean(), dec.mean())
            ostr = '%i,%s,%i,%s,%s,%.6f,%.6f,%.6f,%.6f,%.1f,%.2f,%.0f,%.1f,%.1f,%.5f,%s\n'%(scanID, line, band, sctype, umxs, 
                                               ra.mean(), dec.mean(), l0, b0,
                                               hdr0['LINEFREQ'], hdr0['SYNTFREQ'], hdr0['SYNTMULT'], hdr0['VLSR'], hdr0['IF0'], hdr0['ELEVATON'],dfile)
            f.write(ostr)
            # print(ostr)

print('saved file: ', mtfile)

