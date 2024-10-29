import sys
import glob
import numpy as np
import configparser
from astropy.io import fits
from astropy import constants as const
from pybaselines import Baseline, utils
import matplotlib.pyplot as plt


'''
Goal:
   1)  Read all rows in all given fits files
   2)  Check baseband spectra against threshold per mixer/band for ringing
   3)  For REFs and HOTs, note failed rows in log
   4)  For OTFs, check Tsys from quasi-Level 1 

'''

resgood = np.empty((0,3))   # (std_dev, unixtime, mixer#)  STD_DEV of good rows of baseband spectra
resbad  = np.empty((0,3))   # (std_dev, unixtime, mixer#)  STD_DEV of bad  rows of baseband spectra
ta_std  = np.empty((0,3))   # (tsys,    unixtime, mixer#)  STD_DEV of Level 1
fitnes  = np.empty((0,3))   # (tsys,    unixtime, mixer#)  
tsys    = np.empty((0,3))   # (tsys,    unixtime, mixer#)


def find_mixer(line, mixer):
    #
    if(line=='NII'):
        mixers=[2,3,4,6]
    elif(line=='CII'):
        mixers=[2,3,5,8]

    # return index in mixer list of given mixer number
    if mixer in mixers:
        return mixers.index(mixer)
    else:
        return 5


def doTsys(hot, ref, thot):
    tsky = 35 # callen welton 1.4 THz
    tsky = 46 # callen welton 1.9 THz
    y_factor = np.zeros(hot.shape)
    y_factor = hot/ref

    # apply 30% backend non-linearity
    y_factor = ((y_factor-1)/1.3) + 1

    tsys_mix = (thot - tsky*y_factor[:])/(y_factor[:] - 1.0)

    return (tsys_mix[:].sum(axis=0))/len(tsys_mix[:])


def doTmix(tsys_mix, spec, ref):
    ta_mix = 2.*tsys_mix * (spec - ref)/ref

    return ta_mix[:]


def doStuff(scan, ta_std, fitnes, tsys, resgood, resbad):
    # open fits file
    hdu    = fits.open(scan, mode='update')
    # read header
    header = hdu[0].header
    line   = header['LINE']
    # read data_table
    data   = hdu[1].data
    spec   = data['spec'] 
    scan_type = data['scan_type']
    mixer     = data['MIXER']
    unixtime  = data['UNIXTIME']
    THOT      = data['THOT']
    ROW_FLAG  = data['ROW_FLAG']

    if 'HISTORY' not in header:
        print(scan, " No bad row flags in header .. flagging bad rows")
    else:
        history_list = header['HISTORY']
        if f'noisy rows flagged' in history_list:
            print(scan, " Row flagging already done .. stopping")
            #exit


    if(line=='NII'):
        x0=12
        x1=24
        l0=12
        mixers=[2,3,4,6]
        thresh=[.1, .1, .1, .1]
    elif(line=='CII'):
        x0=25
        x1=50
        l0=25
        mixers=[2,3,5,8]
        thresh=[.01, .1, .1, .0004]

    # Look through all of the spectra to set row_flags for noisy data
    xdata = np.arange(0,l0)
    data  = np.zeros(l0)
    for i in range(len(spec)):

        # compute standard deviation
        ydata = spec[i][x0:x1]
        z = np.polyfit(xdata, ydata, 5)
        p = np.poly1d(z)
        for j in range(l0):
           data[j] = ydata[j] - p(j)

        if(np.std(data)>thresh[find_mixer(line, mixer[i])]):
            ROW_FLAG[i] = 1
            resbad  = np.vstack([resbad,  [[np.std(data), unixtime[i], mixer[i]]]])
        else:
            ROW_FLAG[i] = 0
            resgood  = np.vstack([resgood, [[np.std(data), unixtime[i], mixer[i]]]])

        if(scan_type[i] == 'OTF'):
            # compute Tsys around zero vlsr
            # for the current row, irrespective of the ROW_FLAG, see what the Tmix would be

            # first, masks of ref,refhot,hot rows
            ref_mask = np.logical_and(mixer == mixer[i], scan_type == 'REF')
            refhot_m = np.logical_and(mixer == mixer[i], scan_type == 'REFHOT')
            hot_mask = np.logical_and(mixer == mixer[i], scan_type == 'HOT')

            # for now, skip the exact row, just average all the REF,REFHOT,HOT(s)
            refhot = spec[refhot_m,x0:x1].sum(axis=0)/len(spec[refhot_m,x0:x1])
            ref    = spec[ref_mask,x0:x1].sum(axis=0)/len(spec[ref_mask,x0:x1])
            tsysmix = doTsys(refhot, ref, THOT[i])
            tsys = np.vstack([tsys, [[tsysmix, unixtime[i], mixer[i]]]])

            hot = spec[hot_mask,x0:x1].sum(axis=0)/len(spec[hot_mask,x0:x1])
            ta_mix = doTmix(tsysmix, spec[i][x0:x1]/hot, ref/refhot)

            ta_std_mix = np.std(ta_mix, axis=0)              # std dev of spectral data
            ta_std = np.vstack([ta_std, [[ta_std_mix, unixtime[i], mixer[i]]]])

            fitnes_mix = np.std(data)
            fitnes = np.vstack([fitnes, [[fitnes_mix, unixtime[i], mixer[i]]]])

    # Write the change back to the fits file
    hdu[0].header.add_history('noisy rows flagged')
    hdu.flush()
    return(ta_std, fitnes, tsys, resgood, resbad)


# ConfigParser Object
config = configparser.ConfigParser()

# Read config file for Data Paths
config.read('config.ini')
paths=[]
paths.append(config.get('Paths', 'B1_path'))
paths.append(config.get('Paths', 'B2_path'))

partial = sys.argv[1]
search_files=[]
for path in paths:
   search_files+=sorted(glob.glob(f"{path}/{partial}.fits"))

# Debug
#print("\n".join(search_files))

for file in (search_files):
    ta_std, fitnes, tsys, resgood, resbad = doStuff(file, ta_std, fitnes, tsys, resgood, resbad)

plt.plot(ta_std[ta_std[:, 2] == 2])
plt.plot(ta_std[ta_std[:, 2] == 3])
plt.plot(ta_std[ta_std[:, 2] == 6])
plt.ylim((0, 1000))
