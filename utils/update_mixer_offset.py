import numpy as np
import glob
import os
import time

from astropy.table import Table
from astropy.io import fits
from astropy.table import vstack
from astropy import units as u
from astropy import constants


import struct
from astropy.coordinates import EarthLocation,SkyCoord
from astropy.time import Time
from astropy.coordinates import AltAz
import datetime
import argparse
import configargparse


def get_cal_mixer_offsets(calib_file):
    f = open(calib_file,'r')
    offsets=[]
    azoffs=[]
    eloffs=[]
    bm = []
    for line in f:
        line.strip('\n')
        if line[0] == 'B':
            # read offset data
            txt = line.split('\t')
            bmname=txt[0]
            azoff1 = float(txt[1])
            eloff1 = float(txt[2])
            mtype  = txt[3]
            off1 = [azoff1,eloff1]
            bm.append([bmname,mtype])
            azoffs.append(azoff1)
            eloffs.append(eloff1)
    #return offset in arcmin
    azoffs = np.array(azoffs)*60
    eloffs = np.array(eloffs)*60
    bm = np.array(bm)
    #bm is band / mixer identifier and theory/measured
    # BxMy  THEORY/AS_MEASRURED
    f.close()
    return( azoffs, eloffs , bm )

def unix2time(unixtime, utime):
    return [datetime.datetime.fromtimestamp(a) for a in utime]

def history_test(hdr,history_phrase,verbose=False):
    """
    
    test FITS header history key words or a specfic phrase.
    If phrase exists in header return True,  if no HISTORY or phrase not found return False.
    
    """
    if 'HISTORY' not in hdr:
        if verbose: 
            print(f'HISTORY not in {hdr}: Applying {history_phrase}')
        return(False)
    else: 
        history_list = hdr.get('HISTORY')
        if f'{history_phrase}' in history_list:
            if verbose: 
                print(f'{history_phrase} present')
            return(True)
        else:
            if verbose:
                print(f'{history_phrase} not found')
            return(False)




def update_mixer_offset(foff,hdr0, ra, dec,otime, band, mix, calib_file,verbose=False):
    """
    Function Applies a fractional offset to each mixer for band 1
    Input:  fraction_upate, fraction of original offsets to be added for new offset (azimuth, altitude)
    Input:  ra,dec,otime of sample postion:  in EXT 1 of FITS file
            lat,lon altitude of balloon  in header EXT[0] of 0.7 product
            band: 1 or 2 N+ or C+
            mix:  mixer (1-8) in band
    Output: updated RA/DEC of mixer at unixtime
    Output: hdr0 is updated with HISTORY field 
    
    
    """
    balloon = EarthLocation(lat=hdr0['GON_LAT']*u.deg, lon=hdr0['GON_LON']*u.deg, height=hdr0['GON_ALT']*u.m)

    #obstime = self.unix2time(self.ttime)
    #obstime = unix2time(unixtime)
    faz = float(foff[0])
    falt= float(foff[1])

    aa = AltAz(location=balloon, obstime=otime)

    #coord = SkyCoord(self.ra*u.rad, self.dec*u.rad, frame='icrs')
    coord = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    altaz = coord.transform_to(aa)

    #Open the Mixer offset calibration file
    azoffset, aloffset, bm = get_cal_mixer_offsets(calib_file)
    # get beam offsets:

    azoffm=[]
    aloffm=[]
    #for ix, mx in enumerate(mix):
    mName = f'B{band}M{mix}'
    indx = np.argwhere(bm[:,0] == mName).flatten()

    #use the last entry to get AS_MEASURED values if present
    azoffm.append( azoffset[indx[-1]] )
    aloffm.append( aloffset[indx[-1]] )

    #convert to decimal degrees
    azoffm = np.array(azoffm)/60.0
    aloffm = np.array(aloffm)/60.0
        
        



    # apply beam offsets
    if band == 1:
        naz = (altaz.az.deg  - azoffm  + faz * azoffm) * u.deg
        nalt =(altaz.alt.deg - aloffm  + falt* aloffm) * u.deg
    else:
        #don't update band 2
        naz = altaz.az.deg*u.deg 
        nalt = altaz.alt.deg*u.deg 
        


    # calculate new R.A. and Dec.
    ncc = SkyCoord(AltAz(az=naz, alt=nalt, obstime=otime, location=balloon))
    nradec = ncc.transform_to('icrs')
    return(nradec)

def main(args=None,verbose=True):
    if args==None:
        # Create the input parser
        my_parser = configargparse.ArgumentParser(prog='update_mixer_offset',
                                            usage='%(prog)s filetype directory ',
                                            description='Update RA/DEC positions in GUTSO data file')
        my_parser.version = "Version 0.0.1 (8 Dec 2025) "
        my_parser.add_argument('-v', action='version')

        my_parser.add('-o',
                               metavar='--foff',
                               nargs = 2, 
                               help = 'fractional offset from calibration (AZ, ALT)')
        my_parser.add('-r',
                               metavar='--fileroot',
                               type=str,
                               required=True,
                               help='Name of data file to be updated.')
        my_parser.add('-d',
                               metavar='--directory',
                               required = False,
                               default = '/data/scratch/GUSTO/gusto-datasystem/Data/level1/',
                               type=str,
                               help='Path to data directory')

        my_parser.add('-c',
                               metavar='--calib_file',
                               default = '/data/scratch/GUSTO/gusto-datasystem/calib/cal_offsets.txt',
                               type=str,
                               help='Path to calibration directory')
        args = my_parser.parse_args()
    

    
    #print(my_parser.format_values())
    #print(args)
    fileroot=args.r
    directory =args.d
    calib_file = args.c
    foff = np.array(args.o)
    print(f'{directory}/{fileroot}*.fits')

    to_process_files = sorted(glob.glob(f'{directory}/{fileroot}*.fits'))
    to_process_files = np.array(to_process_files).flatten()
    #print(to_process_files.shape)
    utc0=np.datetime64("1970-01-01T00:00:00",'s') 
    for infile in to_process_files:

        hdu=fits.open(infile,mode = 'update')
        hdr0=hdu[0].header

        # Set the history phrase to be added to FITS header
        history_phrase = f'Mixer Offsets Updated for Band 1'

        in_history = history_test(hdr0,history_phrase,verbose=True)
        
        if not in_history:
            # Collect needed parameter to apply a mixer offset postion
            hdr1=hdu[1].header
            data=hdu[1].data
    
            #extract data arrays from fits file
            mixer = data['mixer']
            unixtime=data['unixtime']
            # 'time' of observation
            otime = np.median(unixtime)
            ra=data['RA']
            dec=data['DEC']
            mixers = data['mixer']

            mixer = np.unique(mixers)
            
    
            band   = int(hdu[0].header['BAND'])
            
            #median time of observation time
            obstime = utc0 + np.timedelta64(int(otime*1000),'ms')

            for mx in mixer:
                
                qmx = data['mixer'] == mx
                mra = data['RA'][qmx]
                mdec= data['DEC'][qmx]
                #apply offsets for each mixer
                nradec = update_mixer_offset(foff,hdr0, mra, mdec, obstime, band, mx, calib_file,verbose=False)
                newra  = nradec.ra.deg
                newdec = nradec.dec.deg
        
                
                #Write back to fits array
                data['RA'][qmx] = np.array(newra)
                data['DEC'][qmx]= np.array(newdec)
    
            #Write back to file
            print(f'Updating {infile}')
            hdu[0].header.add_history(history_phrase)
            hdu.flush()
        else:
            print(f'Offsets already applied: {infile} not updated)')
        hdu.close()   
        #print(to_process_files.shape)
        #print(infile)

if __name__ == "__main__":
    main()

