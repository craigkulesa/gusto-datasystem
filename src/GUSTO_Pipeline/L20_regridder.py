import os
import sys
import numpy as np
from tqdm import trange, tqdm
from datetime import datetime

from scipy import interpolate
from astropy import units as u

from astropy.io import fits
from astropy.io.fits import Header
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic

import matplotlib.pyplot as plt

from .DataIO import loadSDFITS
from .flagdefs import *

from multiprocessing import Process, Queue

#logger = logging.getLogger('pipelineLogger')
num_workers = 4

directory = '/home/obs/data/GUSTO/level1/'

#build a 3D array for regridding: axis 1: GALLONG axis 2: GALLAT axis 3: VELOCITY
#build a second 3D array for counting denominator to divide by when done
spacing = 2/60 # 2 arcmin spatial spacing
vel_spacing = 2.0 # km/s
linename = "CII"

#GPS
#source_name="GPS"
#l_start=332.5
#l_end=363.5
#l_vector = np.arange(l_start,l_end,spacing)
#b_vector = np.arange(-1.5,1.5,spacing)
#v_vector = np.arange(-250,250,vel_spacing)

#NGC6334
#source_name="NGC6334"
#l_start=350.5
#l_end=354
#l_vector = np.arange(l_start,l_end,spacing)
#b_vector = np.arange(-1.5,2.5,spacing)
#v_vector = np.arange(-50,50,vel_spacing)
#t_max = 4

#NGC6357
source_name="NGC6357"
l_start=352.8
l_end=353.5
l_vector = np.arange(l_start,l_end,spacing)
b_vector = np.arange(0,1.5,spacing)
v_vector = np.arange(-50,16,vel_spacing)
t_max = 4

#LMC
#source_name="LMC"
#l_start=278.5
#l_end=281
#l_vector = np.arange(l_start,l_end,spacing)
#b_vector = np.arange(-32.1,-30,spacing)
#v_vector = np.arange(0,500,vel_spacing)

#NGC3603
#source_name="NGC3603"
#l_start=291
#l_end=292.2
#l_vector = np.arange(l_start,l_end,spacing)
#b_vector = np.arange(-1.0,0.0,spacing)
#v_vector = np.arange(-250,116,vel_spacing)
#t_max = 3

l_grid, b_grid, v_grid = np.meshgrid(l_vector,b_vector,v_vector,indexing='ij')
intensity_sum = np.zeros_like(v_grid,dtype='float64')
count_sum = np.zeros_like(v_grid,dtype='float64')

# find the closest l,b
# determine indicies to address
# compact vlsr to v_vector size
# add to temp_sump
# increment count_sum

def worker(input_q,output_q):    
    #print("Worker started")
    worker_intensity_sum = np.zeros_like(v_grid,dtype='float64')
    worker_count_sum = np.zeros_like(v_grid,dtype='float64')
    while True:
        input_filename = input_q.get()
        if input_filename is None:
            #print("Worker terminated")
            break  # Signal to terminate
        ifile = os.path.join(directory,input_filename)
        dsc = int(input_filename.split('_')[2].split('.')[0])
        spec, data, hdr, hdr1 = loadSDFITS(ifile, verbose=False)
        rowFlag = data['ROW_FLAG']
        n_spec, n_pix = spec.shape
        
        # new level 1 files are already in V_lsr, no longer need this code
        ## compute velocity                                                           
        npix    = hdr['NPIX']            
        #IF_pix  = hdr['CRPIX1']
        #IF_val  = hdr['CRVAL1']
        #IF_del  = hdr['CDELT1']
        #IF_freq = (np.arange(npix)-IF_pix)*IF_del+IF_val
        #VLSR    = hdr['VLSR']        
        #
        #IF_vlsr0= hdr['IF0']
        #vlsr    = (IF_vlsr0 - IF_freq)/line_freq*constants.c.value/1.e3 + VLSR # Vlsr in km/s 
        
        # new level 1 V_lsr import
        VLSR_pix  = hdr['CRPIX1']
        VLSR_val  = hdr['CRVAL1']
        VLSR_del  = hdr['CDELT1']
        vlsr = (np.arange(npix)-VLSR_pix)*VLSR_del+VLSR_val
        
        line_freq = hdr['LINEFREQ']
        
        if (linename == "CII") & (line_freq < 1900500):
            continue
        if (linename == "NII") & (line_freq > 1900500):
            continue
        
        
        #osel = np.argwhere((data['scan_type'] == 'OTF') & (data['ROW_FLAG']==0)).flatten()
        #osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & 0x60)==0) & ((data['MIXER']==5) | (data['MIXER']==8)) & (data['rms']<5)).flatten()
        if (linename == "CII"):
            #osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & 0x60)==0) & ((data['MIXER']==8))).flatten()
            osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & 0x60)==0) & ((data['MIXER']==5) | (data['MIXER']==8))).flatten()
        if (linename == "NII"):
            osel = np.argwhere((data['scan_type'] == 'OTF') & ((data['ROW_FLAG'] & 0x60)==0) & ((data['MIXER']==2) | (data['MIXER']==3) | (data['MIXER']==6))).flatten()
        
        if len(osel) <= 0:
            print('WARNING: No OTF spectra available in ',input_filename)
            # logger.warning('No OTF spectra available.')                           
        else:
            spec_OTF = np.squeeze(spec[osel,:])
            data_OTF = np.squeeze(data[osel])
            #mask = np.zeros(spec[osel,:].shape, dtype=bool)
            #indices = np.where(data_OTF['CHANNEL_FLAG'] & ChanFlags.SPUR)
            #mask[indices] = True
            #spec_OTF = np.ma.MaskedArray(spec[osel,:], mask)
            n_OTF, n_otfpix = spec_OTF.shape
            #x=np.arange(n_otfpix)
    
            # Instantiate for ra,dec->l,b transform                                      
            c_ra_dec = SkyCoord(ra=data_OTF['RA']*u.degree, dec=data_OTF['DEC']*u.degree, frame='icrs')
    
            #basecorr = np.zeros(spec_OTF.shape)
            #rmsotf = np.zeros(n_OTF)
            #rf = np.zeros(n_OTF)
            
            # empty arrays to fill                                                        
            
            c_l_b = c_ra_dec.transform_to(Galactic)    # transform to l,b
            
            wrapped_l = c_l_b.l.value[0]
            if wrapped_l<180: wrapped_l += 360
            #iterator.set_description(desc=f"{wrapped_l:.3f}")
            
            if wrapped_l < l_start:
                #print(f"Missed {wrapped_l:.3f}")
                continue
            if wrapped_l > l_end:
                #print(f"Missed {wrapped_l:.3f}")
                continue
            
            #if c_l_b.b.value[0] < c_l_b.b.value[n_OTF-1]:
            #    direction='P';
            #else:
            #    direction='N';
            
            #plt.clf()
            #        print('processing data ...')
            # calculate total power                                             
            for i0 in range(n_OTF):
                wrapped_l = c_l_b.l.value[i0]
                if wrapped_l<180: wrapped_l += 360
                l_index = (l_vector >= wrapped_l-spacing/2) & (l_vector < wrapped_l+spacing/2)
                b_index = (b_vector >= c_l_b.b.value[i0]-spacing/2) & (b_vector < c_l_b.b.value[i0]+spacing/2)
                #full_std_dev = np.ma.std(spec_OTF[i0,:])
                #weight = 1.0/(full_std_dev ** 2)
                weight = 1.0/(data_OTF['rms'][i0]**2)
                #if (l_index[10] & b_index[71]):
                #    print(dsc,i0,data_OTF['MIXER'][i0],hex(data_OTF['ROW_FLAG'][i0]),np.ma.sum(spec_OTF[i0,:]),full_std_dev,weight)
                #    plt.plot(vlsr,spec_OTF[i0])
                #    plt.show()
                for idx,vel in np.ndenumerate(v_vector):
                    v_index = np.abs(vlsr - vel) <= vel_spacing/2
                    worker_intensity_sum[l_index,b_index,idx] += weight*np.ma.sum(spec_OTF[i0,v_index])
                    worker_count_sum[l_index,b_index,idx] += weight*np.ma.count(spec_OTF[i0,v_index])

            #intensity=intensity_sum/count_sum
            #plt.clf()
            #plt.pcolormesh(l_grid[:,:,1],b_grid[:,:,1],intensity[:,:,65])
            #plt.show()
    print("Worker sending results...")
    output_q.put([worker_intensity_sum,worker_count_sum])

    

#filenames = ["CII_10000_24038_L10.fits"]

if __name__ == '__main__':
    send_q = Queue(num_workers)
    recv_q = Queue()
    processes = []
    for i in range(num_workers):
        p = Process(target=worker, args=(send_q,recv_q,))
        processes.append(p)
        p.start()

    iterator = tqdm(os.listdir(directory))
    for input_filename in iterator:
        #for input_filename in filenames:
        # check that the filename ends with .fits
        if not input_filename.endswith(".fits"):
            continue
        send_q.put(input_filename)
        

    #print("Sending termination...")
    for i in range(num_workers):
        send_q.put(None)  # Signal workers to terminate

    #print("Waiting for join...")
    #for p in processes:
    #    p.join()
    #print("Receiving results...")
    for i in range(num_workers):
        [worker_intensity_sum,worker_count_sum] = recv_q.get()
        intensity_sum += worker_intensity_sum
        count_sum += worker_count_sum
        #print(worker_count_sum[worker_count_sum>0].shape)
        
    #print("Waiting for join...")
    for p in processes:
        p.join()

    print("Generating outputs...")
    intensity=intensity_sum/(count_sum*vel_spacing)
#    start_range=50
#    end_range=200

        # open a new blank FITS file
    hdr = fits.Header()
    hdr['NAXIS']   = 3
    hdr['DATAMIN'] = np.nanmin(intensity)
    hdr['DATAMAX'] = np.nanmax(intensity)
    hdr['BUNIT']   = 'K           '

    hdr['CTYPE3']  = 'GLON        '
    hdr['CRVAL3']  = min(l_vector)
    hdr['CDELT3']  = spacing            # 2 arcmin beam
    hdr['CRPIX3']  = 0                  # reference pixel array index
    hdr['CROTA3']  = 0
    hdr['CUNIT3']  = 'deg         '

    hdr['CTYPE2']  = 'GLAT        '
    hdr['CRVAL2']  = min(b_vector)
    hdr['CDELT2']  = spacing            # 2 arcmin beam
    hdr['CRPIX2']  = 0                  # reference pixel array index
    hdr['CROTA2']  = 0
    hdr['CUNIT2']  = 'deg         '

    hdr['CTYPE1']  = 'VEL        '
    hdr['CRVAL1']  = min(v_vector)
    hdr['CDELT1']  = vel_spacing        # 2 km/s
    hdr['CRPIX1']  = 0                  # reference pixel array index
    hdr['CROTA1']  = 0
    hdr['CUNIT1']  = 'km/s       '

    hdr['OBJECT']  = f'{source_name}'
    hdr['GLON']    = min(l_vector)            # Fiducial is arbitrarily (glat,glon) min
    hdr['GLAT']    = min(b_vector)

    # Write the data cube and header to a FITS file
    hdu = fits.PrimaryHDU(data=intensity, header=hdr)
    datestring = datetime.today().strftime('%Y%m%d')
    hdu.writeto(f'{source_name}-{linename}-{datestring}.fits', overwrite=True)

    
    for idx in range(len(v_vector)-1):
        #print(idx,v_vector[idx])
        plt.figure(figsize=(2,2),dpi=100)
        plt.clf()
        output = plt.pcolormesh(l_grid[:,:,1],b_grid[:,:,1],intensity[:,:,idx],vmin=-1,vmax=t_max)
        plt.gca().invert_xaxis()
        plt.xlabel("Galactic Longitude",color='black')
        plt.ylabel("Galactic Latitude",color='black')
        plt.gca().set_facecolor('black')
        if (linename == "CII"):
            #plt.gca().title.set_text("line:CII mixers:8 vel:%.0f to %.0f km/s" % (v_vector[idx], v_vector[idx+1]))
            plt.gca().title.set_text("line:CII mixers:5,8 vel:%.0f to %.0f km/s" % (v_vector[idx], v_vector[idx+1]))
        if (linename == "NII"):
            plt.gca().title.set_text("line:NII mixers:2,3,6 vel:%.0f to %.0f km/s" % (v_vector[idx], v_vector[idx+1]))
        plt.gca().set_aspect('equal', adjustable='datalim')
        plt.colorbar(output)
        #plt.show()
        plt.savefig(f'{source_name}-{linename}-rms-{(idx):03d}.png',dpi=100)

