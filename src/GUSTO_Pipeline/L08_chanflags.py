import os
from datetime import datetime
from astropy.io import fits
from astropy.time import Time
import numpy as np
from multiprocessing import Pool
from functools import partial
from .flagdefs import *
from .DataIO import *


def L08_Pipeline(args, scanRange):
    global commit_info
    dirDataOut = args.path+'level0.8/'
    prefix = ['NII_', 'CII_'] 
    os.makedirs(dirDataOut, exist_ok=True)
    if args.erase:
        clear_folder(dirDataOut)

    commit_info = runGitLog('0.8', 'L08_chanflags.py')  # do this only once     
    path = args.path+'level0.7/'
    sum_files = 0
    
    for band in args.band:
        files = makeFileGlob(path, prefix[int(band)-1], 'fits', scanRange)
        sum_files += len(files)
        if(args.cpus):
            pool = Pool(processes=int(args.cpus))
        else:
            pool = Pool()
        pool.map(partial(doStuff, options=args), files)
    return sum_files
    

def doStuff(scan, options):
    global commit_info

    hdu    = fits.open(scan)
    header = hdu[0].header
    band   = int(header['BAND'])-1 # indexed to 0
    data   = hdu[1].data
    spec   = data['DATA']
    mixer  = data['MIXER']
    nrow   = len(spec)
    ROW_FLAG  = data['ROW_FLAG']
    CHANNEL_FLAG = data['CHANNEL_FLAG']

    x0 = [12, 25]
    x1 = [24, 50]
    l0 = [12, 25]
    mixers = [[2,3,4,6], [2,3,5,8]]
    thresh = [[.1, .1, .1, .1], [.01, .1, .1, .0004]]
    spurs = [[[31,36],[66,69],[100,101],[133,136],[147,152],[165,168]], [[66,70],[134,136],[201,203],[266,272],[294,303],[331,335]]]
    oobs = [[[0,31],[410,511]], [[0,62],[820,1023]]]
    xdata = np.arange(0,l0[band])
    data  = np.zeros(l0[band])
    
    for i in range(nrow):
        # first set the channel flags
        CHANNEL_FLAG[i] = 0 # resetting flags -- this is safe because we are the first user
        for spur in spurs[band]:
            CHANNEL_FLAG[i][spur[0]:spur[1]]  |= ChanFlags.VARIABLE_SPUR
        for oob in oobs[band]:
            CHANNEL_FLAG[i][oob[0]:oob[1]]    |= ChanFlags.OOB 
        # set fringing rowflag by computing standard deviation
        ydata = spec[i][x0[band]:x1[band]]
        z = np.polyfit(xdata, ydata, 5)
        p = np.poly1d(z)
        for j in range(l0[band]):
            data[j] = ydata[j] - p(j)
        try:
            if(np.std(data) > thresh[band][mixers[band].index(mixer[i])]):
                ROW_FLAG[i] |=  (RowFlags.RINGING_BIT1 | RowFlags.RINGING_BIT0)  # set ringing bits
            else:
                ROW_FLAG[i] &= ~(RowFlags.RINGING_BIT1 | RowFlags.RINGING_BIT0)  # clear ringing bits
        except:
            print("L08: Error in thresholding -- this should NOT HAPPEN!")
            
    header['DLEVEL'] = 0.8
    header['COMMENT'] = commit_info
    tred = Time(datetime.now()).fits
    header.add_history('Level 0.8 processed at %s'%(tred))
    
    # Write the changes back to fits as new file
    hdu.writeto(scan.replace("level0.7", "level0.8").replace("L07", "L08"), overwrite=True)
    return None


