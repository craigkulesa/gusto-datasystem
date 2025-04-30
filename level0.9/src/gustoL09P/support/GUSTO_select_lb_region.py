"""Function to get list of scans for selected region.
NGC 3603:

PS: 6586-6778, 19016
OTF: 8791-9115; 9504-9896; 10291-10395; 19012; 19020-19132;     2602-3784 (radec)


RCW 120:

SOFIA coverage:
l = [348, 348.5]
b = [0.29, 0.65]


PS: 9464-9500; 10685-10729; 11139-12095;

OTF
"""

from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import tarfile
import os.path
from tqdm import tqdm

# this file has been created using gustoL09P/support/GUSTO_create_master_table.py
ifile = '/rdat/Projects/GUSTO/Data/master_table_20250109.csv'

def gusto_get_lb_region_list(lrange=None, band=1, mt_file=None, mt=None, verbose=False):
    """Function to get list of scans for l,b-region..


    Parameters
    ----------
    lrange : tuple of floats
        range of data in Galacitic Longitude
    band : int
        receiver band for data selection
    mt_file : str
        Path to master table
    mt : astropy table
        master table provided
        master table can be created using script gustoL09P/support/GUSTO_create_master_table.py
    verbose : boolean
        chatty flag

    Returns
    -------
    list of scans for region for further processing.
    """
    
    if mt_file is None:
        mt_file = '/rdat/Projects/GUSTO/Data/master_table_20250109.csv'
    
    if mt == None:
        mt = Table().read(mt_file, delimiter=',')
        
    if verbose:
        print(mt.dtype)
        print()
    
    cc = SkyCoord(mt['ra']*u.deg, mt['dec']*u.deg, frame='icrs')
    mt['l'] = cc.galactic.l.deg
    mt['b'] = cc.galactic.b.deg
    
    if verbose:
        print(mt[0:20])
    
    scans = mt['scanID'][(lrange[0]<mt['l'])&(lrange[1]>mt['l'])&(mt['band']==band)]
    return scans

# /Users/volkertolls/Projects/GUSTO/Data/level0.9/ACS3_23300_L09.fits
def make_tarfile(tfile, scans, band):
    if band==1:
        root = 'ACS5'
    elif band==2:
        root = 'ACS3'
    else:
        print('Wrong receiver band selected!')
        return None
        
    cregion = '_'.join(tfile.split('/')[-1].split('_')[0:2])
    print(cregion)
    
    with tarfile.open(tfile, "w:gz") as tar:
        for scan in tqdm(scans):
            fname = '%s_%05i_L%s.fits'%(root, scan, level.replace('.',''))
            tar.add(os.path.join(inDir,fname), arcname=os.path.join(cregion,fname))
            #print('added: ', fname)


if __name__ == '__main__':
    
    cregion = 'RCW120'
    ra0 = '17h12m23.2s'
    dec0 = '−38d26m51.2s'
    cc0 = SkyCoord(ra0, dec0, frame='icrs')
    
    print('Region: ', cregion)
    print(cc0.galactic.l.deg, cc0.galactic.b.deg)
    
    lrange = [347.5, 348.75]
    print('l-range: ', lrange)

    
    band = 2
    
    level = '09'
    level = '1.0'
    
    inDir = '/rdat/Projects/GUSTO/Data/level%s/'%(level)
    tfile = ('/rdat/Projects/GUSTO/Data/gusto_region_files/%s_B%i_%.1f_%.1f_L%s'%(cregion, band, lrange[0], lrange[1], level.replace('.',''))).replace('.','p') + '.tar.gz'
    
    scans = gusto_get_lb_region_list(lrange=lrange, band=band, verbose=True)
        
    print(scans.value.size)
    print(scans.value)
    print('creating data file ...')
    
    make_tarfile(tfile, scans, band)
    print('saved region data in: ', tfile)
    
