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

# this file has been created using gustoL09P/support/GUSTO_create_master_table.py
ifile = '/Users/volkertolls/Projects/GUSTO/Data/master_table_20241223.csv'

def gusto_get_lb_region_list(lrange=None, mt_file=None, mt=None, verbose=False):
    """Function to get list of scans for l,b-region..


    Parameters
    ----------
    mt_file : str
        Path to master table
    mt : astropy table
        master table provided
        master table can be created using script gustoL09P/support/GUSTO_create_master_table.py

    Returns
    -------
    list of scans for region for further processing.
    """
    
    if mt_file is None:
        mt_file = '/Users/volkertolls/Projects/GUSTO/Data/master_table_20241223.csv'
    
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
    
    scans = mt['scanID'][(lrange[0]<mt['l'])&(lrange[1]>mt['l'])&(mt['band']==2)]
    return scans

if __name__ == '__main__':
    
    cregion = 'RCW120'
    ra0 = '17h12m23.2s'
    dec0 = '−38d26m51.2s'
    cc0 = SkyCoord(ra0, dec0, frame='icrs')
    
    print('Region: ', cregion)
    print(cc0.galactic.l.deg, cc0.galactic.b.deg)
    
    lrange = [347.5, 348.75]
    print('l-range: ', lrange)

    scans = gusto_get_lb_region_list(lrange=lrange, verbose=True)
    
    print(scans.value.size)
    print(scans.value)
