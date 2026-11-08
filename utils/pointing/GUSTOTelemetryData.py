import numpy as np
import struct
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import EarthLocation,SkyCoord
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import AltAz
import datetime



azoffset_org = np.array([[-0.0459,-0.0466,-0.0470,-0.0464,-0.0155,-0.0156,-0.0157,-0.0157],
                         [0.0155,0.0156,0.0157,0.0157,0.0459,0.0466,0.0470,0.0454],
                         [-0.0529,-0.0445,-0.0354,-0.0283,-0.0137,-0.0056,0.0025,0.1090]])

eloffset_org = np.array([[-0.0466,-0.0155,0.0158,0.0462,-0.0469,-0.0158,0.0157,0.0465],
                         [-0.0469,-0.0158,0.0157,0.0465,-0.0466,-0.0155,0.0158,0.0462],
                         [-0.0572,-0.0177,0.0211,0.0607,-0.0654,-0.0258,0.0130,0.0525]])

fplane = np.zeros((4,), dtype=[('band','i4',4), ('mx','i4',4), ('acs','i4'), ('unit', 'i4')])
fplane[0] = ([2,2,2,2], [2,3,5,8], 3, 4) # B2 M2, M3, M5, M8   ACS 3  UNIT 4
fplane[1] = ([2,3,3,3], [6,2,3,4], 4, 5) # B2 M6, B3(M2 M3 M4) ACS 4  UNIT 5
fplane[2] = ([1,1,1,1], [2,3,4,6], 5, 6) # B1 M2, M3, M4, M6   ACS 5  UNIT 6
fplane[3] = ([3,3,3,3], [5,6,7,8], 6, 7) # B3 M5, M6, M7, M8   ACS 6  UNIT 7

# fiducial pixel is B2M8
azoffset = azoffset_org - azoffset_org[1,7]
eloffset = eloffset_org - eloffset_org[1,7]


class GUSTOTelemetryData:
    
    def __init__(self, tlfile, verbose=False):
        """Initializing class to handle the UDP/telemetry data from GUSTO.
        
        """

        BUFSIZE=156
        tsec = []
        tnsec = []
        alt = []       # altitude in meters
        lon = []       # position on
        lat = []       # Earth (in radian)
        dec = []       # declination (in radian)
        ra = []        # right ascension (in radian)
        l = []
        b = []
        fixNum = []     # star camera measurement number; if valid: >0
        starCount = []  # number of stars recognized
        q_adj = []
        q_attT =[]
        q_att = []
        telescope = []

        # read the database file

        stream = open(tlfile, "rb")
        while True:
            chunk = stream.read(BUFSIZE)
            if not chunk:
                break 
            # example struc (Payload Software ICD 7503-9035):
            # (1705877081, 564653361, 
            # 39443.0, -1.3347253799438477, -1.4172226190567017, 
            # 4.344618320465088, -0.8075085282325745, 
            # 3196, 31, 
            # -7.60478209754778e-06, -3.103714617491704e-05, -7.893887675775035e-05, 0.9999999963737582,   # adjustment quaternion
            # -0.8961781432746596, 0.19305917620925384, -0.13594816575702523, 0.37564742275901214,         # attitudeTaken
            # 0.8993521966113858, -0.23244267644512367, -0.04042201500265141, 0.3681060843295423,          # attitude
            # 0.0019392415296856042, 0.0035148919505684595, 0.9999919424059702)                            # telescope vector       
            data=struct.unpack('<2I5f2I15d', chunk)
            tsec.append(data[0])
            tnsec.append(data[1])
            alt.append(data[2])
            lon.append(data[3])
            lat.append(data[4])
            ra.append(data[5])
            dec.append(data[6])
            fixNum.append(data[7])
            starCount.append(data[8])
            quaternion_adjusted = [data[9],data[10],data[11],data[12]]
            q_adj.append(quaternion_adjusted)
            quaternion_Taken = [data[13],data[14],data[15],data[16]]
            q_attT.append(quaternion_Taken)
            quaternion_attitude = [data[17],data[18],data[19],data[20]]
            q_att.append(quaternion_attitude)
            telescope_vector = [data[21],data[22],data[23]]
            telescope.append(telescope_vector)
            
            #print(data)
            print(struct.unpack('<2I5f2I15d',chunk)) 
        stream.close()

        # determine which pointings are valid pointings
        valid = np.zeros(len(fixNum), dtype=bool)
        valid[np.array(fixNum)>0] = 1

        # applu valid pointings selection to all data
        self.ttime = (np.array(tsec)+np.array(tnsec)/1E9)[valid]
        self.alt = (np.array(alt))[valid]
        self.lon = (np.array(lon))[valid]
        self.lat = (np.array(lat))[valid]
        self.ra = (np.array(ra))[valid]
        self.dec = (np.array(dec))[valid]
        self.tlfile = tlfile
        self.fixNum = (np.array(fixNum))[valid]
        self.starCount = (np.array(starCount))[valid]
        self.q_adj = np.array(q_adj)
        self.q_attT = np.array(q_attT)
        self.q_att = np.array(q_att)
        self.tel_vector = np.array(telescope)
        
        if verbose:
            print('Read udp file: ', tlfile)
            print()
            print('Available data: ', np.unique(keys))


    def getlbPosition(self):
        gc = SkyCoord(self.ra*u.rad, self.dec*u.rad, frame='fk5')

        return gc.galactic.l.deg, gc.galactic.b.deg, self.ttime


    def getEarthLocation(self):
        return self.lon, self.lat, self.alt, self.ttime


    def getRaDec(self, units='rad'):
        if units == 'rad':
            return self.ra, self.dec, self.ttime
        elif units == 'deg':
            # cc = SkyCoord(self.ra*u.rad, self.dec*u.rad, frame='icrs')
            # return cc.ra.deg, cc.dec.deg, self.ttime
            return self.ra/np.pi*180., self.dec/np.pi*180., self.ttime


    def unix2time(self, utime):
        return [datetime.datetime.fromtimestamp(a) for a in utime]


    def getOffsetPos(self, band, mix, ra, dec, rdtime, lon, lat, alt, llatim, verbose=False):
        """Function calculating the mixer pointings using offsets from the fiducial pixel.
        """
        balloon = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=alt*u.m)
    
        obstime = self.unix2time(llatim)
    
        aa = AltAz(location=EarthLocation(lat=lat, lon=lon, height=alt*u.m), obstime=self.unix2time(llatim))
        
        coord = SkyCoord(ra*u.rad, dec*u.rad, frame='icrs')
        altaz = coord.transform_to(aa)
        
        # get beam offsets
        az_bm = azoffset[band, mix]
        alt_bm = eloffset[band, mix]
        
        # apply beam offsets
        naz = altaz.az.deg*u.deg + az_bm*u.deg
        nalt = altaz.alt.deg*u.deg + alt_bm*u.deg
    
        if verbose:
            print(obstime[0])
            print('Coord: ', coord[0], ra[0], dec[0])
            print(SkyCoord(ra[0]*u.rad, dec[0]*u.rad, frame='icrs'))
            print('Az: ', altaz.az.deg[0]*u.deg, az_bm*u.deg, naz[0])
            print('Alt: ', altaz.alt.deg[0]*u.deg, alt_bm*u.deg, nalt[0])
            print()
        
        ncc = SkyCoord(AltAz(az=naz, alt=nalt, obstime=obstime, location=balloon))
        nradec = ncc.transform_to('icrs')
    
        return nradec


    def getlbPositionMx(self, unit, dev, verbose=False):
        
        # translate unit/dev to band/mix
        # translate unit/dev to band/mix
        # unit: 1 - 7 and dev: 1 - 4
        band = fplane[unit-4]['band'][dev-1]
        mix = fplane[unit-4]['mx'][dev-1]
        if verbose:
            print('unit: %i, dev: %i -> band: %i, mx: %i'%(unit, dev, band, mix))

        balloon = EarthLocation(lat=self.lat, lon=self.lon, height=self.alt*u.m)
    
        obstime = self.unix2time(self.llatim)
    
        aa = AltAz(location=balloon, obstime=self.unix2time(self.llatim))
        
        coord = SkyCoord(self.ra*u.rad, self.dec*u.rad, frame='icrs')
        altaz = coord.transform_to(aa)
        
        # get beam offsets
        az_bm = azoffset[band, mix]
        alt_bm = eloffset[band, mix]
        
        # apply beam offsets
        naz = altaz.az.deg*u.deg + az_bm*u.deg
        nalt = altaz.alt.deg*u.deg + alt_bm*u.deg

        # calculate new l an b
        ncc = SkyCoord(AltAz(az=naz, alt=nalt, obstime=self.obstime, location=balloon))
        nlb = ncc.transform_to('galactic')
    
        return nlb.l.deg, nlb.b.deg, obstime


    def getPosBandMixer(self, band, mix, verbose=False):
        """Function calculating the mixer pointings using offsets from the fiducial pixel.
        """
        balloon = EarthLocation(lat=self.lat, lon=self.lon, height=self.alt*u.m)
    
        obstime = self.unix2time(self.ttime)
        aa = AltAz(location=balloon, obstime=obstime)
        
        coord = SkyCoord(self.ra*u.rad, self.dec*u.rad, frame='icrs')
        altaz = coord.transform_to(aa)
        
        # get beam offsets
        az_bm = azoffset[band-1, mix-1]
        alt_bm = eloffset[band-1, mix-1]
        
        # apply beam offsets
        naz = altaz.az.deg*u.deg + az_bm*u.deg
        nalt = altaz.alt.deg*u.deg + alt_bm*u.deg
    
        if verbose:
            print(obstime[0])
            print('Coord: ', coord[0], self.ra[0], self.dec[0])
            print(SkyCoord(self.ra[0]*u.rad, self.dec[0]*u.rad, frame='icrs'))
            print('Az: ', altaz.az.deg[0]*u.deg, az_bm*u.deg, naz[0])
            print('Alt: ', altaz.alt.deg[0]*u.deg, alt_bm*u.deg, nalt[0])
            print()
        
        # calculate new R.A. and Dec.
        ncc = SkyCoord(AltAz(az=naz, alt=nalt, obstime=obstime, location=balloon))
        nradec = ncc.transform_to('icrs')
        
    
        return nradec, self.ttime


    def getlbUnitDev(self, unit, dev, verbose=False):

        # translate unit/dev to band/mix
        band = fplane[int(unit-4)]['band'][int(dev-1)]
        mix = fplane[int(unit-4)]['mx'][int(dev-1)]
        
        nradec, obstime = self.getPosBandMixer(band, mix)
        
        ngal = nradec.transform_to('galactic')
        
        return ngal.l.deg, ngal.b.deg, obstime 


def unitdev2bandmix(unit, dev):
    """Function translating the unit/dev arrangement to the actual bands and mixers.
    """
    unit = np.array((unit,)).flatten()
    dev = np.array((dev,)).flatten()
    
    nn = np.array(unit).size
    bands = np.zeros(nn, dtype=int)
    mixers = np.zeros(nn, dtype=int)
    
    for i in range(nn):
        bands[i] = fplane[int(unit[i]-4)]['band'][int(dev[i]-1)]
        mixers[i] = fplane[int(unit[i]-4)]['mx'][int(dev[i]-1)]
    return bands, mixers


if __name__ == '__main__':
    from GUSTOTelemetryData import GUSTOTelemetryData
    import matplotlib.pyplot as plt

    plt.rcParams['axes.formatter.useoffset'] = False

    test = 3
    
    if test == 1:

        tlfile = '/Users/volkertolls/Projects/GUSTO/Data/Iridium/udps/udp_05800.dat'
        # tlfile = '/Users/volkertolls/Projects/GUSTO/Data/Iridium/udps/udp_01400.dat'
        tlfile = '/Users/volkertolls/Projects/GUSTO/Data/Iridium/udps/udp_07090.dat'
    
    
        gtd = GUSTOTelemetryData(tlfile)
        
        l, b, utime = gtd.getlbPosition()
        
        lmean = np.nanmean(l)
        bmean = np.nanmean(b)
        
        dl = (l - lmean) * 3600.
        db = (b - bmean) * 60.
        
        print()
        print(utime.min())
        print(utime.max())
        print(dl.min(), dl.max())
    
        fig, ax1 = plt.subplots(1, 1, sharex=True, figsize=(6,5))
        
        ax1.scatter(dl, db, marker='*', label='eta Car', linewidth=1.5)
        #ax1.plot(l, b, '*', label='eta Car', linewidth=1.5)
        ax1.set_ylabel('Galactic Latitude (arcmin)')
        ax1.set_xlabel('Galactic Longitude (arcsec)')
        
        fig.tight_layout(pad=0.1)
        #plt.savefig("plot.png", dpi=130.0, pad_inches=0.1, bbox_inches = 'tight')
        plt.show()

    elif test == 2:
        upath = '/Users/volkertolls/Projects/GUSTO/Data/Iridium/udps/'
        ifile = '/Users/volkertolls/Projects/GUSTO/Data/Iridium/udps/udp_09519.dat'
        
        gtd = GUSTOTelemetryData(ifile)
        lon, lat, alt, llatim = gtd.getEarthLocation()
        ra, dec, rdtime = gtd.getRaDec()
        radg, decdg, rdtime = gtd.getRaDec(units='deg')
        
        print(lon[0], lat[0], alt[0], llatim[0])
        print(ra[0], dec[0], rdtime[0])

        nradec = gtd.getOffsetPos(2, 2, ra, dec, rdtime, lon, lat, alt, llatim)
        
        print('Location: ', lon[0], lat[0], alt[0], EarthLocation(lat=lat[0], lon=lon[0], height=alt[0]*u.m))
        print('RA:  ', ra[0], nradec.ra.rad[0])
        print('Dec: ', dec[0], nradec.dec.rad[0])
        print()
        
        
        fig = plt.figure(figsize=(12,6))
        sp1 = plt.subplot()
        
        pl1 = plt.plot(radg, decdg, '*', color='red')
        pl2 = plt.plot(nradec.ra.deg, nradec.dec.deg, 'o', color='blue')
        
        plt.xlabel('R.A. (deg)')
        plt.ylabel('Declination (deg)')
        
        # print(nradec.ra.deg)
        # print(nradec.dec.deg)

    elif test == 3:
        upath = '/Users/volkertolls/Projects/GUSTO/Data/Iridium/udps/'
        ifile = '/Users/volkertolls/Projects/GUSTO/Data/Iridium/udps/udp_09519.dat'
        
        gtd = GUSTOTelemetryData(ifile)
        l, b, tutime = gtd.getlbUnitDev(4, 1)
        
        print(l)
        print(b)
        print(tutime)
        
