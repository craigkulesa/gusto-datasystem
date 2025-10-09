from influxdb import InfluxDBClient
from astropy.io import fits
from astropy.time import Time
from multiprocessing import Pool
from functools import partial
from struct import unpack
from datetime import datetime
from tqdm import tqdm

import os
import numpy as np
import subprocess

from .flagdefs import *
from .DataIO import *
from .Logger import *

logger = logging.getLogger('pipelineLogger')

def makeSequences(input, output):
    logger.info("Making new sequences file")
    newSeq = True
    startID = oldscanID = oldobj = obj = oldseqType= ''
    seqType = 'APS'
    seqNum = 0
    with open(output, 'w') as seqFile:
        with open(input) as f:
            for line in f:
                words = line.split()
                (scanID, obsType, mode, target)=(words[3], words[4], words[5], words[-1])
                if(scanID == oldscanID or int(scanID) == 0): # catch dupes or nulls
                    continue
                if(obsType == 'REF' and newSeq == True): # new target, reset startID
                    newSeq = False
                    startID = int(scanID)
                elif obsType == 'OTF':
                    seqType = obsType
                elif obsType == 'SRC' and mode == 'ABS_OBS':
                    seqType = 'APS'
                elif(obsType == 'REF'):
                    if(int(scanID) == startID + 1):
                        startID = int(scanID)  # back to back REFs, only use the last one
                    else:
                        endID = int(scanID)
                        seqNum+=1
                        print(f"{seqNum:05d}", f"{startID:05d}", f"{endID:05d}", seqType, obj, file=seqFile)
                        startID = endID  # make this REF the starting ref for the next seq                
                if(mode == 'SPIRAL' or mode == 'BASIC'):  # skip these and throw out the sequence
                    newSeq = True
                elif(obsType != 'REF'):
                    obj = target
                    if(obj != oldobj and oldobj != ''):
                        newSeq = True  # New object -- flush old sequence to list and restart
                        if(obsType == 'OTF'):
                            endID = int(scanID)
                            seqType = obsType
                            seqNum+=1
                            print(f"{seqNum:05d}", f"{startID:05d}", f"{endID:05d}", seqType, obj, file=seqFile)
                        elif(obsType == 'SRC' and mode == 'ABS_OBS'):
                            endID = int(scanID)-1
                            seqType = 'APS'
                            seqNum+=1
                            print(f"{seqNum:05d}", f"{startID:05d}", f"{endID:05d}", oldseqType, oldobj, file=seqFile)
                oldscanID = scanID
                oldobj = obj
                oldseqType = seqType

        # end-of-file flush of last sequence
        endID = int(scanID)
        seqType = obsType
        seqNum+=1
        print(f"{seqNum:05d}", f"{startID:05d}", f"{endID:05d}", seqType, obj, file=seqFile)


def checkSequence(fileList):
    seqFlag = 0
    REF=0
    OTF=1
    HOT=2
    APS=3
    typeStr = ['REF', 'OTF', 'HOT', 'SRC']
    type = []
    scanID = []
    isHOTREF = []
    num = [0,0,0,0]
    foundUpperREF = False
    foundLowerREF = False
    
    for x in fileList:
        scanID.append(int(x[-14:-9]))
        type.append(x[-18:-15])
        if x[-18:-15] == 'REF':  
            isHOTREF.append(int(x[-14:-9]))  # add this scanID to list
    for y in range(0,4):
        num[y] = type.count(typeStr[y])

    if(num[REF] == 0):
        logger.info("ERROR: no REFs found")
        seqFlag |= SeqFlags.NOREFS
    elif(num[REF] > 0 and num[REF] < 2 and num[OTF] > 0):
        index = type.index('OTF')
        OTFscanID = scanID[index]
        refList= [i for i,e in enumerate(type) if e == 'REF']
        for item in refList:
            if scanID[item] < OTFscanID:
                foundLowerREF = True
            elif scanID[item] > OTFscanID:
                foundUpperREF = True
        if foundLowerREF == False:
            seqFlag |= SeqFlags.MISSING_LEADING_REF
        elif foundUpperREF == False:
            seqFlag |= SeqFlags.MISSING_TRAILING_REF
        else:
            logger.info("ERROR: inconsistent REF count")
    if(num[OTF] == 0 and num[APS] == 0):
        logger.info("ERROR: no data!")
        seqFlag |= SeqFlags.NODATA
    if(num[HOT] != num[REF] + num[OTF]):
        seqFlag |= SeqFlags.MISSING_HOT
    return seqFlag,isHOTREF


def makeUDP(startID, stopID, dir):
    BUFSIZE=156
    output = []
    for file in range(startID, 1+stopID):
        filename = dir+'udp_'+str(file).zfill(5)+'.dat'
        try:
            with open(filename, "rb") as stream:
                while True:
                    chunk = stream.read(BUFSIZE)
                    if not chunk:
                        break 
                    data=unpack('<2I5f2I15d',chunk)
                    output.append(list(data))  # convert tuple to list to modify the time
        except Exception as error:
            logger.info(error)
    for entry in output:  # convert integer timespec to floating point unixtime
        entry[0] = entry[0]+entry[1]/1.0e9
    if len(output) == 0:
        logger.info("ERROR: No UDP data available for this sequence")
        return np.array([0])
    return np.array(output)


def getInflux(startTime, endTime, queryStr, seqFlag, getAll=False, lookBack=1800):
    startTime = int((startTime-lookBack)*1.0e+09)
    endTime = int((endTime+lookBack/40)*1.0e+09)
    client = InfluxDBClient(host='localhost', port=8086, database='gustoDBlp')
    if getAll == False:
        query = f"SELECT * FROM /^{queryStr}*/ WHERE time > {startTime:d} AND time < {endTime:d} ORDER BY time DESC LIMIT 1"
        results = client.query(query)
    else:
        query = f"SELECT * FROM /^{queryStr}*/ WHERE time > {startTime:d} AND time < {endTime:d} ORDER BY time DESC"
        results = client.query(query, epoch='ms')
    if results.items() == []:  
        logger.info("WARNING:  Nothing from influx query", queryStr, "returned.")
        if queryStr == 'B2_AD590_' or queryStr == 'B1_AD590_':
            seqFlag |= SeqFlags.MAYBE_HUNG_LO
        else:
            seqFlag |= SeqFlags.MISSING_HK
    return results, seqFlag


def getIFinfo(startTime, path):
    try:
        data = np.genfromtxt(path+"IF.txt", delimiter=None, dtype=None, usecols=(0,1,3,5,6,7), encoding=None, names=True)
    except Exception as e:
        raise Exception(f"Error reading CSV file: {e}")

    for index,value in enumerate(data['time']):
        if value > startTime:
            break
    info = (data['IF0'][index-1], data['VLSR'][index-1], data['object'][index-1], data['LO1'][index-1], data['LO2'][index-1])
    return info


# check if this is a known bad scanID given a tuple of bad ranges
def isKnownBad(number, ranges):
    for start, end in ranges:
        if start <= number <= end:
            return True
    return False


# split and rearrange an influx query into a table with time and mixer data in columns
def splitConcatenate(array, dim=4):
    array = np.array(array)
    chunks = np.split(array, dim)
    newarray = np.concatenate(chunks, axis=1)
    chop = list(range(2, newarray.shape[1], 2))
    newarray = np.delete(newarray, chop, axis=1)
    return newarray


def processFITS(data_path, input_files, output_file, bandNum, pointingStream, seqID, seqFlag, listREF, catName):
    RAD2DEG = 57.295779513 
    C2K = 273.15
    CDELT = [5000.0/511.0, 5000.0/1023.0] 
    MULT = [108, 144]

    with fits.open(input_files[0]) as hdul1:
        nrows1 = hdul1[1].data.shape[0]
        columns = hdul1[1].columns
        header = hdul1[0].header
        if header['DLEVEL'] != 0.5:
            logger.info("ERROR: Input data is not at level 0.5.  Exiting.")
            return
        data = {col.name: [] for col in columns}
        for col in columns:
          data[col.name].extend(hdul1[1].data[col.name])
    
        for i in range(1, len(input_files)):
            with fits.open(input_files[i]) as hdul2:
                for col in columns:
                    data[col.name].extend(hdul2[1].data[col.name])
                    
    nrows = len(data[columns[0].name])
    
    # The level 0.7 data is now merged! Now we modify table and header contents
    # first, let's get interpolated RA, DEC for every row
    UNIXTIME = data['UNIXTIME']
    if len(pointingStream) > 1:
        (udpTime, udpRA, udpDEC) = (pointingStream[:,0], pointingStream[:,5], pointingStream[:,6])
        data['RA'] = [RAD2DEG*np.interp(t, udpTime, udpRA) for t in UNIXTIME]
        data['DEC'] = [RAD2DEG*np.interp(t, udpTime, udpDEC) for t in UNIXTIME]
        
    # Now let's update the "fast" housekeeping telemetry in the binary table
    MIXER = data['MIXER']
    TYPE = data['scan_type']
    ROWFLAG = data['ROW_FLAG']
    PSAT = data['PSat']
    VMON = data['Vmon']
    IMON = data['Imon']
    GMON = data['Gmon']
    SCANID = data['scanID']
    pList = ((2,3,4,6), (2,3,5,8))
    pIdx = (0,0,1,2,3,0,4,0,0,0), (0,0,1,2,0,3,0,0,4,0)
    index = 0
    psat = []
    vmon = []
    imon = []
    gmon = []

    # because the data are not in time order, we need to make an array of HK data first
    results,seqFlag = getInflux(min(UNIXTIME), max(UNIXTIME), 'PSatI_B'+str(bandNum)+'M', seqFlag, getAll=True)
    for i in pList[bandNum-1]:
        qStr = 'PSatI_B'+str(bandNum)+'M'+str(i)
        points = results.get_points(measurement=qStr)
        for point in points:
            psat.append((float(point['time']/1000.0), point['cur']))
    results,seqFlag = getInflux(min(UNIXTIME), max(UNIXTIME), 'bias', seqFlag, getAll=True)
    for i in pList[bandNum-1]:
        qStr = 'biasVltB'+str(bandNum)+'M'+str(i)
        points = results.get_points(measurement=qStr)
        for point in points:
            vmon.append((float(point['time']/1000.0), point['volts']))
        qStr = 'biasCurB'+str(bandNum)+'M'+str(i)
        points = results.get_points(measurement=qStr)
        for point in points:
            imon.append((float(point['time']/1000.0), point['cur']))
    results,seqFlag = getInflux(min(UNIXTIME), max(UNIXTIME),'B'+str(bandNum)+'_GMONI_', seqFlag, getAll=True)
    for i in pList[bandNum-1]:
        qStr = 'B'+str(bandNum)+'_GMONI_'+str(i)
        points = results.get_points(measurement=qStr)
        for point in points:
            gmon.append((float(point['time']/1000.0), point['cur']))

    skipHKinsert = False
    try:
        psat = splitConcatenate(psat)
        vmon = splitConcatenate(vmon)
        imon = splitConcatenate(imon)
        gmon = splitConcatenate(gmon)
    except Exception as e:
        skipHKinsert = True
        logger.info("WARNING: ", e, "Some columns will have zeroed bias data.")

    # Now we can loop through all the rows and assign the closest HK in time.
    # While we are here, update data types and row flags.
    index=0
    imonRange = (31.0, 41.0)
    for item in TYPE:
        if item == 'HOT' and SCANID[index] in listREF:
            TYPE[index] = 'REFHOT'
        curTime = UNIXTIME[index]
        if skipHKinsert == False:
            closest = np.argmin(np.abs(psat[:,0] - curTime))
            PSAT[index] = psat[closest][pIdx[bandNum-1][MIXER[index]]]
            closest = np.argmin(np.abs(vmon[:,0] - curTime))
            VMON[index] = vmon[closest][pIdx[bandNum-1][MIXER[index]]]-0.35
            closest = np.argmin(np.abs(imon[:,0] - curTime))
            IMON[index] = imon[closest][pIdx[bandNum-1][MIXER[index]]]
            closest = np.argmin(np.abs(gmon[:,0] - curTime))
            GMON[index] = gmon[closest][pIdx[bandNum-1][MIXER[index]]]
        if bandNum == 2 and isKnownBad(SCANID[index], B2Unlocked):
            ROWFLAG[index] |= RowFlags.LO_SYNTH_UNLOCKED
        if IMON[index] < min(imonRange) or IMON[index] > max(imonRange):
            ROWFLAG[index] |= RowFlags.MIXER_MISPUMPED
            if IMON[index] > 60 or IMON[index] < 0.0:
                ROWFLAG[index] |= RowFlags.MIXER_UNPUMPED
        index += 1
    
    
    # Now onto the primary HDU, lots to fill in here
    header['CUNIT1'] = ('MHz', 'Spectral unit')
    header['CRPIX1'] = (0.0, 'Index location')
    header['CRVAL1'] = (0.0, 'Start of spectra (MHz)')
    header['CDELT1'] = (CDELT[bandNum-1], 'Channel width (MHz)')

    hktemp_names = ("CRADLE02","CRYCSEBK","CRYOPORT","CALMOTOR","CRADLE03","QAVCCTRL",
                    "COOLRTRN","FERADIAT", "CRYCSEFT","CRADLE04","THOT","OAVCCTRL",
                    "COOLSUPL","CRADLE01","EQUILREF","SECONDRY")
    hktemp_descs = ("Cradle 2 temp", "Cryostat Back temp", "Cryostat pumpout port temp", 
		   "Calibration flip mirror motor temp", "Cradle 3 temp",
                   "QCL AVC Cryocooler CTRL temp", "Cooling Loop Return temp",
                   "Front End Radiator temp",  "Cryostat Left Side temp", "Cradle 4 temp",
                   "Calibration load temp", "OVCS AVC Cryocooler CTRL temp",
                   "Cooling Loop Supply temp",  "Cradle 1 temp", "Equilibar Reference temp",
                   "Secondary temp")
    results,seqFlag = getInflux(min(UNIXTIME), max(UNIXTIME), "HK_TEMP", seqFlag)
    for i in range(1,17):
        qStr = "HK_TEMP"+str(i)
        points = results.get_points(measurement=qStr)
        for point in points:
            if qStr == 'HK_TEMP1':
                header['HKSCANID'] = (int(point['scanID']), 'scanID of slow housekeeping data')
            header[hktemp_names[i-1]] = (C2K+point['temp'],  hktemp_descs[i-1])

    lotemp_names = ("UNUSED","B1_SYNTH","UNUSED","UNUSED", "B1M5_AMP","UNUSED","UNUSED","UNUSED","B1_PWR_1","B1_PWR_2","B1_PWR_3","B1_PWR_4","B2_UCTRL","B2MLTDRV","UNUSED","UNUSED", "UNUSED","B2AVA183","B1M5MULT","B2M5_AMP","B2_PWR_1","B2_PWR_2", "B2_PWR_3","B2_PWR_4")
    lotemp_descs = ("UNUSED", "B1 LO Synthesizer", "UNUSED", "UNUSED", "B1 LO Spacek amplifier Ch5", "UNUSED", "UNUSED", "UNUSED", "B1 LO Pwr Box 1",  "B1 LO Pwr Box 2", "B1 LO Pwr Box 3", "B1 LO Pwr Box 4", "B2 MK66FX uCTRL", "B2 Mult Driver", "UNUSED", "UNUSED", "UNUSED", "B2 LO X-band Amplifier",  "B1 LO final tripler Ch5",  "B2 LO Spacek amplifier Ch5",  "B2 LO Pwr Box 1", "B2 LO Pwr Box 2", "B2 LO Pwr 3", "B2 LO Pwr 4")
    results,seqFlag = getInflux(min(UNIXTIME), max(UNIXTIME), "B1_AD590_", seqFlag)
    for i in range(0,12):
        qStr = "B1_AD590_"+str(i)
        points = results.get_points(measurement=qStr)
        for point in points:
            if lotemp_names[i] != "UNUSED":
                header[lotemp_names[i]] = (C2K+point['temp'],  lotemp_descs[i])

    results,seqFlag = getInflux(min(UNIXTIME), max(UNIXTIME), "B2_AD590_", seqFlag)
    for i in range(12,24):
        qStr = "B2_AD590_"+str(i-12)
        points = results.get_points(measurement=qStr)
        for point in points:
            if lotemp_names[i] != "UNUSED":
                header[lotemp_names[i]] = (C2K+point['temp'],  lotemp_descs[i])

    cryo_names = ("IS", "IVCS", "LNA", "MIXER","OS", "OVCS", "QCL", "TANK")
    cryo_descs = ("temp inner shield (K)", "temp inner vapor shield (K)", "temp LNA (K)",
                  "temp mixers (K)", "temp outer shield (K)", "temp outer vapor shield (K)",
                  "temp QCL (K)", "temp liquid-He tank (K)")

    results,seqFlag = getInflux(min(UNIXTIME), max(UNIXTIME), "DT670", seqFlag)
    i=0
    for name in cryo_names:
        qStr = 'DT670_'+name
        points = results.get_points(measurement=qStr)
        for point in points:
            header['T_'+cryo_names[i]] = (point['cryo'], cryo_descs[i])
        i += 1

    if len(pointingStream) > 1:
        header['GON_ALT'] = (np.mean(pointingStream[:,2]), 'Gondola Altitude (m)')
        header['GON_LON'] = (RAD2DEG*np.mean(pointingStream[:,3]), 'Gondola Longitude (deg)')
        header['GON_LAT'] = (RAD2DEG*np.mean(pointingStream[:,4]), 'Gondola Latitude (deg)')

    results,seqFlag = getInflux(min(UNIXTIME)-86400, max(UNIXTIME)+30, "gonVel", seqFlag)
    points = results.get_points(measurement='gonVel')
    for point in points:
        header['GON_VEL'] = (point['vel'], 'Gondola Velocity to LSR (km/s)')
        
    info = getIFinfo(np.mean(UNIXTIME), data_path)
    header['OBJECT'] = (catName, 'Name of the target object')
    header['IF0'] = (info[0], 'IF Frequency (MHz) of catalog velocity')
    header['SYNTFREQ'] = (info[bandNum+2], 'Synthesizer frequency (MHz)')
    header['SYNTMULT'] = (MULT[bandNum-1], 'Synthesizer multiplier')
    header['VLSR'] = (info[1], 'Commanded catalog velocity (km/s LSR)')
    
    # modify or add some last second items
    tred = Time(datetime.now()).fits
    header['DLEVEL'] = 0.7
    header['SEQ_ID'] = (int(seqID), 'Sequence ID')
    header['SEQ_FLAG'] = seqFlag
    header['COMMENT'] = commit_info
    header.add_history('Level 0.7 processed at %s'%(tred))

    # now write it out
    hdu = fits.BinTableHDU.from_columns(columns, nrows=nrows, name="DATA_TABLE")
    for colname in columns.names:
        hdu.data[colname][:] = data[colname]

    hduList = fits.HDUList([fits.PrimaryHDU(header=header), hdu])
    hduList.writeto(output_file, overwrite=True)
    return


def scanSequenceFile(input, options):
    inRange = []
    with open(input) as f:
        for line in f:
            words = line.split()
            (seqID, startID, endID, obsType, catName) = words[0:5]
            if(int(options.scanid[0]) < int(endID) and int(options.scanid[1]) >= int(startID)):
                inRange.append(words)
    return inRange


def L07_Pipeline(args):
    global commit_info
    dirDataOut = args.path + 'level0.7/'
    input = args.path + 'sequences.txt'

    commit_info = runGitLog('0.7', 'L07_telemetry.py')  # lookup git commit info only once
    
    if not os.path.exists(input) or os.path.getsize(input) == 0:
        makeSequences(args.path+"dataLog.txt", args.path + sequencesFile)
    else:
        logger.info("Sequences file seemingly exists, skipping step...")

    os.makedirs(dirDataOut, exist_ok=True)
    if args.erase:
        clear_folder(dirDataOut)

    # scan to find the elements in the right scanID range
    inRange = scanSequenceFile(input, args)
    # now multi-process the processed list
    if(args.cpus):
        n_procs = int(float(args.cpus)/2)  # make it half of what we else everywhere else
        if n_procs < 1:  # but not less than 1
            n_procs=1
        pool = Pool(processes=n_procs)
        logger.info('Number of cores used for processing: %i\n'%(n_procs))
    else:
        pool = Pool()
    list(tqdm(pool.imap_unordered(partial(processSequence, args), inRange), total=len(inRange), colour='yellow', leave=False))
    return len(inRange)


def processSequence(options, line):
    global commit_info
    inputPrefix = ['ACS5_', 'ACS3_']
    outputPrefix = ['NII_', 'CII_']
    dirUDP = options.path + 'udp/'
    dirDataIn = options.path + 'level0.5/'
    dirDataOut = options.path + 'level0.7/'
    (seqID, startID, endID, obsType, catName) = line[0:5]
    options.scanid[0] = startID 
    options.scanid[1] = endID
    
    for band in options.band:
        logger.info(f"Processing Band {band} sequence {seqID} from {startID} - {endID}")
        scanRange = [int(x) for x in options.scanid]
        fileList = makeFileGlob(dirDataIn, inputPrefix[int(band)-1], 'fits', scanRange)
        seqFlag,listREF = checkSequence(fileList)
        if(seqFlag < SeqFlags.NOREFS):  # acceptable, process it
            pointingStream = makeUDP(int(startID), int(endID), dirUDP)
            output_file = dirDataOut+outputPrefix[int(band)-1]+seqID+'_'+startID+'_L07.fits'
            processFITS(options.path, fileList, output_file, int(band), pointingStream, seqID, seqFlag, listREF, catName)
        else:
            # skip it because seqFlag says it's unusable
            logger.info(f"Sequence {seqID} NOT OK, flag is {seqFlag}")

