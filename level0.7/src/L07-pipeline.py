#!/usr/bin/env python3.11
from influxdb import InfluxDBClient
from astropy.io import fits
from multiprocessing import Pool
from functools import partial
import glob
import struct
from datetime import datetime
import sys
import os.path
import numpy as np
import configargparse
import argparse
sys.path.append("../../common/")
import flagdefs

acsPrefix = ['ACS5_', 'ACS3_']
outputPrefix = ['NII_', 'CII_']
bandPrefix = ['B1/', 'B2/']
RAD2DEG = 57.295779513 
C2K = 273.15
CDELT = [5000.0/511.0, 5000.0/1023.0] 
MULT = [108, 144]
sequencesFile = 'sequences.txt'


def flatten(xss):
    return [x for xs in xss for x in xs]


def makeSequences(input, output):
    print("Making new sequences file")
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

        
def makeFileGlob(startID, endID, bandNum, dirDataIn):
    fileList = []
    for scanID in range(startID, endID+1):
        fileList.append(glob.glob(dirDataIn+bandPrefix[bandNum-1]+acsPrefix[bandNum-1]+'*_'+str(scanID).zfill(5)+'.fits'))
    return(flatten(fileList))


def checkSequence(fileList):
    seqFlag = 0
    REF=0
    OTF=1
    HOT=2
    typeStr = ['REF', 'OTF', 'HOT']
    type = []
    scanID = []
    isHOTREF = []
    num = [0,0,0]
    
    for x in fileList:
        scanID.append(int(x[-10:-5]))
        type.append(x[-14:-11])
        if x[-14:-11] == 'REF':  
            isHOTREF.append(int(x[-10:-5]))  # add this scanID to list
    for y in range(0,3):
        num[y] = type.count(typeStr[y])

    if(num[REF] < 2):
        seqFlag |= flagdefs.SeqFlags.MISSING_REF
        if(num[REF] == 0):
            print("ERROR: no REFs found")
            seqFlag |= flagdefs.SeqFlags.NOREFS
    if(num[OTF] == 0):
        print("ERROR: no data!")
        seqFlag |= flagdefs.SeqFlags.NODATA
    if(num[HOT] != num[REF] + num[OTF]):
        seqFlag |= flagdefs.SeqFlags.MISSING_HOT
    return seqFlag,isHOTREF


def makeUDP(startID, stopID, dir):
    BUFSIZE=156
    output = []
    for file in range(startID, 1+stopID):
        filename = dir+'udp_'+str(file).zfill(5)+'.dat'
        stream = open(filename, "rb")

        while True:
            chunk = stream.read(BUFSIZE)
            if not chunk:
                break 
            data=struct.unpack('<2I5f2I15d',chunk)
            output.append(list(data))  # convert tuple to list to modify the time
        stream.close()

    for entry in output:  # convert integer timespec to floating point unixtime
        entry[0] = entry[0]+entry[1]/1.0e9
    return np.array(output)


def getInflux(startTime, endTime, queryStr):
    startTime = int((startTime-600)*1.0e+09)
    endTime = int(endTime*1.0e+09)
    client = InfluxDBClient(host='localhost', port=8086, database='gustoDBlp')
    query = f"SELECT * FROM /^{queryStr}*/ WHERE time > {startTime:d} AND time < {endTime:d} ORDER BY time DESC LIMIT 1"
    return client.query(query)


def getIFinfo(startTime):
    try:
        data = np.genfromtxt(options.inpath+"IF.txt", delimiter=None, dtype=None, usecols=(0,1,3,5,6,7), encoding=None, names=True)
    except Exception as e:
        raise Exception(f"Error reading CSV file: {e}")

    for index,value in enumerate(data['time']):
        if value > startTime:
            idx = index-1
            break
    info = (data['IF0'][idx], data['VLSR'][idx], data['object'][idx], data['LO1'][idx], data['LO2'][idx])
    return info


def processFITS(input_files, output_file, bandNum, pointingStream, seqFlag, listREF):
    with fits.open(input_files[0]) as hdul1:
        nrows1 = hdul1[1].data.shape[0]
        columns = hdul1[1].columns
        header = hdul1[0].header
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
    RA = data['RA']
    DEC = data['DEC']
    (udpTime, udpRA, udpDEC) = (pointingStream[:,0], pointingStream[:,5], pointingStream[:,6])
    for i in range(nrows):  # interpolate for RA and DEC
        RA[i] = RAD2DEG*np.interp(UNIXTIME[i], udpTime, udpRA)
        DEC[i] = RAD2DEG*np.interp(UNIXTIME[i], udpTime, udpDEC)

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
    prevItem = None
    newHOT = False
    index = 0
    psat = []
    vmon = []
    imon = []
    gmon = []
    fuzz = 15.0
    # because the data are not in order, we need to make an array of HK data first
    for item in TYPE:
        curTime = UNIXTIME[index]
        if item == 'HOT' and prevItem != 'HOT':
            newHOT=True
            hotTime = curTime
        elif item == 'HOT' and prevItem == 'HOT' and curTime-hotTime > 30:
            newHOT=True
            hotTime = curTime
        prevItem = item
        if newHOT:
            psat.append(curTime)
            vmon.append(curTime)
            imon.append(curTime)
            gmon.append(curTime)
            results=getInflux(curTime-fuzz, curTime+fuzz, 'PSatI_B'+str(bandNum)+'M')
            for i in pList[bandNum-1]:
                qStr = 'PSatI_B'+str(bandNum)+'M'+str(i)
                points = results.get_points(measurement=qStr)
                for point in points:
                    psat.append(point['cur'])
            results=getInflux(curTime-fuzz, curTime+fuzz, 'bias')
            for i in pList[bandNum-1]:
                qStr = 'biasVltB'+str(bandNum)+'M'+str(i)
                points = results.get_points(measurement=qStr)
                for point in points:
                    vmon.append(point['volts'])
                qStr = 'biasCurB'+str(bandNum)+'M'+str(i)
                points = results.get_points(measurement=qStr)
                for point in points:
                    imon.append(point['cur'])
            results=getInflux(curTime-fuzz, curTime+fuzz, 'B'+str(bandNum)+'_GMONI_')
            for i in pList[bandNum-1]:
                qStr = 'B'+str(bandNum)+'_GMONI_'+str(i)
                points = results.get_points(measurement=qStr)
                for point in points:
                    gmon.append(point['cur'])
            newHOT=False
        index += 1

    psat = np.array(psat)
    psat = psat.reshape(int(psat.shape[0]/5.0), 5)
    vmon = np.array(vmon)
    vmon = vmon.reshape(int(vmon.shape[0]/5.0), 5)
    imon = np.array(imon)
    imon = imon.reshape(int(imon.shape[0]/5.0), 5)
    gmon = np.array(gmon)
    gmon = gmon.reshape(int(gmon.shape[0]/5.0), 5)
    
    # Now we can loop through all the rows and assign the closest HK in time.
    # While we are here, update data types and row flags.
    index=0
    imonRange = (30.0, 40.0)
    for item in TYPE:
        if item == 'HOT' and SCANID[index] in listREF:
            TYPE[index] = 'REFHOT'
        curTime = UNIXTIME[index]
        closest = np.argmin(np.abs(psat[:,0] - curTime))
        PSAT[index] = psat[closest][pIdx[bandNum-1][MIXER[index]]]
        VMON[index] = vmon[closest][pIdx[bandNum-1][MIXER[index]]]-0.35
        IMON[index] = imon[closest][pIdx[bandNum-1][MIXER[index]]]
        GMON[index] = gmon[closest][pIdx[bandNum-1][MIXER[index]]]
        if IMON[index] < min(imonRange) or IMON[index] > max(imonRange):
            ROWFLAG[index] |= flagdefs.RowFlags.MIXER_MISPUMPED
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
    results=getInflux(min(UNIXTIME), max(UNIXTIME), "HK_TEMP")
    if results.items() == []:
        print("Nothing returned!")


    for i in range(1,17):
        qStr = "HK_TEMP"+str(i)
        points = results.get_points(measurement=qStr)
        for point in points:
            if qStr == 'HK_TEMP1':
                header['HKSCANID'] = (int(point['scanID']), 'scanID of slow housekeeping data')
            header[hktemp_names[i-1]] = (C2K+point['temp'],  hktemp_descs[i-1])

    lotemp_names = ("UNUSED","B1_SYNTH","UNUSED","UNUSED", "B1M5_AMP","UNUSED","UNUSED","UNUSED","B1_PWR_1","B1_PWR_2","B1_PWR_3","B1_PWR_4","B2_UCTRL","B2MLTDRV","UNUSED","UNUSED", "UNUSED","B2AVA183","B1M5MULT","B2M5_AMP","B2_PWR_1","B2_PWR_2", "B2_PWR_3","B2_PWR_4")
    lotemp_descs = ("UNUSED", "B1 LO Synthesizer", "UNUSED", "UNUSED", "B1 LO Spacek amplifier Ch5", "UNUSED", "UNUSED", "UNUSED", "B1 LO Pwr Box 1",  "B1 LO Pwr Box 2", "B1 LO Pwr Box 3", "B1 LO Pwr Box 4", "B2 MK66FX uCTRL", "B2 Mult Driver", "UNUSED", "UNUSED", "UNUSED", "B2 LO X-band Amplifier",  "B1 LO final tripler Ch5",  "B2 LO Spacek amplifier Ch5",  "B2 LO Pwr Box 1", "B2 LO Pwr Box 2", "B2 LO Pwr 3", "B2 LO Pwr 4")
    results=getInflux(min(UNIXTIME), max(UNIXTIME), "B1_AD590_")
    if results.items() == []:
        print("Nothing returned!")

    for i in range(0,12):
        qStr = "B1_AD590_"+str(i)
        points = results.get_points(measurement=qStr)
        for point in points:
            if lotemp_names[i] != "UNUSED":
                header[lotemp_names[i]] = (C2K+point['temp'],  lotemp_descs[i])

    results=getInflux(min(UNIXTIME), max(UNIXTIME), "B2_AD590_")
    if results.items() == []:
        print("Nothing returned!")

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

    results=getInflux(min(UNIXTIME), max(UNIXTIME), "DT670")
    if results.items() == []:
        print("Nothing returned!")
    i=0
    for name in cryo_names:
        qStr = 'DT670_'+name
        points = results.get_points(measurement=qStr)
        for point in points:
            header['T_'+cryo_names[i]] = (point['cryo'], cryo_descs[i])
        i += 1

    header['GON_ALT'] = (np.mean(pointingStream[:,2]), 'Gondola Altitude (m)')
    header['GON_LON'] = (RAD2DEG*np.mean(pointingStream[:,3]), 'Gondola Longitude (deg)')
    header['GON_LAT'] = (RAD2DEG*np.mean(pointingStream[:,4]), 'Gondola Latitude (deg)')

    info = getIFinfo(np.mean(UNIXTIME))
    header['OBJECT'] = (info[2], 'Name of the target object')
    header['IF0'] = (info[0], 'IF Frequency (MHz) of catalog velocity')
    header['SYNTFREQ'] = (info[bandNum+2], 'Synthesizer frequency (MHz)')
    header['SYNTMULT'] = (MULT[bandNum-1], 'Synthesizer multiplier')
    header['VLSR'] = (info[1], 'Commanded catalog velocity (km/s LSR)')
    
    # modify some last second items
    now = datetime.now()
    header['PROCTIME'] = now.strftime("%Y%m%d_%H%M%S")
    header['DLEVEL'] = 0.7
    header['SEQ_FLAG'] = seqFlag

    # now write it out
    hdu = fits.BinTableHDU.from_columns(columns, nrows=nrows)
    for colname in columns.names:
        hdu.data[colname][:] = data[colname]

    hduList = fits.HDUList([fits.PrimaryHDU(header=header), hdu])
    hduList.writeto(output_file, overwrite=True)
    return


def clear_folder(folder_path):
    print('Erasing contents of '+folder_path)
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")



def scanSequenceFile(input, options):
    inRange = []
    with open(input) as f:
        for line in f:
            words = line.split()
            (seqID, startID, endID, obsType, catName) = words[0:5]
            if(int(options.scanid[0]) < int(endID) and int(options.scanid[1]) >= int(startID)):
                inRange.append(words)
    return inRange



def processSequence(options, line):
    bandNum=int(options.band)
    start = int(options.scanid[0])
    stop = int(options.scanid[1])
    dirUDP = options.inpath + 'udp/'
    dirDataIn = options.inpath + 'level0.5/'
    dirDataOut = options.outpath + 'level0.7/' + bandPrefix[bandNum-1]

    (seqID, startID, endID, obsType, catName) = line[0:5]
    print("Processing sequence", seqID, "from", startID, "-", endID)
    fileList = makeFileGlob(int(startID), int(endID), bandNum, dirDataIn)
    seqFlag,listREF = checkSequence(fileList)
    if(seqFlag < flagdefs.SeqFlags.NOREFS):  # acceptable, process it
        print("seq ", seqID, " OK, flag is ", seqFlag)            
        pointingStream = makeUDP(int(startID), int(endID), dirUDP)
        output_file = dirDataOut+outputPrefix[bandNum-1]+seqID+'_'+startID+'.fits'
        processFITS(fileList, output_file, bandNum, pointingStream, seqFlag, listREF)
    else:
        # skip it because seqFlag says it's unusable
        print("seq", seqID, "NOT OK, flag is", seqFlag)



if __name__ == '__main__':
    p = configargparse.ArgParser(default_config_files=['./L07-config', '~/.L07-config'])
    p.add('-c', '--config', required=False, is_config_file=True, help='config file path')
    p.add('-e', '--erase', required=False, action=argparse.BooleanOptionalAction, help='erase contents of output folder before starting')
    p.add('-j', '--cpus', required=False, help='set number of CPUs to use')
    p.add('-b', '--band', required=True, help='process Band 1 or 2 data')
    p.add('-i', '--inpath', required=True, help='path to input udp and level0.5 folders')
    p.add('-o', '--outpath', required=True, help='path to output Level 0.7 files')
    p.add('-s', '--scanid', required=True, help='scanID range', nargs=2)
    options = p.parse_args()
    
    print(options)
    print(p.format_values())    # useful for logging where different settings came from

    bandNum=int(options.band)        
    dirDataOut = options.outpath + 'level0.7/' + bandPrefix[bandNum-1]    
    input = options.outpath + sequencesFile
    
    if not os.path.exists(input):
        makeSequences(options.inpath+"dataLog.txt", options.outpath + sequencesFile)
    else:
        print("Sequences file exists, skipping step...")

    if not os.path.exists(dirDataOut):
        os.makedirs(dirDataOut)
        print(f"Directory {dirDataOut} created.")
    else:
        print(f"Reusing directory {dirDataOut}")
        if options.erase:
            clear_folder(dirDataOut)

    # scan to find the elements in the right scanID range
    inRange = scanSequenceFile(input, options)
    # now multiprocess the processed list
    if(options.cpus):
        pool = Pool(processes=int(options.cpus))
    else:
        pool = Pool()
    pool.map(partial(processSequence, options), inRange)
