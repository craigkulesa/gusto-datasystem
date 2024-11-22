"""This module just holds initial parameters that are going 
to be used in the GUSTO pielines.

The initial parameters and values are collected in an excel spreadsheet.
The goal is to leave the parameters in the spreadsheet and update this
module from time to time automatically using a python function (tbd).

History:
Version 0.4 was a major rewrite of the class to enable better
    access to the data.
"""


__version__ = '0.4'
__date__ = '9/19/2024'
__updated__ = '20240919'
__author__ = 'V. Tolls, CfA | Harvard & Smithsonian'


import astropy.units as u
import astropy.constants as c
import pandas as pd
import numpy as np
from numpy.lib import recfunctions as rfn
from dataclasses import dataclass, asdict
from pprint import pprint
import os
import math
import inspect
import datetime
import warnings
import pkg_resources
warnings.filterwarnings('ignore')
import logging
log = logging.getLogger(__name__)


fheader, tail = os.path.split(inspect.stack()[0][1])
__pyfile__ = tail
__path__ = fheader


func_fullname = inspect.stack()[0][1]
func_path, func_name = os.path.split(func_fullname)


class GL09PParameters:
    
    parlist = []  # Class attribute to track instance names.

    def __init__(self, verbose=False, parfile=None):
        """
        """
        # exec("self.{} = '{}'".format(vname, vvalue))

        self.padd('__date__', __date__, '', 'created')
        self.padd('__version__', __version__, '', 'class version')
        self.padd('__updated__', __updated__, '', 'last update')
        self.padd('__pyfile__', __pyfile__, '', 'file containing this class')
        self.padd('__path__', __path__, '', 'path to this file')
        self.padd('__author__', __author__, '', 'author')
                

        
        if parfile is None:
            ifile = pkg_resources.resource_filename('gustoL1', 'Data/GUSTO_Baseline_Data_draft2.xlsx')
        else:
            ifile = parfile
        self.padd('__parfile__', ifile, '', 'file containing parameters')
            
        # with pd.ExcelFile(func_path+ifile) as gpxlsx:
        with pd.ExcelFile(ifile) as gpxlsx:
            # load all the worksheets (wss) into a list
            wss = {}
            sheetNames = gpxlsx.sheet_names
            for sheet_name in sheetNames:
                wss[sheet_name] = gpxlsx.parse(sheet_name)
            n_wss = len(wss)
            if verbose:
                print('n_wss: ', n_wss, type(wss))
                print()

            # =====================================================
            # work on the first sheet
            df1 = wss['GUSTO_PreLaunch']
            sz1 = df1.shape
            # we are polling through each row of the sheet and
            # eliminate the empty rows and rows starting with a "#"
            # entries per row:
            #      Parameter, n_Parms, Unit, Par1, Par2, Par3, Info
            # padd(vname, vvalue, vunit, vinfo)
            for i in range(sz1[0]):
                row = df1.loc[i].values
                test = False
                if isinstance(row[0], float):
                    # if type(row[0])==type(float('nan')):
                    test = float(row[0]) != float('nan')
                if isinstance(type(row[0]), type('nan')):
                    test = row[0][0] == '#'
                    
                if not test:
                    if row[1] == 1:
                        if type(row[2])==type(2.1):
                            row2 = ''
                        else:
                            if 'nan' in row[2]:
                                row2 = ''
                            else:
                                row2 = row[2]
                        
                        if isinstance(type(row[6]), type('')):
                            if row[6]=='nan':
                                row6 = ''
                            else:
                                row6 = row[6]
                        else:
                            row6 = ''
                        self.padd(row[0], row[3], row2, row6)
                    elif row[1] == 3:
                        if isinstance(type(row[2]), type('')):
                            if row[2][0]=='nan':
                                row2 = ''
                            else:
                                row2 = row[2]
                        else:
                            row2 = ''
                        if isinstance(type(row[6]), type('')):
                            if row[6][0]=='nan':
                                row6 = ''
                            else:
                                row6 = row[6]
                        else:
                            row6 = ''
                        self.padd(row[0], np.array([row[3], row[4], row[5]]), 
                                  row2, row6)

            # =====================================================
            # work on the detector sheets
            n_dpix = 8

            idtype = [('cband', 'a32'), ('band', int, n_dpix),
                      ('pidx', int, n_dpix),
                      ('az', float, n_dpix), ('el', float, n_dpix),
                      ('detpix', int, n_dpix), ('FWHM', float, n_dpix),
                      ('AFWHMas', float, n_dpix), ('AFWHMam', float, n_dpix), 
                      ('EFWHMas', float, n_dpix), ('EFWHMam', float, n_dpix),
                      ('eta_beam', float, n_dpix), 
                      ('eta_point', float, n_dpix), ('eta_ap', float, n_dpix),
                      ('Tsys', float, n_dpix), ('Trms', float, n_dpix)]
            udtype = [('cbandUnit', 'a32'), ('bandUnit', 'a32'),
                      ('pidxUnit', 'a32'),
                      ('azUnit', 'a32'), ('elUnit', 'a32'),
                      ('detpixUnit', 'a32'), ('FWHMUnit', 'a32'),
                      ('AFWHMasUnit', 'a32'), ('AFWHMamUnit', 'a32'), 
                      ('EFWHMasUnit', 'a32'), ('EFWHMamUnit', 'a32'),
                      ('eta_beamUnit', 'a32'), 
                      ('eta_pointUnit', 'a32'), ('eta_apUnit', 'a32'),
                      ('TsysUnit', 'a32'), ('TrmsUnit', 'a32')]

            # read the detector parameters
            df = wss['Det_Offsets']
            sz = df.shape
            detpar = np.zeros((24,sz[1]))
            for i in range(detpar.shape[0]):
                detpar[i,:] = df.loc[i+1].values
            if verbose:
                print('sz: ', sz) 
                print('columns: ', df.columns.values)
                print()
                # print('df type:', type(df))
                # print(detpar)
            # column names
            urow = df.loc[0].values
            detparNames = np.ndarray((urow.size), dtype='a32')
            detparNames = [idtype[i][0] for i in range(urow.size)]
            print('detparNames: ', detparNames)
            # get the units for the columns
            detparUnits = np.ndarray((urow.size), 'a32')
            for i in range(urow.size):
                if type(urow[i]) == type(float):   #.encode()=='nan':
                    print('nan')
                    detparUnits[i] = '' # u.Unit()
                else:
                    # detparUnits[i] = u.Unit(str(urow[i]))
                    detparUnits[i] = u.Unit(urow[i])
            print('Units: ', detparUnits)
#             # loop over rows
#             for i in range(1, sz[0]):
#                 row = df.loc[i].values
#                 #if verbose:
#                 print(i, row)
#                 # loop over column index
#                 for j in range(1, row.size):
#                     detpar[int(row[0]-1)][j][(i-1)%8] = row[j]
#             self.detpar = detpar  
            
#             print()
#             print(detparNames)
#             print('el: ', detpar[0][2])
#             print('az: ', detpar[0][3])
#             print()
            self.cadd('rec')
            self.rec.addclp('cii')
            self.rec.addclp('nii')
            self.rec.addclp('oi')
            n_dets = 8
            for i in range(n_dets):
                self.rec.cii.padd(detparNames[i+1], detpar[:8,i], detparUnits[i])
                self.rec.nii.padd(detparNames[i+1], detpar[8:16,i], detparUnits[i])
                self.rec.oi.padd(detparNames[i+1], detpar[16:24,i], detparUnits[i])
#             self.nii = self.recBand(self, 0, detparNames, detpar[:8,:], detparUnits)
#             self.cii = self.recBand(self, 1, detparNames, detpar[8:16,:], detparUnits)
#             self.oi  = self.recBand(self, 2, detparNames, detpar[16:24,:], detparUnits)
            

            # =====================================================
            # read the telescope parameters with units
            pentries = 1
            tidtype = [('cband', 'U4', pentries), 
                       ('refDesignFreq', float, pentries),
                       ('refDesignWave', float, pentries),
                       ('primaryDiameter', float, pentries),
                       ('usePrimaryDiameter', float, pentries),
                       ('centralObscuration', float, pentries),
                       ('telescopeFNumber', float, pentries),
                       ('telescopeMagnification', float, pentries),
                       ('designEdgeTaper', float, pentries),
                       ('edgeTaperAperture', float, pentries),
                       ('apertureEfficiency', float, pentries),
                       ('mainLobeFWHM', float, pentries),
                       ('firstNull', float, pentries),
                       ('secondLobeHeight', float, pentries),
                       ('diffEncircledEnergy', float, pentries),
                       ('errorBeam', float, pentries),
                       ('focusSensitivity', float, pentries),
                       ('secondReflectionCoeff', float, pentries),
                       ('cassegrainBeamWaist', float, pentries),
                       ('cassegrainNumericalAperture', float, pentries)]
            # work-around: keeping units separate
            tUidtype = [('empty', 'U4'),
                        ('refDesignFreqUnit', 'a32'),
                        ('refDesignWaveUnit', 'a32'),
                        ('primaryDiameterUnit', 'a32'),
                        ('usePrimaryDiameterUnit', 'a32'),
                        ('centralObscurationUnit', 'a32'),
                        ('telescopeFNumberUnit', 'a32'),
                        ('telescopeMagnificationUnit', 'a32'),
                        ('designEdgeTaperUnit', 'a32'),
                        ('edgeTaperApertureUnit', 'a32'),
                        ('apertureEfficiencyUnit', 'a32'),
                        ('mainLobeFWHMUnit', 'a32'),
                        ('firstNullUnit', 'a32'),
                        ('secondLobeHeightUnit', 'a32'),
                        ('diffEncircledEnergyUnit', 'a32'),
                        ('errorBeamUnit', 'a32'),
                        ('focusSensitivityUnit', 'a32'),
                        ('secondReflectionCoeffUnit', 'a32'),
                        ('cassegrainBeamWaistUnit', 'a32'),
                        ('cassegrainNumericalApertureUnit', 'a32')]

            telpar = np.recarray([3,], dtype=tidtype)
            tnames = telpar.dtype.names[1:]
            telparunits = np.recarray([1,], dtype=tUidtype)
            tunames = telparunits.dtype.names[1:]
            tdf = wss['Telescope']
            sz = tdf.shape
            if verbose:
                print('Detectors parameters:')
                print('sz: ', sz) 
                print('columns: ', tdf.columns.values)
                print()
                print('Telescope parameters:')
            # loop over rows
            for i,tname,tuname in zip(range(0, sz[0]-2),tnames,tunames):
                trow = tdf.loc[i].values
                trsz = trow.size
                telpar[tname] = trow[2:trsz]
                telparunits[tuname] = trow[1]
                if verbose:
                    print(i, tname, tuname, telpar[tname], telparunits[tuname])
            self.telpar = telpar
            self.telparunits = telparunits

            # =====================================================
            # read the history sheet
            df = wss['History1']
            if verbose:
                print()
                # print('History worksheet: ', df.columns.values)
                print('History worksheet: ')
            sz = df.shape
            hist = np.recarray((sz[0],), dtype=[('Version', object), ('Date', object), ('By', 'U64'), ('Info','U256')] )
            if (sz[0] > 0) & (sz[1] > 0):
                for i in range(0, sz[0]):    # loop over the rows
                    if verbose:
                        print('sel by index: ', df.loc[i].values)
                    hist[i-1]['Version'] = df.loc[i].values[0]
                    hist[i-1]['Date'] = ('%s'%df.loc[i].values[1])[0:10]
                    hist[i-1]['By'] = df.loc[i].values[2]
                    hist[i-1]['Info'] = df.loc[i].values[3]
                    # print('sel by position: ', df.iloc[i,0])
            else:
                if verbose:
                    print('empty')
            self.hist = hist[:sz[0]-1]
            if verbose:
                print('hist2: \n', hist)
    
        
    def padd(self, pname, pvalue, punit='', pinfo='', level=0, levelnames=''):
        """Adding a parameter
        
        Input:
        ------
        
        """
        self.parlist.append(pname)
        setattr(self, pname, Par(pname, pvalue, punit, pinfo))
          
    
    def cadd(self, cname):
        """Adding an empty ReCclass
        
        Input:
        ------
        
        """
        self.parlist.append(cname)
        setattr(self, cname, RecClass())
        
        
    def getParList(self):
        #return self.parlist
        return self.__dict__

                
                
                
@dataclass
class Par(GL09PParameters):
    name: str
    value: float
    unit: str
    info: str



class RecClass:
    """Class to create nested parameter classes
    """

    def __init__(self): 
        pass

    def addclp(self, clname):
        """Add a parameter class with name clname
        """
        setattr(self,clname, Recpars())

    def padd(self, pname, pvalue, punit='', pinfo=''):
        """Add a parameter to the class
        """
        setattr(self, pname, Par(pname, pvalue, punit, pinfo))


class Recpars:

    def __init__(self): 
        self.par1 = 20

    def padd(self, pname, pvalue, punit='', pinfo=''):
        """Add a parameter to the class
        """
        setattr(self, pname, Par(pname, pvalue, punit, pinfo))

    def addclp(self, clname):
        """Add a parameter class with name clname
        """
        setattr(self,clname, Recpars())




if __name__ == '__main__':

    parfile = '/Users/volkertolls/git/gusto_misc/data/GUSTO_Baseline_Data_draft2.xlsx'
    gp = GL09PParameters(verbose=True, parfile=parfile)
    # gp = GL09PParametersDev(verbose=False)

    #print('\n\n\n')
    print()
    print('Parameter keys(): ')
    #pprint(gp.__dict__)
    print(gp.rec.cii.az.value)
    print()
    print('CII detector parameters:')
    pprint(gp.rec.cii.__dict__)

    print()
    print('Class History: ')
    print('Version: ', gp.__version__)
    print('Date created: ', gp.__date__)
    print('Date last updated', gp.__updated__)
    print('Class file: ', gp.__pyfile__)
    print('Parameter file loaded: ', gp.__parfile__)

    print()
    print('List all parameters: ')
    pprint(gp.getParList())
