Debug mode of Level 0.9 Pipeline
================================

The debug mode of the Level 0.9 pipeline is invoked using either the command line switches -d or --debug, 
or by setting "debug = True" in the configuration file. 

In debug mode, the pipeline creates additional data products to enable easier debugging. The following data 
products will be available in their respective extensions. The data are in tables and the list below shows 
all the table columns with their respective numpy dtype format.

Extension 1: data1:

	('MIXER', '>i4')
	('NINT', '>i4')
	('UNIXTIME', '>f8')
	('NBYTES', '>i4')
	('CORRTIME', '>i4')
	('INTTIME', '>f4')
	('ROW_FLAG', '>i4')
	('Ihigh', '>i4')
	('Qhigh', '>i4')
	('Ilow', '>i4')
	('Qlow', '>i4')
	('Ierr', '>i4')
	('Qerr', '>i4')
	('VIhi', '>f4')
	('VQhi', '>f4')
	('VIlo', '>f4')
	('VQlo', '>f4')
	('scanID', '>i4')
	('subScan', '>i4')
	('scan_type', 'S6')
	('THOT', '>f4')
	('RA', '>f4')
	('DEC', '>f4')
	('filename', 'S48')
	('PSat', '>f4')
	('Imon', '>f4')
	('Gmon', '>f4')
	('spec', '>f4', (512,))
	('CHANNEL_FLAG', '>i2', (512,))


Extension 2: data2:

	('Tsyseff', '>f8', (512,))
	('hcorr', '>f8', (512,))
	('spref', '>f8', (512,))
	('spref2', '>f8', (512,))
	('frac', '>f8', (2,))
	('spec', '>f4', (512,))
	('tant', '>f8', (512,))
	('ringvar', '>f8')])


Extension 3: data3:

	('hots', '>f8', (512,))
	('htime', '>f8')
	('hmixer', '>f8')
	('htint', '>f8')


Extension 4: data4

	('Tsys', '>f4', (2, 512))
	('tsmix', '>i4')


Extension 5: data5:

	('Ta', '>f8', (512,)

