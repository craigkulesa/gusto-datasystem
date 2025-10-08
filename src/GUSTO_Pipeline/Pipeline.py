"""
This is the GUSTO Pipeline run script.
"""
from .Logger import *
from .L05_lagstospectra import *
from .L07_telemetry import *
from .L08_chanflags import *
from .L09_calibrate import *
from .L10_pointing import *
from .Configuration import *
from datetime import datetime

def runGUSTO(verbose=False):
    r"""Function executing the GUSTO pipeline.
    1) read the command line and configuration parameters
    2) initialize the pipeline
       - analyze which parts of the pipeline are executed
    3) execute the selected pipeline parts
       - calling the appropriate sub-pipeline functions
    """
    print(f'{datetime.now().strftime("%d-%m-%Y %H:%M:%S")} Executing GUSTO Pipeline.')
    if verbose: 
        print('\n%s: Initializing GUSTO Pipeline'%(time.strftime("%c")))

    args = getConfiguration(verbose=verbose)
    verbose = args.verbose
    
    # initialize logging:
    logDir = args.path + "log"
    os.makedirs(logDir, exist_ok=True)
    
    loglevel = 'INFO'
    logfile = os.path.join(logDir,'pipeline_%s.log'%(time.strftime("%Y%m%d%H%M%S")))
    logger = init_logger(loglevel=loglevel, logfile=logfile, loggername='pipelineLogger')
    logger.warning('Started logging to '+logfile)  # this will go to console, info messages only go to logfile

    levels = ['0.5', '0.7', '0.8', '0.9', '1.0']
    isl = args.level[0]
    if isl not in levels:
        print('Please check start data level ...')
    iarg = levels.index(isl)
    if len(args.level) > 1:
        iel = args.level[1]
    else:
        iel = isl
    if iel not in levels:
        print('Please check end data level!\nPossible levels are: %s'%(levels[iarg:]))
    earg = levels.index(iel)
    exlevels = levels[iarg:earg+1]
    
    print('\nExecuting pipeline levels: ', exlevels)
    scanRange = [int(x) for x in args.scanid]
    print('Scan range for data processing: ', scanRange)
    print('#############################################################################')
    
    # this is the part that forks to the various pipeline levels
    if '0.5' in exlevels:
        print('Executing Level 0.5: generating power spectra from lags and saving as SDFITS')
        n_files = L05_Pipeline(args, scanRange)
        print('Level 0.5 to 0.7 done. ', n_files, ' lag files were processed.\n')        
        print('#############################################################################')
        
    if '0.7' in exlevels:
        print('Executing Level 0.7: generating telemetry headers and making sequence files')
        n_files = L07_Pipeline(args)  # it will generate scanRange on its own
        print('Level 0.7 to 0.8 done. ', n_files, ' sequences were processed.\n')
        print('#############################################################################')
        
    if '0.8' in exlevels:
        print('Executing Level 0.8: channel flags for spurs and row flags for bad fringing')
        n_files = L08_Pipeline(args, scanRange)
        print('Level 0.8 to 0.9 done. ', n_files, ' SDFITS files were processed.\n')
        print('#############################################################################')
        
    if '0.9' in exlevels:
        stt = time.time()
        print('\n%s: Executing 0.9: calibrating spectra'%(time.strftime("%c")))
        n_files = L09_Pipeline(args, scanRange, verbose=verbose)
        ent = time.time()
        print('Level 0.9 pipeline done.\n')
        print('Execution time: %.2fh  %.2fm   %.2fs'%((ent-stt)/3600.,(ent-stt)/60.,ent-stt))
        print('Processed %i files.'%(n_files))
        print()
        print('#############################################################################')
    
    if '1.0' in exlevels:
        print('Executing Level 1.0 pipeline: coordinate corrections')
        n_files = L10_Pipeline(args, scanRange, verbose=verbose)
        print('Level 0.9 to 1.0 done. ', n_files, 'SDFITS files were processed.\n')
        print('#############################################################################')
    print()
    print(f'{datetime.now().strftime("%d-%m-%Y %H:%M:%S")} Done running GUSTO Pipeline.')
    

if __name__ == '__main__':
    runGUSTO()
