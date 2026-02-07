"""
This is the GUSTO Pipeline run script.
"""
from .Logger import *
from .Configuration import *
from datetime import datetime
import importlib

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

    # initialize logging
    logDir = args.path + "log"
    os.makedirs(logDir, exist_ok=True)
    logfile = os.path.join(logDir,'pipeline_%s.log'%(time.strftime("%Y%m%d%H%M%S")))
    logger = init_logger(logfile=logfile, loggername='pipelineLogger')
    logger.info('Started logging to '+logfile)  # this will go to console, info messages only go to logfile

    # now figure out what levels we need to run
    levels = ['0.5', '0.7', '0.9', '1.0']
    funcNames = ['L05_Pipeline', 'L07_Pipeline', 'L09_Pipeline', 'L10_Pipeline']
    modules = ['.L05_lagstospectra', '.L07_telemetry', '.L09_calibrate', '.L10_pointing']
    
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
    
    logger.info(f'Executing pipeline levels: {exlevels}')
    scanRange = [int(x) for x in args.scanid]
    logger.info(f'Scan range for data processing: {scanRange}')
    print('#############################################################################\n')

    for index, step in enumerate(levels):
        if step in exlevels:
            logger.info(f'Starting GUSTO Pipeline Level {step}')
            stt = time.time()
            pipelineModule = importlib.import_module(modules[index], package='GUSTO_Pipeline')
            pipelineDef = getattr(pipelineModule, funcNames[index])
            n_files = pipelineDef(args, scanRange, verbose=verbose)
            ent = time.time()
            logger.info(f'Pipeline step done. {n_files} sequences or files were processed.')
            logger.info('Execution time: %.2fh  %.2fm   %.2fs'%((ent-stt)/3600.,(ent-stt)/60.,ent-stt))
            print('\n#############################################################################\n')
    print(f'{datetime.now().strftime("%d-%m-%Y %H:%M:%S")} Done running GUSTO Pipeline.')

    
if __name__ == '__main__':
    runGUSTO()
