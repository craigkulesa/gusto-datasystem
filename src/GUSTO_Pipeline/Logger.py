"""Logging functionality for gusto pipeline
"""

import logging
import time
import os


def init_logger(loglevel='DEBUG', logfile=None, loggername='pipelineLogger'):
    # create logger
    logger = logging.getLogger(loggername)
    logger.setLevel(loglevel)
    
    # create file handler which logs even debug messages
    if logfile is None:
        logfile = 'GUSTO_%s.log'%(time.strftime("%Y%m%d%H%M%S"))
    fh = logging.FileHandler(logfile)
    fh.setLevel(loglevel)

    # create console handler with a high log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)
    
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger
