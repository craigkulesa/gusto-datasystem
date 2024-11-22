"""Logging functionality for gusto L09 pipline
"""

import logging
import time
import os
import inspect


__version__ = 0.11
__date__ = '20240919'
__updated__ = '20240919'
fheader, tail = os.path.split(inspect.stack()[0][1])
__pyfile__ = tail
__path__ = fheader

def init_logger(loglevel='DEBUG', logfile=None, loggername='GL09PLogger'):
    # create logger 
    logger = logging.getLogger(loggername)
    logger.setLevel(loglevel)
    
    # create file handler which logs even debug messages
    if logfile is None:
        logfile = 'GL09P_%s.log'%(time.strftime("%Y%m%d%H%M%S"))
    fh = logging.FileHandler(logfile)
    fh.setLevel(loglevel)
    
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(loglevel)
    
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger
    
    
    
if __name__ == '__main__':
    
    class Auxiliary:
        def __init__(self):
            self.logger = logging.getLogger('spam_application.auxiliary.Auxiliary')
            self.logger.info('creating an instance of Auxiliary')
        def do_something(self):
            self.logger.info('doing something')
            a = 1 + 1
            self.logger.info('done doing something')
        
    def some_function():
        module_logger = logging.getLogger('GL09PLogger')
        module_logger.info('received a call to "some_function"')
    
    logger = init_logger(loglevel='DEBUG', logfile=None, loggername='GL09PLogger')
    
    
    logger.info('creating an instance of auxiliary_module.Auxiliary')
    a = Auxiliary()
    logger.info('created an instance of auxiliary_module.Auxiliary')
    logger.info('calling auxiliary_module.Auxiliary.do_something')
    a.do_something()
    logger.info('finished auxiliary_module.Auxiliary.do_something')
    logger.info('calling auxiliary_module.some_function()')
    some_function()
    logger.info('done with auxiliary_module.some_function()')