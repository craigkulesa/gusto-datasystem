#!/usr/bin/env python
"""
This is the executable for the GUSTO L2 Pipeline.
"""

__date__ = '9/21/2024'
__updated__ = '20240921'
__version__ = '0.1'
__author__ = 'V. Tolls, CfA | Harvard & Smithsonian'

from joblib import Parallel, delayed
from joblib import Memory
import glob
import numpy as np
import time
import pkg_resources
import parsl
import sys
import os
import logging
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pathlib import Path
from pprint import pprint
import logging
log = logging.getLogger(__name__)


def GL095Pipeline(cfi, scanRange, verbose=False):
    """Function processing the Level 0.9 data and perform a
    baseline correction. 


    Parameters
    ----------
    scanRange : int

    Returns
    -------

    """
    # 1: read in the data from file
    # 2: apply bseline correction to each spectrum
    # 3: save back to files with only setting a flag about the baseline correction,
    # but all other parameters are passed through
    print('   GL095Pipeline to be implemented')