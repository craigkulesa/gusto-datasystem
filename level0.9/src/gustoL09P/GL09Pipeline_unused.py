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


def GL09Pipeline(cfi, scanRange, verbose=False):
    """Function processing the Level 0.8 data. Input are uncalibrated 
    REF, HOT, and OTF spectra and output are calibrated OTF spectra


    Parameters
    ----------
    param1 : int

    Returns
    -------

    """
    print('   GL09Pipeline to be implemented')