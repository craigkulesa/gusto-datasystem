"""
SPHEREx L4 Ices Simulator execution script.
"""

from .GL09Pipeline import runGL09P
from datetime import datetime


def execGL09P():
    print(f'{datetime.now().strftime("%d-%m-%Y %H:%M:%S")} Executing GUSTO Level 0.9 Pipeline.')

    runGL09P(verbose=True)

    print(f'{datetime.now().strftime("%d-%m-%Y %H:%M:%S")} Done running GUSTO Level 0.9 Pipeline.')
