
# Filename: fits.py
# Aim: to re-write the FAST fits (ver:0.7) to the prefered fits data format.
# TODO: far from finished ...

import sys,csv
import numpy as np
from astropy.time import Time
from astropy.io import fits as pyfits

class LBFITS():
    def init(filein):

        try:
            ihdul = pyfits.open(filein)
        except FileNotFoundError:
            print("ERROR: LBFITS load FILE NOT FOUND - %S" % fname)
            return -1
        

        return 0
