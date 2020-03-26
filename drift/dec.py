#! /home/bliu/anaconda3/bin/python
# Filename: dec.py
# Aim: to get dec from the drift data path
# file_str = '/data/fast/drift/Dec+0500/20180919/Dec+0500_drifting-M01_W_0001.fits'

def get_dec(fname):

    try:
        decM = float(str.split(fname,'/')[-1][6:8])/60.
        decD = float(str.split(fname,'/')[-1][3:6])
    except ValueError:
        print("Error: dec get from path")
        return 0

    if decD > 0:
        dec = decD + decM
    else:
        dec = decD - decM

    return dec
