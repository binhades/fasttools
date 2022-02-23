# -*-coding: utf-8 -*-
# Filename: obslog.py
# Aim: to read information from observation logs.

import csv

def rotlog(file_log, file_obs, delimiter=','):
    ''' Rotation Log
    '''
    flist=[]
    rot=[]

    with open(file_log,'rt') as filein:

        reader = csv.DictReader(filein,delimiter=delimiter)
        headers= reader.fieldnames
        for row in reader:
            rot.append(float(row['rot']))
            flist.append(row['filename'].strip())
    try:
        ind = flist.index(file_obs)
        return rot[ind]
    except ValueError:
        return None

