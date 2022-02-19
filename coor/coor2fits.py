
# Filename: coor2fits.py
# Aim: to write the telescope pointing to the fits data file.

import sys,csv
import numpy as np
from astropy.time import Time
from astropy.io import fits

def csv_load(filename,delimiter=',', beam=1):

    mjd=[]
    ra=[]
    dec=[]
    az=[]
    el=[]

    beamstr = '{:02d}'.format(beam)

    with open(filename,'rt') as filein:

        reader = csv.DictReader(filein,delimiter=delimiter)
        headers= reader.fieldnames
        for row in reader:
            mjd.append(float(row['mjd']))
            ra.append(float(row['ra_m'+beamstr]))
            dec.append(float(row['dec_m'+beamstr]))
            az.append(float(row['az_m'+beamstr]))
            el.append(float(row['el_m'+beamstr]))

    return np.array([mjd,ra,dec,az,el])

def utc2mjd(obs_utc):

    t = Time(obs_utc,scale='utc',format='isot')

    return t.mjd

def radec2fits(hdul,beam=1,file_coor='coor_table.csv',delimiter=','):

    tab = csv_load(file_coor,delimiter=delimiter,beam=beam)
    tab_mjd = tab[0,:]
    tab_ra = tab[1,:]
    tab_dec = tab[2,:]
   
    obs_utc = hdul[1].data['DATE-OBS'] + 0.5*hdul[1].data['EXPOSURE'][0]
    obs_mjd = utc2mjd(obs_utc)

    for i in range(len(obs_mjd)):
        ind = (np.abs(tab_mjd-obs_mjd[i])).argmin()
        hdul[1].data['OBJ_RA'][i] = tab_ra[ind]
        hdul[1].data['OBJ_DEC'][i] = tab_dec[ind]

    hdul.flush()

    return 0

def azel2fits(hdul,beam=1,file_coor='coor_table.csv',delimiter=','):
    
    tab = csv_load(file_coor,delimiter=delimiter,beam=beam)
    tab_mjd = tab[0,:]
    tab_az = tab[3,:]
    tab_el = tab[4,:]
   
    obs_utc = hdul[1].data['DATE-OBS'] + 0.5*hdul[1].data['EXPOSURE'][0]
    obs_mjd = utc2mjd(obs_utc) 

    obs_az = []
    obs_el = []

    for i in range(len(obs_mjd)):
        ind = (np.abs(tab_mjd-obs_mjd[i])).argmin()
        obs_az.append(tab_az[ind])
        obs_el.append(tab_el[ind])

    if 'AZ' in hdul[1].columns.names:
        hdul[1].data['AZ'] = obs_az
        hdul[1].data['EL'] = obs_el
    else:
        col_az = fits.Column(name='AZ',format='1D',array=obs_az)
        col_el = fits.Column(name='EL',format='1D',array=obs_el)
        newcols = hdul[1].columns + fits.ColDefs([col_az,col_el])
        newhdu = fits.BinTableHDU.from_columns(newcols)
        hdul[1].data = newhdu.data

    hdul.flush()

    return 0

