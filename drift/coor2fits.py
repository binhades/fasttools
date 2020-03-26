
# Filename: coor2fits.py
# Aim: to write the telescope pointing to the fits data file.

import csv
from astropy.io import fits
from ..coor.coor2fits import csv_load

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


def radec2fits(hdul,beam=1,tab):

    tab_mjd = tab[0,:]
    tab_ra = tab[1,:]
    tab_dec = tab[2,:]
    if 'OBJ_RA' in hdul[1].columns.names:
        hdul[1].data['OBJ_RA']  = tab_ra
        hdul[1].data['OBJ_DEC'] = tab_dec
    else:
        col_ra = fits.Column(name='OBJ_RA',format='1D',array=tab_ra)
        col_dec= fits.Column(name='OBJ_DEC',format='1D',array=tab_dec)
        newcols = hdul[1].columns + fits.ColDefs([col_az,col_el])
        newhdu = fits.BinTableHDU.from_columns(newcols)
        hdul[1].data = newhdu.data

    return 0

def azel2fits(hdul,beam=1,tab):
    
    tab_mjd = tab[0,:]
    tab_az = tab[3,:]
    tab_el = tab[4,:]
   
    if 'AZ' in hdul[1].columns.names:
        hdul[1].data['AZ'] = tab_az
        hdul[1].data['EL'] = tab_el
    else:
        col_az = fits.Column(name='AZ',format='1D',array=tab_az)
        col_el = fits.Column(name='EL',format='1D',array=tab_el)
        newcols = hdul[1].columns + fits.ColDefs([col_az,col_el])
        newhdu = fits.BinTableHDU.from_columns(newcols)
        hdul[1].data = newhdu.data

    return 0

def coor2fits(hdul,beam=1,file_coor='coor_table.csv',delimiter=','):

    tab = csv_load(file_coor,delimiter=delimiter,beam=beam)
    radec2fits(hdul,tab):
    azel2fits(hdul,tab)

    #hdul.flush() # TODO

    return 0
