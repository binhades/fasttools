
# Filename: coor2fits.py
# Aim: to write the telescope pointing to the fits data file.

from astropy.io import fits
from ..coor.coor2fits import csv_load

def radec2fits(hdul,tab):

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

def azel2fits(hdul,tab):
    
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
    radec2fits(hdul,tab)
    azel2fits(hdul,tab)

    hdul.flush()

    return 0
