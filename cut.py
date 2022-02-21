
# Filename: cut.py
# Aim: to cut the fits data for ON/OFF observation.

import numpy as np
from astropy.io import fits
from copy import deepcopy

def cut_data2(data,ind_st,ind_end):

    row = ind_end - ind_st

    col1  = fits.Column(name='OBSNUM'  ,format='1K'      ,row=row,array=data['OBSNUM'][ind_st:ind_end])
    col2  = fits.Column(name='SCAN'    ,format='1K'      ,row=row,array=data['SCAN'    ][ind_st:ind_end])
    col3  = fits.Column(name='OBSTYPE' ,format='16A'     ,row=row,array=data['OBSTYPE' ][ind_st:ind_end])
    col4  = fits.Column(name='QUALITY' ,format='1L'      ,row=row,array=data['QUALITY' ][ind_st:ind_end])
    col5  = fits.Column(name='UTOBS'   ,format='1D'      ,row=row,array=data['UTOBS'   ][ind_st:ind_end])
    col6  = fits.Column(name='DATE-OBS',format='24A'     ,row=row,array=data['DATE-OBS'][ind_st:ind_end])
    col7  = fits.Column(name='OBJ-RA'  ,format='1D'      ,row=row,array=data['OBJ-RA'  ][ind_st:ind_end])
    col8  = fits.Column(name='OBJ-DEC' ,format='1D'      ,row=row,array=data['OBJ-DEC' ][ind_st:ind_end])
    col9  = fits.Column(name='OFF_RA'  ,format='1D'      ,row=row,array=data['OFF_RA'  ][ind_st:ind_end])
    col10 = fits.Column(name='OFF_DEC' ,format='1D'      ,row=row,array=data['OFF_DEC' ][ind_st:ind_end])
    col11 = fits.Column(name='TSYS'    ,format='1D'      ,row=row,array=data['TSYS'    ][ind_st:ind_end])
    col12 = fits.Column(name='EXPOSURE',format='1D'      ,row=row,array=data['EXPOSURE'][ind_st:ind_end])
    col13 = fits.Column(name='NCHAN'   ,format='1K'      ,row=row,array=data['NCHAN'   ][ind_st:ind_end])
    col14 = fits.Column(name='FREQ'    ,format='1D'      ,row=row,array=data['FREQ'    ][ind_st:ind_end])
    col15 = fits.Column(name='CHAN_BW' ,format='1D'      ,row=row,array=data['CHAN_BW' ][ind_st:ind_end])
    col16 = fits.Column(name='BEAM_EFF',format='1D'      ,row=row,array=data['BEAM_EFF'][ind_st:ind_end])
    col17 = fits.Column(name='PRESSURE',format='1D'      ,row=row,array=data['PRESSURE'][ind_st:ind_end])
    col18 = fits.Column(name='TAMBIENT',format='1D'      ,row=row,array=data['TAMBIENT'][ind_st:ind_end])
    col19 = fits.Column(name='WINDSPD' ,format='1D'      ,row=row,array=data['WINDSPD' ][ind_st:ind_end])
    col20 = fits.Column(name='WINDDIR' ,format='1D'      ,row=row,array=data['WINDDIR' ][ind_st:ind_end])
    col21 = fits.Column(name='DATA'    ,format='262144E' ,row=row,array=data['DATA'    ][ind_st:ind_end])
    col22 = fits.Column(name='CALSTAT ',format='1L'      ,row=row,array=data['CALSTAT '][ind_st:ind_end])
    cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15,col16,col17,col18,col19,col20,col21,col22])

    return cols

def cut_data(hdr,data,ind_st,ind_end):

    cols = []
    row = ind_end - ind_st

    for i in range(1, hdr['NAXIS2'] + 1):

        name = hdr['TTYPE{:d}'.format(i)]
        form = hdr['TFORM{:d}'.format(i)]

        if name == 'DATA':
            dim = hdr['TDIM{:d}'.format(i)]
            col = fits.Column(name=name,format=form,row=row,dim=dim,array=data[name][ind_st:ind_end])
        else:
            col = fits.Column(name=name,format=form,row=row,array=data[name][ind_st:ind_end])

        cols.append(col)

    return fits.ColDefs(cols)


def cutSigRef(hdul, on_st, on_end, off_st, off_end):

    prim_hdu = fits.PrimaryHDU(header=hdul[0].header)

    hdr_sig = deepcopy(hdul[1].header)
    hdr_ref = deepcopy(hdul[1].header)

    hdr_sig['NAXIS2'] = on_end - on_st
    hdr_ref['NAXIS2'] = off_end - off_st

    cols_on = cut_data(hdul.data,on_st,on_end)
    tab_hdu_on = fits.BinTableHDU.from_columns(cols_on)
    tab_hdu_on.header = hdr_sig
    hdul_on = fits.HDUList([prim_hdu,tab_hdu_on])
    
    cols_off = cut_data(hdul.data,off_st,off_end)
    tab_hdu_off = fits.BinTableHDU.from_columns(cols_off)
    tab_hdu_off.header = hdr_ref
    hdul_off = fits.HDUList([prim_hdu,tab_hdu_off])

    return hdul_on, hdul_off

