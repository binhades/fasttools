
# Filename: cal.py
# Aim: to write the telescope cal to the fits data file.

import numpy as np
from astropy.io import fits

def getcal(calon,caloff,expose,leng):

    len_on =  int(calon/expose)
    len_off =  int(caloff/expose)
    period = int(len_on+len_off)

    N = int(leng/period)+1
    cal_p = np.zeros(period,dtype=bool)
    cal_p[0:len_on] = True
    calarr = np.tile(cal_p,N)
    calstat = calarr[0:leng]

    return calstat

def cal2fits(hdul,calon=1,caloff=1,expose=1,caltype='high'):

    leng = hdul[1].header['NAXIS2']
    calstat = getcal(calon,caloff,expose,leng)
    hdul[1].header['CALTYPE'] = caltype
    if 'CALSTAT' in hdul[1].columns.names:
        hdul[1].data['CALSTAT'] = calstat
    else:
        new_col = fits.Column(name='CALSTAT',format='1L',array=calstat)
        newhdu = fits.BinTableHDU.from_columns(hdul[1].columns + fits.ColDefs([new_col]))
        hdul[1].data = newhdu.data
    hdul.flush()

    return 1

