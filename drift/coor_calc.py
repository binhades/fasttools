
# Filename: coor_calc.py
# Aim: to calculate the coordinates of the drift beams.

from numpy import zeros 
from numpy import array 
from ..time import utc2lst
from ..coor.coor_table import beam_coors, data2dict, csv_write # the last one is used for whold moudal 
from ..coor.coor_conv import radec2azel

def coor_beam_0(obs_utc, dec=0, hr=0):

    lst = utc2lst(obs_utc)
    ra = lst - hr
    ra_arr = ra.deg

    dec_arr = zeros(ra_arr.shape)
    dec_arr[:] = dec

    return array([ra_arr,dec_arr]) # shape: [2, leng]

def coor_beam_n(ra_c,dec_c,mjd,rot):
    beam = zeros((20,2,len(ra_c)))

    for i in range(len(ra_c)):
        beam[:,:,i] = beam_coors(beam_cx=ra_c[i],beam_cy=dec_c[i],theta=rot)

    ra = beam[:,0,:]
    dec = beam[:,1,:]
    az,el = radec2azel(ra,dec,mjd)
    coor_dict = data2dict(mjd,ra,dec,az,el)

    return coor_dict
