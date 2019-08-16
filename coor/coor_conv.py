#! /home/bliu/anaconda3/bin/python
# Filename: coor_conv.py
# Aim: to convert the XYZ to RaDec

import sys,os,glob,datetime,csv
import numpy as np
from astropy.coordinates import *
from astropy import units as u
from astropy.time import Time
    
def xyz2azel(x,y,z):
    # check if x,y,z are scalars or 1d array with same length.
    R = np.sqrt(x**2+y**2+z**2)
    az = np.degrees(np.arctan2(-x, -y))

    if az.size > 1:

        az[az<0] +=360

    else:

        if az < 0:
            az += 360

    el = np.degrees(np.arcsin(-z/R))

    return az, el

def radec2azel(ra,dec,mjd):

    fast_lon = 106.85645657571428*u.deg
    fast_lat = 25.6534387*u.deg
    fast_alt = 1138.72178*u.m
    fast_loc = EarthLocation(lat=fast_lat, lon=fast_lon, height=fast_alt)

    obswl = 21 * u.cm
    temperature = 25 * u.deg_C
    relative_humidity = 0.5
    pressure = 90000 * u.Pa
    obstime = Time(mjd,format='mjd')

    coor_icrs = SkyCoord(ra=ra*u.degree,dec=dec*u.degree,frame='icrs')
    coor_altaz = coor_icrs.transform_to(AltAz(obstime=obstime, obswl=obswl,\
                         relative_humidity=relative_humidity, pressure=pressure,\
                         temperature=temperature, location=fast_loc)) 
    alt = coor_altaz.alt.degree
    az  = coor_altaz.az.degree

    return az,alt

def azel2radec(az,el,mjd):

    fast_lon = 106.85645657571428*u.deg
    fast_lat = 25.6534387*u.deg
    fast_alt = 1138.72178*u.m
    fast_loc = EarthLocation(lat=fast_lat, lon=fast_lon, height=fast_alt)

    obswl = 21 * u.cm
    temperature = 25 * u.deg_C
    relative_humidity = 0.5
    pressure = 90000 * u.Pa
    obstime = Time(mjd,format='mjd')
    
    coor_azel = SkyCoord(az=az*u.degree, alt=el*u.degree, obstime=obstime,\
                         relative_humidity=relative_humidity, pressure=pressure,\
                         temperature=temperature, obswl=obswl, location=fast_loc, frame='altaz')
    coor_icrs = coor_azel.transform_to(ICRS) 
    ra = coor_icrs.ra.degree
    dec= coor_icrs.dec.degree

    return ra,dec

def coor_conv(mjd,x,y,z):

    az,el = xyz2azel(x,y,z)         # units in degrees
    ra,dec = azel2radec(az,el,mjd)

    return ra,dec,az,el

#----------------------------------
if __name__ == '__main__':
    if len(sys.argv)<4:

        print('input x,y,z')

    else:

        x = float(sys.argv[1])
        y = float(sys.argv[2])
        z = float(sys.argv[3])
        coor_conv(x,y,z)
