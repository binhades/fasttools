# Filename: coor_conv.py
# Aim: to convert the XYZ to RaDec

# Frame: HADec needs Astropy version 5.0+

from astropy.coordinates import SkyCoord, AltAz #, HADec
from astropy import units as u
from astropy.time import Time
from ..location import get_fast_location

#from astropy.utils import iers 
# to replace the default url to : http://maia.usno.navy.mil/ser7/finals2000A.all
#iers.conf.iers_auto_url =  'http://jlrat.bao.ac.cn/~bliu/doc/finals2000A.all'
    
def xyz2azel(x,y,z):
    import numpy as np
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

    fast_loc = get_fast_location()

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

    fast_loc = get_fast_location()

    obswl = 21 * u.cm
    temperature = 25 * u.deg_C
    relative_humidity = 0.5
    pressure = 90000 * u.Pa
    obstime = Time(mjd,format='mjd')
    
    coor_azel = SkyCoord(az=az*u.degree, alt=el*u.degree, obstime=obstime,\
                         relative_humidity=relative_humidity, pressure=pressure,\
                         temperature=temperature, obswl=obswl, location=fast_loc, frame='altaz')
    coor_icrs = coor_azel.transform_to('icrs') 
    ra = coor_icrs.ra.degree
    dec= coor_icrs.dec.degree

    return ra,dec

def coor_conv(mjd,x,y,z):

    az,el = xyz2azel(x,y,z)         # units in degrees
    ra,dec = azel2radec(az,el,mjd)

    return ra,dec,az,el

# ---------------------------------
# Frame: HADec needs Astropy version 5.0+
#def radec2hadec(ra,dec,mjd):
#
#    fast_loc = get_fast_location()
#
#    obswl = 21 * u.cm
#    temperature = 25 * u.deg_C
#    relative_humidity = 0.5
#    pressure = 90000 * u.Pa
#    obstime = Time(mjd,format='mjd')
#
#    coor_icrs = SkyCoord(ra=ra*u.degree,dec=dec*u.degree,frame='icrs')
#    coor_hadec = coor_icrs.transform_to(HADec(obstime=obstime, obswl=obswl,\
#                         relative_humidity=relative_humidity, pressure=pressure,\
#                         temperature=temperature, location=fast_loc)) 
#    ha  = coor_hadec.ha.degree
#    dec = coor_hadec.dec.degree
#
#    return ha, dec
#
#def azel2hadec(az,el,mjd):
#
#    fast_loc = get_fast_location()
#
#    obswl = 21 * u.cm
#    temperature = 25 * u.deg_C
#    relative_humidity = 0.5
#    pressure = 90000 * u.Pa
#    obstime = Time(mjd,format='mjd')
#    
#    coor_azel = SkyCoord(az=az*u.degree, alt=el*u.degree, obstime=obstime,\
#                         relative_humidity=relative_humidity, pressure=pressure,\
#                         temperature=temperature, obswl=obswl, location=fast_loc, frame='altaz')
#    coor_hadec = coor_azel.transform_to(HADec(obstime=obstime, obswl=obswl,\
#                         relative_humidity=relative_humidity, pressure=pressure,\
#                         temperature=temperature, location=fast_loc)) 
#    ha = coor_hadec.ha.degree
#    dec= coor_hadec.dec.degree
#
#    return ha,dec
# ---------------------------------
