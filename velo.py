
# Filename: time.py
# Aim: to convert utc to mjd or lst at location of the FAST site

import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.io import fits
#from astropy.coordinates import *
from astropy.coordinates import SkyCoord
from .location import get_fast_location
# Schonrich et al.(2010) Barycentric velocity relative to the LSR, which is defined in
# Galactic (right-handed) cartesian velocity components (U,V,W) = (11.1,12.24,7.25) km/s

# v_lsr = v_bary+U*cos(l)cos(b)+V*sin(l)*cos(b)+W*sin(b) # l,b are the Galactic coordinates.

def getVframeLSR(ra,dec,mjd):
    fast_loc = get_fast_location()
    obstime = Time(mjd,format='mjd')
    sc = SkyCoord(ra=ra*u.degree,dec=dec*u.degree,frame='icrs')
    b=sc.galactic.b.rad
    l=sc.galactic.l.rad
    barycorr = sc.radial_velocity_correction(obstime=obstime,location=fast_loc)
    vlsr=barycorr.value/1000.0+11.1*np.cos(l)*np.cos(b)+12.24*np.sin(l)*np.cos(b)+7.25*np.sin(b)
    return vlsr # in km/s

def vframe2fits(hdul, ra, dec, mjd, update=True):

    vframe = getVframeLSR(ra,dec,mjd)

    if 'VFRAME' in hdul[1].columns.names:
        hdul[1].data['VFRAME'] = vframe
    else:
        new_col = fits.Column(name='VFRAME',format='1D',array=vframe)
        newhdu = fits.BinTableHDU.from_columns(hdul[1].columns + fits.ColDefs([new_col]))
        hdul[1].data = newhdu.data

    if update:
        hdul.flush()

    return vframe

def freq2velo(freq,line_rest_freq): # freq in MHz

    rest_freq = line_rest_freq * u.MHz
    frequency = freq * u.MHz
    relativistic_equiv = u.doppler_relativistic(rest_freq)

    velo = frequency.to(u.km/u.s,equivalencies=relativistic_equiv)

    return velo.value # km/s

def velo2freq(velo,line_rest_freq): #velo in km/s;  freq in MHz

    rest_freq = line_rest_freq * u.MHz
    velocity = velo * u.km/u.s
    relativistic_equiv = u.doppler_relativistic(rest_freq)

    freq = velocity.to(u.MHz,equivalencies=relativistic_equiv)

    return freq.value # freq in MHz
