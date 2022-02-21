
# Filename: location.py
# Aim: to return the location of FAST site

from astropy import units as u
from astropy.coordinates import EarthLocation

def get_fast_location():

    fast_lon = 106.85645657571428 * u.deg
    fast_lat = 25.6534387 * u.deg
    fast_alt = 1138.72178 * u.m

    fast_loc = EarthLocation(lat=fast_lat, lon=fast_lon, height=fast_alt)

    return fast_loc
    
