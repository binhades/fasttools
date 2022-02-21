
# Filename: time.py
# Aim: to convert utc to mjd or lst at location of the FAST site

def fast_time(obs_utc,scale='utc',format='isot'):
    from astropy.time import Time
    from .location import get_fast_location

    fast_loc = get_fast_location()
    t = Time(obs_utc,scale=scale,format=format,location=fast_loc)

    return t

def utc2mjd(obs_utc,scale='utc',format='isot'):

    t = fast_time(obs_utc,scale=scale,format=format)

    return t.mjd

def utc2lst(obs_utc,scale='utc',format='isot'):

    t = fast_time(obs_utc,scale=scale,format=format)
    lst = t.sidereal_time('mean')

    return lst

def mjd2fits(hdul,update=True):

    from astropy.io import fits

    obs_utc = hdul[1].data['DATE-OBS']
    obs_mjd = utc2mjd(obs_utc) + 0.5*hdul[1].data['EXPOSURE'][0]/(24.0*3600.0)

    if 'MJD' in hdul[1].columns.names:
        hdul[1].data['MJD'] = obs_mjd
    else:
        new_col = fits.Column(name='MJD',format='1D',array=obs_mjd)
        newhdu = fits.BinTableHDU.from_columns(hdul[1].columns + fits.ColDefs([new_col]))
        hdul[1].data = newhdu.data

    if update:
        hdul.flush()
    return obs_mjd
        
