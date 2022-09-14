# Filename: coor_conv.py
# Aim: to convert the XYZ to RaDec

# Frame: HADec needs Astropy version 5.0+

from astropy.coordinates import SkyCoord, AltAz #, HADec
from astropy import units as u
from astropy.time import Time
from ..location import get_fast_location
import numpy as np
import math

#from astropy.utils import iers 
# to replace the default url to : http://maia.usno.navy.mil/ser7/finals2000A.all
#iers.conf.iers_auto_url =  'http://jlrat.bao.ac.cn/~bliu/doc/finals2000A.all'
    
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

def epoch_J2000_to_obs(ra,dec,mjd):

    obstime = Time(mjd,format='mjd')
    coor_icrs = SkyCoord(ra=ra*u.degree,dec=dec*u.degree,frame='icrs')
    coor_obs = coor_icrs.transform_to(FK5(equinox=obstime))
    ra_obs = coor_obs.ra.degree
    dec_obs= coor_obs.dec.degree

    return ra_obs,dec_obs

def epoch_obs_to_J2000(ra,dec,mjd):

    obstime = Time(mjd,format='mjd')
    coor_obs = SkyCoord(ra=ra*u.degree,dec=dec*u.degree,frame=FK5, equinox=obstime)
    coor_icrs = coor_obs.transform_to('icrs')
    ra_J2000 = coor_icrs.ra.degree
    dec_J2000= coor_icrs.dec.degree

    return ra_J2000,dec_J2000



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
   
def coor_conv_multibeam(mjd,beam_id,x0,y0,z0,Yaw,Pitch,Roll,multibeamAngle):

    beam_idx = np.arange(19)+1

    x,y,z = kypara2xyz_MultiBeam(beam_id,x0,y0,z0, Yaw, Pitch, Roll, multibeamAngle)         # units in degrees
    x = np.array(x).astype(np.float64) # must convert, otherwise ufunc error.
    y = np.array(y).astype(np.float64)
    z = np.array(z).astype(np.float64)
    az,el = xyz2azel(x,y,z) # units in degrees
    ra,dec = azel2radec(az,el,mjd)

    return ra,dec,az,el

# Following code: adapted from HiFAST, which adapted from FAST_ResultDataInversTool.exe 
    
def kypara2xyz(nB, globalCenterX, globalCenterY, globalCenterZ, globalYaw, globalPitch, globalRoll, multibeamAngle):
    """
    return Az: Azimuth in radian, El: Elevation
    """
    feedLocalCoord = [  [ 0.0  , 0.0  , 0.0   ],
                        [ -0.27, 0.0, 0.0],
                        [ -0.135, 0.233826859, 0.0],
                        [ 0.135, 0.233826859, 0.0],
                        [ 0.27, 0.0, 0.0],
                        [ 0.135, -0.233826859, 0.0],
                        [ -0.135, -0.233826859, 0.0],
                        [ -0.54, 0.0, 0.0],
                        [ -0.405, 0.233826859, 0.0],
                        [ -0.27, 0.467653718, 0.0],
                        [ 0.0, 0.467653718, 0.0],
                        [ 0.27, 0.467653718, 0.0],
                        [ 0.405, 0.233826859, 0.0],
                        [ 0.54, 0.0, 0.0],
                        [ 0.405, -0.233826859, 0.0],
                        [ 0.27, -0.467653718, 0.0],
                        [ 0.0, -0.467653718, 0.0],
                        [ -0.27, -0.467653718, 0.0],
                        [ -0.405, -0.233826859, 0.0] ]

    #1.多波束转角带来的旋转
    rotationMatrixMultiBeam= CalMultiBeamRotationMatrix(multibeamAngle)

    #2.Stewart下平台带来的旋转、
    rotationMatrixPlatform= CalPlatformRotationMatrix(globalYaw, globalPitch, globalRoll)

    #3.计算旋转后的位置
    useFeed = feedLocalCoord[nB-1]
    posRelative = VectorTransform(rotationMatrixPlatform, VectorTransform(rotationMatrixMultiBeam, useFeed))

    #4.计算全局坐标
    posAbsolute = posRelative + np.array([globalCenterX,  globalCenterY, globalCenterZ])

    return posAbsolute[0],posAbsolute[1],posAbsolute[2]

kypara2xyz_MultiBeam = np.frompyfunc(kypara2xyz, 8, 3)


#other functions
#从多波束转角计算旋转矩阵
def CalMultiBeamRotationMatrix(multibeamAngle):          
    """
    multibeamAngle: scalar  float64  Radian 
    """    
    ca= math.cos(multibeamAngle) 
    sa= math.sin(multibeamAngle) 
    rotationMatrixMultiBeam= np.zeros((3,3), dtype='float64')
    
    rotationMatrixMultiBeam[0,0] = ca 
    rotationMatrixMultiBeam[0,1] = -sa 
    rotationMatrixMultiBeam[0,2] = 0.0 
    rotationMatrixMultiBeam[1,0] = sa 
    rotationMatrixMultiBeam[1,1] = ca 
    rotationMatrixMultiBeam[1,2] = 0.0 
    rotationMatrixMultiBeam[2,0] = 0.0 
    rotationMatrixMultiBeam[2,1] = 0.0 
    rotationMatrixMultiBeam[2,2] = 1.0 
    return rotationMatrixMultiBeam

#从偏航、俯仰、翻滚计算旋转矩阵。
def CalPlatformRotationMatrix(globalYaw, globalPitch, globalRoll):
    """
    globalYaw, globalPitch, globalRoll: scalar  float64
    """
    
    cy = math.cos(globalYaw) 
    sy = math.sin(globalYaw) 
    cp = math.cos(globalPitch) 
    sp = math.sin(globalPitch) 
    cr = math.cos(globalRoll) 
    sr = math.sin(globalRoll) 
    
    rotationMatrixPlatform = np.zeros((3,3), dtype='float64')
    rotationMatrixPlatform[0,0] = cy * cp 
    rotationMatrixPlatform[0,1] = cy * sp * sr - sy * cr 
    rotationMatrixPlatform[0,2] = sy * sr + cy * sp * cr 
    rotationMatrixPlatform[1,0] = sy * cp 
    rotationMatrixPlatform[1,1] = cy * cr + sy * sp * sr 
    rotationMatrixPlatform[1,2] = sy * sp * cr - cy * sr 
    rotationMatrixPlatform[2,0] = -sp 
    rotationMatrixPlatform[2,1] = cp * sr 
    rotationMatrixPlatform[2,2] = cp * cr 
        
    return rotationMatrixPlatform

def VectorTransform(matrixA, vectorB):
    """
    """ 
    
    m11 = matrixA[0, 0] * vectorB[0] + matrixA[0, 1] * vectorB[1] + matrixA[0, 2] * vectorB[2]
    m21 = matrixA[1, 0] * vectorB[0] + matrixA[1, 1] * vectorB[1] + matrixA[1, 2] * vectorB[2]
    m31 = matrixA[2, 0] * vectorB[0] + matrixA[2, 1] * vectorB[1] + matrixA[2, 2] * vectorB[2]

   
    return  np.array([m11, m21, m31])
