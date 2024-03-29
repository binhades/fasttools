
# Filename: coor_table.py
# Aim: to calculate the telescope pointing from the XYZ file

import csv
import numpy as np
import datetime 
from astropy.time import Time
from .coor_conv import *
from ..obslog import rotlog

def beam_coors(beam_cx=0,beam_cy=0,theta=0,fod=0.5390185): # updated 20220914

    #beam seperate unit in degree. 270mm distance to anglar distance
    dist = np.degrees(np.arctan(0.270/(300*fod))) # f/d=0.4611, but scaled to match the measured beam seperation, then need 0.5390185

    beam = np.zeros((20,2))
    beam[0,0]  = beam_cx
    beam[0,1]  = beam_cy
    beam[1,0]  = beam[0,0] + 0
    beam[1,1]  = beam[0,1] + 0

    beam[2,1]  = beam[0,1] + dist*np.sin(theta+np.pi*0/3)
    beam[2,0]  = beam[0,0] + dist*np.cos(theta+np.pi*0/3)/np.cos(np.radians(beam[2,1]))
    beam[3,1]  = beam[0,1] + dist*np.sin(theta+np.pi*5/3)
    beam[3,0]  = beam[0,0] + dist*np.cos(theta+np.pi*5/3)/np.cos(np.radians(beam[3,1]))
    beam[4,1]  = beam[0,1] + dist*np.sin(theta+np.pi*4/3)
    beam[4,0]  = beam[0,0] + dist*np.cos(theta+np.pi*4/3)/np.cos(np.radians(beam[4,1]))
    beam[5,1]  = beam[0,1] + dist*np.sin(theta+np.pi*3/3)
    beam[5,0]  = beam[0,0] + dist*np.cos(theta+np.pi*3/3)/np.cos(np.radians(beam[5,1]))
    beam[6,1]  = beam[0,1] + dist*np.sin(theta+np.pi*2/3)
    beam[6,0]  = beam[0,0] + dist*np.cos(theta+np.pi*2/3)/np.cos(np.radians(beam[6,1]))
    beam[7,1]  = beam[0,1] + dist*np.sin(theta+np.pi*1/3)
    beam[7,0]  = beam[0,0] + dist*np.cos(theta+np.pi*1/3)/np.cos(np.radians(beam[7,1]))

    beam[8,1]  = beam[0,1] + 2*dist*np.sin(theta+np.pi*0/3)
    beam[8,0]  = beam[0,0] + 2*dist*np.cos(theta+np.pi*0/3)/np.cos(np.radians(beam[8,1]))
    beam[10,1] = beam[0,1] + 2*dist*np.sin(theta+np.pi*5/3)
    beam[10,0] = beam[0,0] + 2*dist*np.cos(theta+np.pi*5/3)/np.cos(np.radians(beam[10,1]))
    beam[12,1] = beam[0,1] + 2*dist*np.sin(theta+np.pi*4/3)
    beam[12,0] = beam[0,0] + 2*dist*np.cos(theta+np.pi*4/3)/np.cos(np.radians(beam[12,1]))
    beam[14,1] = beam[0,1] + 2*dist*np.sin(theta+np.pi*3/3)
    beam[14,0] = beam[0,0] + 2*dist*np.cos(theta+np.pi*3/3)/np.cos(np.radians(beam[14,1]))
    beam[16,1] = beam[0,1] + 2*dist*np.sin(theta+np.pi*2/3)
    beam[16,0] = beam[0,0] + 2*dist*np.cos(theta+np.pi*2/3)/np.cos(np.radians(beam[16,1]))
    beam[18,1] = beam[0,1] + 2*dist*np.sin(theta+np.pi*1/3)
    beam[18,0] = beam[0,0] + 2*dist*np.cos(theta+np.pi*1/3)/np.cos(np.radians(beam[18,1]))

    beam[9,1]  = beam[0,1] + 2*np.cos(np.pi*1/6)*dist*np.sin(theta+np.pi*11/6)
    beam[9,0]  = beam[0,0] + 2*np.cos(np.pi*1/6)*dist*np.cos(theta+np.pi*11/6)/np.cos(np.radians(beam[9,1]))
    beam[11,1] = beam[0,1] + 2*np.cos(np.pi*1/6)*dist*np.sin(theta+np.pi*9/6)
    beam[11,0] = beam[0,0] + 2*np.cos(np.pi*1/6)*dist*np.cos(theta+np.pi*9/6)/np.cos(np.radians(beam[11,1]))
    beam[13,1] = beam[0,1] + 2*np.cos(np.pi*1/6)*dist*np.sin(theta+np.pi*7/6)
    beam[13,0] = beam[0,0] + 2*np.cos(np.pi*1/6)*dist*np.cos(theta+np.pi*7/6)/np.cos(np.radians(beam[13,1]))
    beam[15,1] = beam[0,1] + 2*np.cos(np.pi*1/6)*dist*np.sin(theta+np.pi*5/6)
    beam[15,0] = beam[0,0] + 2*np.cos(np.pi*1/6)*dist*np.cos(theta+np.pi*5/6)/np.cos(np.radians(beam[15,1]))
    beam[17,1] = beam[0,1] + 2*np.cos(np.pi*1/6)*dist*np.sin(theta+np.pi*3/6)
    beam[17,0] = beam[0,0] + 2*np.cos(np.pi*1/6)*dist*np.cos(theta+np.pi*3/6)/np.cos(np.radians(beam[17,1]))
    beam[19,1] = beam[0,1] + 2*np.cos(np.pi*1/6)*dist*np.sin(theta+np.pi*1/6)
    beam[19,0] = beam[0,0] + 2*np.cos(np.pi*1/6)*dist*np.cos(theta+np.pi*1/6)/np.cos(np.radians(beam[19,1]))

    return beam

def csv_load(file_xyz,delimiter=',', delta=True):

    mjd=[]
    x=[]
    y=[]
    z=[]
    rot_SDP=[] # measure
    delta_rot = [] # measure - theory
    Yaw=[]
    Pitch=[]
    Roll=[] 

    with open(file_xyz,'rt') as filein:

        reader = csv.DictReader(filein,delimiter=delimiter)
        headers= reader.fieldnames
        for row in reader:
            x.append(float(row['SDP_PhaPos_X']))
            y.append(float(row['SDP_PhaPos_Y']))
            z.append(float(row['SDP_PhaPos_Z']))
            rot_SDP.append(float(row['SDP_AngleM']))
            Yaw.append(float(row['SDP_SwtDPose_Y']))
            Pitch.append(float(row['SDP_SwtDPose_P']))
            Roll.append(float(row['SDP_SwtDPose_R']))
            delta_rot.append(float(row['SDP_AngleM']) - float(row['TRP_AngleM']))
            ctime = row['SysTime']
            mjd.append(ctime2mjd(ctime))

    if delta:
        return np.array([mjd,x,y,z,delta_rot,Yaw,Pitch,Roll])
    else:
        return np.array([mjd,x,y,z,rot_SDP,Yaw,Pitch,Roll])


def csv_write(rows,fileout,delimiter=','):

    with open(fileout,'w',newline='') as csvfile:
        fieldnames = ['mjd','ra_m01','dec_m01','az_m01','el_m01',\
                            'ra_m02','dec_m02','az_m02','el_m02',\
                            'ra_m03','dec_m03','az_m03','el_m03',\
                            'ra_m04','dec_m04','az_m04','el_m04',\
                            'ra_m05','dec_m05','az_m05','el_m05',\
                            'ra_m06','dec_m06','az_m06','el_m06',\
                            'ra_m07','dec_m07','az_m07','el_m07',\
                            'ra_m08','dec_m08','az_m08','el_m08',\
                            'ra_m09','dec_m09','az_m09','el_m09',\
                            'ra_m10','dec_m10','az_m10','el_m10',\
                            'ra_m11','dec_m11','az_m11','el_m11',\
                            'ra_m12','dec_m12','az_m12','el_m12',\
                            'ra_m13','dec_m13','az_m13','el_m13',\
                            'ra_m14','dec_m14','az_m14','el_m14',\
                            'ra_m15','dec_m15','az_m15','el_m15',\
                            'ra_m16','dec_m16','az_m16','el_m16',\
                            'ra_m17','dec_m17','az_m17','el_m17',\
                            'ra_m18','dec_m18','az_m18','el_m18',\
                            'ra_m19','dec_m19','az_m19','el_m19']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)

    return 0


def ctime2mjd(ctime):

    t0 = datetime.datetime.strptime(ctime,'%Y-%m-%d %H:%M:%S.%f')
    dt = datetime.timedelta(hours=8)
    tz = datetime.timezone(dt)
    t1 = t0.replace(tzinfo=tz)
    tt = Time(t1,scale='utc')

    return tt.mjd

def data2dict(mjd,ra,dec,az,el):
    rows = []
    for i in range(len(mjd)):
        row = {'mjd':mjd[i],'ra_m01':ra[1,i],'dec_m01':dec[1,i],'az_m01':az[1,i],'el_m01':el[1,i],\
                            'ra_m02':ra[2,i],'dec_m02':dec[2,i],'az_m02':az[2,i],'el_m02':el[2,i],\
                            'ra_m03':ra[3,i],'dec_m03':dec[3,i],'az_m03':az[3,i],'el_m03':el[3,i],\
                            'ra_m04':ra[4,i],'dec_m04':dec[4,i],'az_m04':az[4,i],'el_m04':el[4,i],\
                            'ra_m05':ra[5,i],'dec_m05':dec[5,i],'az_m05':az[5,i],'el_m05':el[5,i],\
                            'ra_m06':ra[6,i],'dec_m06':dec[6,i],'az_m06':az[6,i],'el_m06':el[6,i],\
                            'ra_m07':ra[7,i],'dec_m07':dec[7,i],'az_m07':az[7,i],'el_m07':el[7,i],\
                            'ra_m08':ra[8,i],'dec_m08':dec[8,i],'az_m08':az[8,i],'el_m08':el[8,i],\
                            'ra_m09':ra[9,i],'dec_m09':dec[9,i],'az_m09':az[9,i],'el_m09':el[9,i],\
                            'ra_m10':ra[10,i],'dec_m10':dec[10,i],'az_m10':az[10,i],'el_m10':el[10,i],\
                            'ra_m11':ra[11,i],'dec_m11':dec[11,i],'az_m11':az[11,i],'el_m11':el[11,i],\
                            'ra_m12':ra[12,i],'dec_m12':dec[12,i],'az_m12':az[12,i],'el_m12':el[12,i],\
                            'ra_m13':ra[13,i],'dec_m13':dec[13,i],'az_m13':az[13,i],'el_m13':el[13,i],\
                            'ra_m14':ra[14,i],'dec_m14':dec[14,i],'az_m14':az[14,i],'el_m14':el[14,i],\
                            'ra_m15':ra[15,i],'dec_m15':dec[15,i],'az_m15':az[15,i],'el_m15':el[15,i],\
                            'ra_m16':ra[16,i],'dec_m16':dec[16,i],'az_m16':az[16,i],'el_m16':el[16,i],\
                            'ra_m17':ra[17,i],'dec_m17':dec[17,i],'az_m17':az[17,i],'el_m17':el[17,i],\
                            'ra_m18':ra[18,i],'dec_m18':dec[18,i],'az_m18':az[18,i],'el_m18':el[18,i],\
                            'ra_m19':ra[19,i],'dec_m19':dec[19,i],'az_m19':az[19,i],'el_m19':el[19,i]}
        rows.append(row)

    return rows

def multi_beam_offset(ra_c,dec_c,mjd,rot):

    beam = np.zeros((20,2,len(ra_c)))

    for i in range(len(ra_c)):
        beam[:,:,i] = beam_coors(beam_cx=ra_c[i],beam_cy=dec_c[i],theta=rot[i])

    ra = beam[:,0,:]
    dec = beam[:,1,:]
    az,el = radec2azel(ra,dec,mjd)
    
    coor_dict = data2dict(mjd,ra,dec,az,el)

    return coor_dict


def multi_beam_offset2(ra_c,dec_c,mjd,rot):
    # beam_coors returns the beam offset in the epoch of observing DATE, not J2000.
    # so we should transfer the ra_c, dec_c to observing epoch before get beam coords.
    # the error is about 0.5 ~ 5 arcsec << 1/20 beam

    beam = np.zeros((20,2,len(ra_c)))

    ra_c_obs, dec_c_obs = epoch_J2000_to_obs(ra_c, dec_c, mjd)
    for i in range(len(ra_c)):
        beam[:,:,i] = beam_coors(beam_cx=ra_c_obs[i],beam_cy=dec_c_obs[i],theta=rot[i])

    ra_obs = beam[:,0,:]
    dec_obs = beam[:,1,:]

    ra, dec = epoch_obs_to_J2000(ra_obs,dec_obs,mjd)
    az,el = radec2azel(ra,dec,mjd)
    
    coor_dict = data2dict(mjd,ra,dec,az,el)

    return coor_dict

def coor_table(file_xyz,fileout='coor_table.csv', delimiter=',', rot_log=None, rot0=0.0, cb=False, nowt=False):
    # we need set rot0 for RA and DE scans respectively.
    # also need to get the delta_rot = (mesure - theory) from csv_load, to exclude feild rotation of the receiver.

    if rot_log is not None:
        rot0 = rotlog(rot_log,file_xyz)
    if rot0 is None:
        rot0 = 0.0

    data = csv_load(file_xyz,delimiter=delimiter,delta=True)

    mjd = data[0,:]
    x   = data[1,:]
    y   = data[2,:]
    z   = data[3,:]
    ra_c,dec_c,az_c,el_c = coor_conv(mjd,x,y,z)
    rot = rot0*np.pi/180 + data[4,:] # delta = (measure - theory)

    if not cb:
        coor_dict = multi_beam_offset(ra_c,dec_c,mjd, rot)
        coor = coor_dict
        if not nowt:
            csv_write(coor_dict,fileout)
    else:
        coor = [mjd, ra_c, dec_c, az_c, el_c]

    return coor

def coor_table_multibeam(file_xyz, fileout='coor_table.csv', delimiter=',', nowt=False):
    # 2022.09.14
    # updated based on the code adapted from HiFAST

    data = csv_load(file_xyz,delimiter=delimiter,delta=False)

    mjd            = data[0,:]
    x              = data[1,:]
    y              = data[2,:]
    z              = data[3,:]
    multibeamAngle = data[4,:]
    Yaw            = data[5,:]
    Pitch          = data[6,:]
    Roll           = data[7,:]

    leng = mjd.shape[0]
    ra = np.zeros((20,leng))
    dec= np.zeros((20,leng))
    az = np.zeros((20,leng))
    el = np.zeros((20,leng))

    for i in range(19):
        beam_id = i+1
        ra[beam_id,:], dec[beam_id,:], az[beam_id,:], el[beam_id,:] = coor_conv_multibeam(mjd, beam_id, x, y, z, Yaw, Pitch, Roll, multibeamAngle)

    coor_dict = data2dict(mjd,ra,dec,az,el)
    if not nowt:
        csv_write(coor_dict,fileout)

    return coor_dict

def field_rotation(az, dec, el=None, mjd=None):

    fast_lat = 25.6534387
    rot = np.arcsin(np.sin(az*np.pi/180)/np.cos(dec*np.pi/180)*np.cos(fast_lat*np.pi/180))
    #----------------------------------------------------
    #ha, dec = azel2hadec(az,el,mjd) # not yet working
    #rot = np.arcsin(np.sin(ha*np.pi/180)/np.cos(el*np.pi/180)*np.cos(fast_lat*np.pi/180))
    #----------------------------------------------------
    return rot

