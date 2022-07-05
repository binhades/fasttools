
# Filename: beam.py
# Aim: to write the telescope beam to the fits data file.

import numpy as np
import os, csv

def get_beam_from_filename(fname,type='raw'):
    if type == 'raw':
        ind = fname.find('_W_')
    elif type == 'drift':
        ind = fname.find('_M')
    beam = int(fname[ind+2:ind+4])
    return beam 

def beam2fits(hdul,beam=1):

    hdul[1].header.set('BEAMID',beam,'Beam number of FAST 19-beams',before='EXTVER')
    hdul.flush()
    return 1

def unit2sdfits(hdul,unit='K',comment=None):
    ind = hdul[1].columns.names.index('DATA')
    hdul[1].header.set(f'TUNIT{ind+1}', unit, comment)
    hdul.flush()

def unit2fits(hdul,unit='K',comment=None):
    hdul[1].header.set('CUNIT2', unit, comment)
    hdul.flush()
    return 1

def beam_size(beam=1,frequency=1400.,plot=False):

    bid = int(beam)
    if bid < 1 or bid > 19:
        print("BeamID {:d} is not correct".format(bid))
        return None
    if frequency < 1000 or frequency > 1500:
        print("Frequency {:f} is out of range".format(frequency))
        return None
    file_beam = os.path.join(os.path.dirname(__file__),'beam_size.csv')
    with open(file_beam,'rt') as filein:
        reader = csv.reader(filein, delimiter=',')
        header = reader.__next__()
        freq_arr = np.array(header[1:]).astype(np.float)
        lambda_arr = 299792458/(1.e6*freq_arr)
        for row in reader:
            if row[0] == 'M{:02d}'.format(bid):
                beam_size = np.array(row[1:]).astype(np.float)

    bc = np.polyfit(lambda_arr,beam_size,3) 
    bf = np.poly1d(bc)
    lambda0 = 299792458/(frequency*1e6)

    if plot:
        import matplotlib.pyplot as plt
        lambda1 = np.arange(lambda_arr[-1],lambda_arr[0],0.001)
        plt.plot(lambda_arr*100, beam_size,'o')
        plt.plot(lambda1*100, bf(lambda1))
        plt.xlabel('$\lambda$ (cm)')
        plt.ylabel('arcmin')
        plt.title('Beam: M{:02d}'.format(bid))
        plt.show()

    return bf(lambda0)

def aperture_eff(beam=1, ZA=0, frequency=1400.,plot=False):

    def eff4za(za, a, b, c, d):
        return np.select([za<26.4, za>=26.4],[a*za+b,c*za+d])
    eff4za_ufunc = np.frompyfunc(eff4za,5,1)

    bid = int(beam)
    if bid < 1 or bid > 19:
        print("BeamID {:d} is not correct".format(bid))
        return None
    if frequency < 1000 or frequency > 1500:
        print("Frequency {:f} is out of range".format(frequency))
        return None
    file_aper = os.path.join(os.path.dirname(__file__),'aper_effi_para.csv')
    with open(file_aper,'rt') as filein:
        reader = csv.reader(filein, delimiter=',')
        header = reader.__next__()
        freq_arr = np.array(header[1:]).astype(np.float)
        for row in reader:
            if row[0] == 'M{:02d}a'.format(bid):
                c_a = np.array(row[1:]).astype(np.float)
            if row[0] == 'M{:02d}b'.format(bid):
                c_b = np.array(row[1:]).astype(np.float)
            if row[0] == 'M{:02d}c'.format(bid):
                c_c = np.array(row[1:]).astype(np.float)
            if row[0] == 'M{:02d}ae'.format(bid):
                e_a = np.array(row[1:]).astype(np.float)
            if row[0] == 'M{:02d}be'.format(bid):
                e_b = np.array(row[1:]).astype(np.float)
            if row[0] == 'M{:02d}ce'.format(bid):
                e_c = np.array(row[1:]).astype(np.float)
 
    pca = np.polynomial.chebyshev.chebfit(freq_arr,c_a,5,w=1/e_a) 
    pcb = np.polynomial.chebyshev.chebfit(freq_arr,c_b,5,w=1/e_b) 
    pcc = np.polynomial.chebyshev.chebfit(freq_arr,c_c,5,w=1/e_c) 

    a0 = np.polynomial.chebyshev.chebval(frequency, pca) * 1e-4
    b0 = np.polynomial.chebyshev.chebval(frequency, pcb) * 1e-1
    c0 = np.polynomial.chebyshev.chebval(frequency, pcc) * 1e-2
    d0 = b0+26.4*(a0-c0)

    eff = eff4za_ufunc(ZA,a0,b0,c0,d0)
    eff = eff.astype(np.float)

#    if ZA < 26.4:
#        eff=a0*ZA+b0
#    else:
#        eff=c0*ZA+d0

    if plot:
        import matplotlib.pyplot as plt

        za_arr = np.arange(0,40,0.1)
        ne_arr = eff4za_ufunc(za_arr,a0,b0,c0,d0)
        ne_arr.astype(np.float)
        plt.plot(za_arr, ne_arr)
        plt.show()

        #a1 = np.polynomial.chebyshev.chebval(freq_arr, pca) 
        #plt.errorbar(freq_arr, c_a, yerr=e_a,linestyle='--')
        #plt.plot(freq_arr, a1)
        #plt.show()

    return eff
# ------------------------------------------------------------------------
# OUTDATED, NOT USE 
#def aperture_eff(beam=1,ZA=0,freq=None): # apply the frequency at 1400
#    # Antenna Aperture Efficiency
#    if beam == 1:
#        a = -4.66e-4
#        b = 6.57e-1
#        c = -1.24e-2
#    elif beam in [2,3,4,5,6,7]:
#        a = 0.02e-4
#        b = 6.31e-1
#        c = -1.32e-2
#    elif beam in [8,10,12,14,16,18]:
#        a = 2.66e-4
#        b = 5.80e-1
#        c = -1.38e-2
#    elif beam in [9,11,13,15,17,19]:
#        a = 12.00e-4
#        b = 5.53e-1
#        c = -1.31e-2
#    d = b+26.4*(a-c)
#
#    if ZA < 26.4:
#        eff=a*ZA+b
#    else:
#        eff=c*ZA+d
#
#    return eff
# ------------------------------------------------------------------------


def Sv_to_Tmb(Sv, beam=1, frequency=1400.,theta=None):
    lam = 299792458./(1e6*frequency)
    if theta is None:
        theta = beam_size(beam=beam,frequency=frequency)/60. * np.pi/180.
    else:
        theta = theta/60. * np.pi/180
    omega = np.pi/(4.*np.log(2)) * (theta)**2
    r_tmb_sv = lam**2 / (2*1.38e-23*omega) * 1e-26

    return r_tmb_sv * Sv

def Tmb_to_Sv(Tmb, beam=1, frequency=1400.,theta=None):
    lam = 299792458./(1e6*frequency)
    if theta is None:
        theta = beam_size(beam=beam,frequency=frequency)/60. * np.pi/180.
    else:
        theta = theta/60. * np.pi/180
    omega = np.pi/(4.*np.log(2)) * (theta)**2
    r_sv_tmb = (2*1.38e-23*omega)/lam**2 * 1e26

    return r_sv_tmb * Tmb

def main_beam_eff(beam=1, ZA=0, frequency=1400.): # apply the frequency at 1400
    """ Antenna Main-beam Efficiency
        I adopt the equation given by Ronald J Maddalena to calculate 
        MB Efficiency from Aperture Efficiency
        https://www.gb.nrao.edu/~rmaddale/Research/TheoreticalBeamEfficiency.pdf
     """
    D = 300 # m
    n_R = 1.0
    lam = 299792458./(1e6*frequency)
    theta = beam_size(beam=beam,frequency=frequency)/60. * np.pi/180.
    ape_eff = aperture_eff(beam=beam, ZA=ZA, frequency=frequency)  
    mb_eff = 0.8899 * ape_eff / n_R * theta**2 * D**2 / lam**2

    return mb_eff

def FAST_Gain0(): # K/Jy
    #gain0 = 1e-26/(2*1.38e-23) * np.pi*(150.0**2) 
    gain0 = 25.61
    return 25.61

def antenna_gain(beam=1, ZA=0, frequency=1400.): # K/Jy
    gain0 = FAST_Gain0()
    ape_eff = aperture_eff(beam=beam,ZA=ZA,frequency=frequency)
    return ape_eff*gain0

def beam_coors(beam_cx=0,beam_cy=0,theta=0,fod=0.5365,id=None):
    # copied from .coor/coor_table.py
    #beam seperate unit in arcmin. 270mm distance to anglar distance
    #-------------------------
    # add 20210331
    # As far as I know, the horn-dist is 270mm and the fod is 0.4621.
    # This give the beam distance of 6.695 arcmin, not matching with Jiang2020.
    # According to Jiang2020, we set beam dist to 5.76, from (270mm && 0.5356fod)
    # But now I see the updated value of dist = 5.74 from the FAST offical train workshop (20191113)
    # this need: (270mm && 0.5390fod) or (231.47mm && 0.4621fod)
    #-------------------------
    # add 20220705
    # If the RA axis is to the left (9876543210), theta should be positve for clockwise rotation.
    # If the RA axis is to the right (012345678), theta should be positve for counterclockwise rotation.
    # I have checked the following code is corrected for this sake.
    #-------------------------

    dist = np.degrees(np.arctan(0.270/(300*fod))) # f/d=0.4611, but scaled to match the measured beam seperation, then need 0.5365
    theta = np.deg2rad(theta)

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

    if id is None:
        return beam
    else:
        return beam[int(id),:]
