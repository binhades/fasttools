
# Filename: beam.py
# Aim: to write the telescope beam to the fits data file.

import numpy as np
import csv

def get_beam_from_filename(fname,type='raw'):
    if type == 'raw':
        ind = fname.find('_W_')
    elif type == 'drift':
        ind = fname.find('_M')
    beam = int(fname[ind+2:ind+4])
    return beam 

def beam2fits(hdul,beam=1):
    if (not ('BEAMNUM' in hdul[1].header)):
        hdul[1].header.append(('BEAMNUM',beam))
    else:
        hdul[1].header['BEAMNUM'] = beam

    hdul.flush()
    return 1

def unit2fits(hdul,unit='K',comment='K'):
    if (not ('BUNIT' in hdul[1].header)):
        hdul[1].header.append(('BUN299792458/(frequency*1e6IT',unit,comment))
    else:
        hdul[1].header['BEAMNUM'] = (beam,comment)

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
    with open('./beam_size.csv','rt') as filein:
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
    bid = int(beam)
    if bid < 1 or bid > 19:
        print("BeamID {:d} is not correct".format(bid))
        return None
    if frequency < 1000 or frequency > 1500:
        print("Frequency {:f} is out of range".format(frequency))
        return None
    with open('./aper_effi_para.csv','rt') as filein:
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

    if ZA < 26.4:
        eff=a0*ZA+b0
    else:
        eff=c0*ZA+d0

    if plot:
        import matplotlib.pyplot as plt
        def get_ne(za, a0, b0, c0, d0):
            if za < 26.4:
                return a0*za+b0
            else:
                return c0*za+d0

        za_arr = np.arange(0,40,0.1)
        ne_arr = []
        for za in za_arr:
            ne_arr.append(get_ne(za,a0,b0,c0,d0))
        ne_arr = np.array(ne_arr)
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


def Sv_to_Tmb(Sv, beam=1, frequency=frequency):
    lam = 299792458./(1e6*frequency)
    theta = beam_size(beam=beam,frequency=frequency)/60. * np.pi/180.
    omega = np.pi/4.*np.log(2) * (theta)**2
    Tmb = lam**2 * Sv / (2*1.38e-23*omega)

    return Tmb

def Tmb_to_Sv(Tmb, beam=1, frequency=frequency):
    lam = 299792458./(1e6*frequency)
    theta = beam_size(beam=beam,frequency=frequency)/60. * np.pi/180.
    omega = np.pi/4.*np.log(2) * (theta)**2
    Sv = Tmb * (2*1.38e-23*omega)/lam**2

    return Sv

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

def antenna_gain(beam=1, ZA=0, frequency=1400.): # K/Jy
    gain0=25.6
    ape_eff = aperture_eff(beam=beam,ZA=ZA,frequency=frequency)

    return ape_eff*gain0


