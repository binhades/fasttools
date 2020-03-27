
# Filename: beam.py
# Aim: to write the telescope beam to the fits data file.

def beam2fits(hdul,beam=1):
    if (not ('BEAMNUM' in hdul[1].header)) or hdul[1].header['BEAMNUM'] == -1:
        hdul[1].header.append(('BEAMNUM',beam))
        hdul.flush()
        return 1
    else:
        return 0

def beam_eff(beam=1,ZA=0,freq=None): # apply the frequency at 1400
    if beam == 1:
        a = -4.66e-4
        b = 6.57e-1
        c = -1.24e-2
    elif beam in [2,3,4,5,6,7]:
        a = 0.02e-4
        b = 6.31e-1
        c = -1.32e-2
    elif beam in [8,10,12,14,16,18]:
        a = 2.66e-4
        b = 5.80e-1
        c = -1.38e-2
    elif beam in [9,11,13,15,17,19]:
        a = 12.00e-4
        b = 5.53e-1
        c = -1.31e-2
    d = b+26.4*(a-c)

    if ZA < 26.4:
        eff=a*ZA+b
    else:
        eff=c*ZA+d

    return eff

def beam_gain(beam_eff,gain0=25.2): # K/Jy

    return beam_eff*gain0

def get_beam_from_filename(fname,type='raw'):
    if type == 'raw':
        ind = fname.find('_W_')
    elif type == 'drift':
        ind = fname.find('_M')
    beam = int(fname[ind+2:ind+4])
    return beam 
