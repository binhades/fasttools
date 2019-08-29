
# Filename: beam.py
# Aim: to write the telescope beam to the fits data file.

def beam2fits(hdul,beam=1):

    if (not ('BEAMNUM' in hdul[1].header)) or hdul[1].header['BEAMNUM'] == -1:

        hdul[1].header.append(('BEAMNUM',beam))
        hdul.flush()
        return 1

    else:
    
        return 0

