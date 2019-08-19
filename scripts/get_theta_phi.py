import numpy as np

def raw2det(detector='mos1', rawx=1,rawy=1,ccdnr=1):
    #
    # PURPOSE:
    #   direct conversion of RAWX, RAWY to DETX,DETY
    # INPUTS:
    #   detector - the XMM detector, can be 'mos1', 'mos2' or 'pn'
    #           NOTE: only implemented for mos1 at the moment
    #   rawx, rawy - the RAW pixel coordinates, run from 1 to 600 for MOS
    #   ccdnr - the CCD number, run from 1 (central CCD) to 7
    #
    # HISTORY:
    #   Created 19 Aug 2019, Ivan Valtchanov, XMM SOC.
    #
    # Note: validated by a comparison with ecoordconv XMM-SAS task
    #
    if (detector != 'mos1'):
        print (f"Not yet implemented for {detector}")
        return None, None
    #
    # consistency checks
    #
    if (ccdnr < 1 or ccdnr > 7):
        print (f"CCD number for {detector} can only be in [1,7]")
        return None,None
    if (rawx < 1 or rawx > 600 or rawy < 1 or rawy > 600):
        print (f"RAWXY valuse for {detector} can only be in [1,600]")
        return None,None
    # CCD rotation signs
    signs = {1: [1,1], 2: [-1,1], 3: [-1,1], 4: [-1,1], 5: [1,-1], 6: [1,-1], 7: [1,-1]}
    #
    # reference DETXY values at RAWX=1 and RAWY=1
    #
    detxy0 = {1: [-6593.5733,-6593.5733], \
              2: [13099.951,-20147.005], \
              3: [19856.315,-6898.8552], \
              4: [13163.421,6428.8926], \
              5: [-13003.676,19693.992], \
              6: [-19770.078,6494.8849], \
              7: [-13117.889,-6864.6255]}

    # empirical slope coefficient
    mcoef = 22.01539
    # additional empirical offset coefficient for non-central CCDs
    coef = 9.219/600.0
    dxraw = rawx-1
    dyraw = rawy-1
    detx = signs[ccdnr][0]*mcoef*dxraw + detxy0[ccdnr][0] + dxraw*coef
    dety = signs[ccdnr][1]*mcoef*dyraw + detxy0[ccdnr][1] + dyraw*coef
    if (ccdnr > 1):
        # Swap RAWX,RAWY for non-central CCDs
        dxraw = (rawy-1)
        dyraw = (rawx-1)
        detx = signs[ccdnr][0]*mcoef*dxraw + detxy0[ccdnr][0] + signs[ccdnr][0]*dxraw*coef
        dety = signs[ccdnr][1]*mcoef*dyraw + detxy0[ccdnr][1] + signs[ccdnr][1]*dyraw*coef
    #
    return (detx,dety)

def get_theta_phi(detector='mos1',rawx=1,rawy=1,ccdnr=1,which='optical'):
    #
    # PURPOSE:
    #   derive the off-axis distance (theta) and the azymuthal angle (phi)
    # INPUTS:
    #   detector - the XMM detector (only implemented for MOS1 at the moment)
    #   rawx, rawy - the RAW pixel coordinates
    #   ccdnr - the CCD number
    #   which - which centre to use for reference:
    #   'optical' - (default) the optical axis centre in RAWX,RAWY coordinates
    #           taken from XMM_MISCDATA_0022.CCF
    #   'centre' - RAWX,RAWY = (300,300), the centre of CCD1
    #
    # OUTPUTS:
    #   theta in arcsec, phi is in radians
    # HISTORY:
    #   Created 19 Aug 2019, Ivan Valtchanov, XMM SOC.
    #
    # Note: validated by a comparison with ecoordconv XMM-SAS task
    #
    if (detector != 'mos1'):
        print (f"Not yet implemented for {detector}")
        return None, None
    #
    # consistency checks
    #
    if (ccdnr < 1 or ccdnr > 7):
        print (f"CCD number for {detector} can only be in [1,7]")
        return None,None
    if (rawx < 1 or rawx > 600 or rawy < 1 or rawy > 600):
        print (f"RAWXY valuse for {detector} can only be in [1,600]")
        return None,None
    #
    if (which == 'optical'):
        adetx0,adety0 = raw2det(detector=detector,rawx=305,rawy=291,ccdnr=1)
    else:
        adetx0,adety0 = raw2det(detector=detector,rawx=300,rawy=300,ccdnr=1)
    #
    adetx,adety = raw2det(detector=detector,rawx=rawx,rawy=rawy,ccdnr=ccdnr)
    #
    # now the theta in arcsec
    theta = np.sqrt((adety-adety0)**2 + (adetx-adetx0)**2)*0.05
    phi = (np.arctan2(adety-adety0,adetx-adetx0) - np.pi/2.0) % (2*np.pi)
    #
    return theta,phi
