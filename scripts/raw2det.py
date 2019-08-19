

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
