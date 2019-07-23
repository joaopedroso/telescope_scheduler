"""
instance: Define observation's start/end time and positions

Functions:
    - mk_obs_time(): returns tuple with starting and ending time
    - mk_obs_set(): return list with SkyCoord positions to observe

Rewrite these functions for specifying a different instance.
"""
import math
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import Latitude, Longitude
from constants import *


### define Tomo-e's observing regions (fixed in the RA/Dec coordinate)
###    (RA, Dec) <==> (b, l) conversion is also calculated.
def tomoe_fov_radec(ra_start, dec_start, tomoe_fov, overlap):
    tmpra_list = []
    tmpdec_list = []
    tmplll_list = []
    tmpbbb_list = []
    tmpfov_list = []

    fovnum = 0
    ### Dec >= 0
    tmpdec = dec_start
    while tmpdec >= -0.01 and tmpdec < 90.0:
        # print (tmpdec)
        for nnn in np.arange(int((360.0 + overlap) / ((tomoe_fov - overlap) / math.cos(tmpdec * math.pi / 180.0)) + 1)):
            if nnn == 0:
                tmpra = ra_start + nnn * (tomoe_fov)
            else:
                tmpra = ra_start + nnn * (tomoe_fov - overlap) / math.cos(tmpdec * math.pi / 180.0)
            if tmpra < 360.0 and tmpdec < 90.0:
                tmpra_list.append(tmpra)
                tmpdec_list.append(tmpdec)
                tmpfov_list.append(fovnum)
                c_icrs = SkyCoord(ra=tmpra * u.degree, dec=tmpdec * u.degree, frame='icrs')
                lb_galactic = c_icrs.galactic
                tmplll = lb_galactic.l.degree
                tmpbbb = lb_galactic.b.degree
                tmplll_list.append(tmplll)
                tmpbbb_list.append(tmpbbb)
                fovnum += 1
        tmpdec += (tomoe_fov - overlap) * math.sqrt(3.0) / 2.0
    ### Dec < 0
    tmpdec = dec_start - (tomoe_fov - overlap) * math.sqrt(3.0) / 2.0
    while tmpdec >= -60.0:
        # print (tmpdec)
        for nnn in np.arange(int((360.0 + overlap) / ((tomoe_fov - overlap) / math.cos(tmpdec * math.pi / 180.0)) + 1)):
            if nnn == 0:
                tmpra = ra_start + nnn * (tomoe_fov)
            else:
                tmpra = ra_start + nnn * (tomoe_fov - overlap) / math.cos(tmpdec * math.pi / 180.0)
            if tmpra < 360.0:
                tmpra_list.append(tmpra)
                tmpdec_list.append(tmpdec)
                tmpfov_list.append(fovnum)
                c_icrs = SkyCoord(ra=tmpra * u.degree, dec=tmpdec * u.degree, frame='icrs')
                lb_galactic = c_icrs.galactic
                tmplll = lb_galactic.l.degree
                tmpbbb = lb_galactic.b.degree
                tmplll_list.append(tmplll)
                tmpbbb_list.append(tmpbbb)
                fovnum += 1
        tmpdec -= (tomoe_fov - overlap) * math.sqrt(3.0) / 2.0

    return tmpra_list, tmpdec_list, tmpfov_list, fovnum, tmplll_list, tmpbbb_list


def mk_obs_time():
    """
    Define time interval for the observations

    :returns: tuple with starting and ending time
    """
    utcoffset = +9 * u.hour  # JST
    time_start = Time('2018-11-21T17:00:00') - utcoffset
    time_end = Time('2018-11-22T07:00:00') - utcoffset
    return time_start, time_end


def mk_obs_set():
    """
    Build set of locations to observe

    :returns: list with SkyCoord positions to observe
    """
    ### inital value
    ra_start = 0.0
    dec_start = 0.0

    ### tomoe field-of-view radius [deg] (and overlaps between adjacent regions)
    tomoerad = 4.4
    tomoe_fov = tomoerad * 2.0
    overlap = 0.223256  # 0.1

    ### Tomo-e regions defined here
    ra_list, dec_list, fovid_list, fovnum, lll_list, bbb_list = \
        tomoe_fov_radec(ra_start, dec_start, tomoe_fov, overlap)

    ### preparation for the conversion (RA,Dec)==>(El,Az)
    ra = Longitude(ra_list, unit=u.deg)
    dec = Latitude(dec_list, unit=u.deg)
    sky = SkyCoord(ra, dec, frame='icrs')

    return sky


if __name__ == "__main__":
    time_start, time_end = mk_obs_time()
    sky = mk_obs_set()
    print("set of observation locations:")
    for i in range(len(sky)):
        print("{}\t({},{})".format(i, sky[i].ra.degree, sky[i].dec.degree))
    print("starting {}, finishing {}".format(time_start, time_end))
