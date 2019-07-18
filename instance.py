import math
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import Latitude, Longitude

# for getting moon coordinates:
from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_body

from constants import *
from astro_tools import calc_altaz, is_hidden_by_moon


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
        for nnn in np.arange(int((360.0+overlap)/((tomoe_fov-overlap)/math.cos(tmpdec*math.pi/180.0))+1)):
            if nnn == 0:
                tmpra = ra_start + nnn * (tomoe_fov)
            else:
                tmpra = ra_start + nnn * (tomoe_fov - overlap)/math.cos(tmpdec*math.pi/180.0)
            if tmpra < 360.0 and tmpdec < 90.0:
                tmpra_list.append(tmpra)
                tmpdec_list.append(tmpdec)
                tmpfov_list.append(fovnum)
                c_icrs = SkyCoord(ra=tmpra*u.degree, dec=tmpdec*u.degree, frame='icrs')
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
        for nnn in np.arange(int((360.0+overlap)/((tomoe_fov-overlap)/math.cos(tmpdec*math.pi/180.0))+1)):
            if nnn == 0:
                tmpra = ra_start + nnn * (tomoe_fov)
            else:
                tmpra = ra_start + nnn * (tomoe_fov - overlap)/math.cos(tmpdec*math.pi/180.0)
            if tmpra < 360.0:
                tmpra_list.append(tmpra)
                tmpdec_list.append(tmpdec)
                tmpfov_list.append(fovnum)
                c_icrs = SkyCoord(ra=tmpra*u.degree, dec=tmpdec*u.degree, frame='icrs')
                lb_galactic = c_icrs.galactic
                tmplll = lb_galactic.l.degree
                tmpbbb = lb_galactic.b.degree
                tmplll_list.append(tmplll)
                tmpbbb_list.append(tmpbbb)
                fovnum += 1
        tmpdec -= (tomoe_fov - overlap) * math.sqrt(3.0) / 2.0

    return tmpra_list, tmpdec_list, tmpfov_list, fovnum, tmplll_list, tmpbbb_list



# define time interval for the observations
def mk_obs_time():
    ### time
    utcoffset = +9 * u.hour  # JST
    time_start = Time('2018-11-21T17:00:00') - utcoffset
    time_end = Time('2018-11-22T07:00:00') - utcoffset
    return time_start,time_end



# build set of locations to observe
def mk_obs_set():
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

    ### preparetion for the conversion (RA,Dec)==>(El,Az)
    ra = Longitude(ra_list, unit=u.deg)
    dec = Latitude(dec_list, unit=u.deg)
    sky = SkyCoord(ra, dec, frame='icrs')

    # exclude points that are not observable (checking with intervals of 15 minutes)
    time_start, time_end = mk_obs_time()
    t = time_start  # time cursor, from begin to end of night
    visible = set()
    while t <= time_end:
        alt, az = calc_altaz(t, kiso, sky)
        for k in range(len(sky)):
            if k in visible:
                continue
            if alt[k] >= MIN_OBS_ANGLE:
                visible.add(k)
        t += 15 * u.minute
    visible = sorted(list(visible))

    # remove locations hidden by the moon
    # moon location
    t = time_start + (time_end - time_start) / 2
    with solar_system_ephemeris.set('builtin'):
        moon = get_body('moon', t, kiso)
    visible = [i for i in visible if not is_hidden_by_moon(moon, sky[i], kiso, t, within_angle=ANGLE_FROM_MOON)]

    ### preparetion for the conversion (RA,Dec)==>(El,Az)
    K = range(len(ra_list))
    ra = Longitude([ra_list[k] for k in K if k in visible], unit=u.deg)
    dec = Latitude([dec_list[k] for k in K if k in visible], unit=u.deg)
    sky = SkyCoord(ra, dec, frame='icrs')

    return sky, time_start, time_end



if __name__ == "__main__":
    sky, time_start, time_end = mk_obs_set()
    print("observation set of locations:")
    for i in range(len(sky)):
        print("{}\t({},{})".format(i, sky[i].ra.degree, sky[i].dec.degree))
