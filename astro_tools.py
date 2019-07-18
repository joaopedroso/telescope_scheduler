import math
import numpy as np
from astropy.coordinates import EarthLocation, AltAz

from constants import *

### convert RA/Dec to Alt/Az
def calc_altaz(time, loc, obj):
    altuse = []
    azuse = []
    tmpaltaz = obj.transform_to(AltAz(obstime=time,location=loc))
    tmpalt = tmpaltaz.alt.degree
    tmpaz = tmpaltaz.az.degree + 90.0
#	print (tmpaz)
    for azazaz in np.arange(len(tmpaz)):
        if tmpaz[azazaz] > 360.0:
            tmptmpaz = tmpaz[azazaz] - 360.0
        else:
            tmptmpaz = tmpaz[azazaz]
        azuse.append(tmptmpaz)
#	print (azuse)
    altuse = tmpalt
    return altuse, azuse
#	return tmpalt, tmpaz

# determine if a position is hidden by the moon
def is_within_angle(moon_alt, moon_az, point_alt, point_az, within_angle):
    """paremeters:
    - point: location of the point (SkyCoord) whose visibility we want to check
    - location: observation point (EarthLocation) (e.g., kiso)
    - obstime: observation time
    - within_angle: defines zone around moon center to avoid
    """
    angle = math.acos(
        math.cos(moon_alt) * math.cos(point_alt)
        + math.sin(moon_alt) * math.sin(point_alt) * math.cos(moon_az - point_az)
    )

    return angle * 180 / math.pi < within_angle


def is_hidden_by_moon(moon, point, location, obstime, within_angle=0.25):
    """paremeters:
    - point: location of the point (SkyCoord) whose visibility we want to check
    - location: observation point (EarthLocation) (e.g., kiso)
    - obstime: observation time
    - within_angle: defines zone around moon center to avoid
    """

    # with solar_system_ephemeris.set('builtin'):
    # 	moon = get_body('moon', obstime, kiso)

    moon_from_loc = moon.transform_to(AltAz(obstime=obstime, location=location))
    moon_alt = math.pi / 2 - moon_from_loc.alt.rad
    moon_az = moon_from_loc.az.rad

    point_from_loc = point.transform_to(AltAz(obstime=obstime, location=location))
    point_alt = math.pi / 2 - point_from_loc.alt.rad
    point_az = point_from_loc.az.rad

    angle = math.acos(
        math.cos(moon_alt) * math.cos(point_alt)
        + math.sin(moon_alt) * math.sin(point_alt) * math.cos(moon_az - point_az)
    )

    return angle * 180 / math.pi < within_angle

