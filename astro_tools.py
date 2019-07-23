"""
astro_tools: Useful tools for telescope scheduling based on astropy
"""

import math
import numpy as np
from astropy.coordinates import EarthLocation, AltAz


def calc_altaz(time, loc, obj):
    """
    Convert RA/Dec to Alt/Az

    :param time: time of the observation
    :param loc: array with positions
    :param obj: array with positions in SkyCoord coordinated
    :return: tuple with arrays of positions in alt and az coordinates
    """
    azuse = []
    tmpaltaz = obj.transform_to(AltAz(obstime=time, location=loc))
    tmpalt = tmpaltaz.alt.degree
    tmpaz = tmpaltaz.az.degree + 90.0
    for azazaz in np.arange(len(tmpaz)):
        if tmpaz[azazaz] > 360.0:
            tmptmpaz = tmpaz[azazaz] - 360.0
        else:
            tmptmpaz = tmpaz[azazaz]
        azuse.append(tmptmpaz)
    altuse = tmpalt
    return altuse, azuse


def is_within_angle(obj_alt, obj_az, point_alt, point_az, within_angle):
    """Determine if a position is hidden by some object

    :param obj_alt: interfering object's position (alt)
    :param obj_az: interfering object's position (az)
    :param point_alt: position of point being observed (alt)
    :param point_az: position of point being observed (az)
    :param within_angle: defines zone around center of interfering object to avoid
    :returns: True if there is interference, False otherwise
    """
    angle = math.acos(
        math.cos(obj_alt) * math.cos(point_alt)
        + math.sin(obj_alt) * math.sin(point_alt) * math.cos(obj_az - point_az)
    )

    return angle * 180 / math.pi < within_angle


def is_hidden_by_moon(moon, point, location, obstime, within_angle=0.25):
    """Determine if a position is hidden by the moon

    :param moon: approximate coordinates of the moon, as returned by astropy.coordinates.get_body
    :param point: sky position being observed (in SkyCoord coordinates)
    :param location: EarthLocation from which observation is made (eg, Kiso)
    :param obstime: observation time
    :param within_angle: defines zone around center of interfering object to avoid
    :returns: True if there is interference, False otherwise
    """

    # `moon` should be obtained with:
    # with solar_system_ephemeris.set('builtin'):
    # 	moon = get_body('moon', obstime, telescope_pos)

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
