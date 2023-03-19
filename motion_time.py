"""
motion_time: Define telescope's movement speed
"""

from constants import DOME_SPEED, DOME_DELAY, TELESCOPE_SPEED, TELESCOPE_DELAY


def motion_time(sky1, sky2):
    """
    Returns the time required for moving the telescope from SkyCoord `sky1` to `sky2`

    :param sky1: initial position
    :param sky2: target position
    :return:
    :rtype: float
    """
    ra = abs(sky1.ra.degree - sky2.ra.degree)
    ra = min(ra, 360 - ra)
    telescope_time = TELESCOPE_SPEED * ra + TELESCOPE_DELAY

    dec = abs(sky1.dec.degree - sky2.dec.degree)
    dome_time = DOME_SPEED * dec + DOME_DELAY

    return max(telescope_time, dome_time)
