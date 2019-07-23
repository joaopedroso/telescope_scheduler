from constants import SPEED

def move_time(sky1, sky2):
    """
    Returns the time required for moving the telescope from SkyCoord `sky1` to `sky2`

    :rtype: float
    """
    ra = abs(sky1.ra.degree - sky2.ra.degree)
    ra = min(ra, 360 - ra)
    dec = abs(sky1.dec.degree - sky2.dec.degree)
    return max(ra, dec) * SPEED
