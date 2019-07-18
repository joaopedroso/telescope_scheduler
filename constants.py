from astropy.coordinates import EarthLocation
import astropy.units as u

# Kiso location
lat_kiso = 35.783
lng_kiso = +137.616
height_kiso = 1130
kiso = EarthLocation(lat=lat_kiso*u.deg, lon=lng_kiso*u.deg, height=height_kiso*u.m)

# some observational constants:
EXPOSURE = 4 * 6 + 3 * 8   # 4 exposures, 3 movements (for using the 4 sensors
SPEED = 14./9   # telescope moving speed
MIN_OBS_ANGLE = 40.   # used to determine observable positions
MIN_OBS_GAP = 1.5 * u.hour    # minimum interval between observations of the same location
ANGLE_FROM_MOON = 30

# code-related constants
INFINITY = float("inf")
MAGICN3 = 2  # number of additional visits wrt minimum for rejecting
MAGICN4 = 5  # number of seconds to wait if not visiting position is available

LOG = False  # whether to print verbose output
EPS = 1.e-4  # tolerance for floating point operations
