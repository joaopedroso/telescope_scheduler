"""
obs_data: module to hold and preprocess observation data

classes:
    - ObsData(): convenience class for holding all information for one observation night
"""

import astropy.units as u
from astropy.coordinates import SkyCoord
from constants import *
from astro_tools import calc_altaz
from move_time import move_time


class ObsData:
    def __init__(self, sky, time_start, time_end):
        self.sky = sky.copy()
        self.K = range(len(sky))
        self.time_start = 0   # set reference time to begin of the experiment (for code speedup)
        self.time_end = (time_end - time_start).value * 24 * 3600
        self.obs_gap_seconds = MIN_OBS_GAP.value * 3600
        self.total_obs_seconds = (time_end - time_start).value * 24 * 3600

        # ivisible: subset of initially visible 'sky' indices
        # tvisible[k]: time at which k appears in the sky (in seconds after starting time)
        # thidden[k]: time at which k disappears
        self.ivisible, self.tvisible, self.thidden = self.preprocess(time_start, time_end)

        self.visible = list(self.tvisible)
        self.move_time = self.calc_times()

    def preprocess(self, time_start, time_end):
        # prepare data for the optimization model
        t = time_start  # time cursor, from begin to end of night (different of self.time_start = 0.)
        tvisible = {}  # tvisible[k] = time at which k appears in the sky
        thidden = {}  # thidden[k] = time at which k disappears
        alt, az = calc_altaz(t, kiso, self.sky)
        visible = [k for k in self.K if alt[k] >= MIN_OBS_ANGLE]
        ivisible = visible.copy()
        new_visible = None
        while t <= time_end:
            if new_visible != None:
                changed = set(new_visible) ^ set(visible)  # set of locations for which visibility changed
                assert len(changed) > 0
                visible = new_visible
            for k in visible:
                if k not in tvisible:
                    tvisible[k] = (t - time_start).value * 24 * 3600
            for k in tvisible:
                if k not in visible and k not in thidden:
                    thidden[k] = (t - time_start).value * 24 * 3600

                    # hunt search for exact moment when some point changes visibility
            step = 1
            while True:
                t_ub = t + step * u.second
                alt, az = calc_altaz(t_ub, kiso, self.sky)
                tmp_visible = [k for k in self.K if alt[k] >= MIN_OBS_ANGLE]
                new_visible = tmp_visible
                if tmp_visible != visible:
                    # went too far, start bisection
                    while True:
                        step //= 2
                        if step == 0:
                            t = t_ub
                            break
                        t_med = t + step * u.second
                        alt, az = calc_altaz(t_med, kiso, self.sky)
                        tmp_visible = [k for k in self.K if alt[k] >= MIN_OBS_ANGLE]
                        if tmp_visible != visible:
                            t_ub = t_med
                            new_visible = tmp_visible
                        else:
                            t = t_med
                    break
                # otherwise, still within bounds -> double step ahead
                t += step * u.second
                step *= 2

        for k in tvisible:
            if k not in thidden:
                thidden[k] = (time_end - time_start).value * 24 * 3600
        return ivisible, tvisible, thidden

    def calc_times(self):
        dSKY = {}
        for k in self.visible:
            for k_ in self.visible:
                dSKY[k, k_] = move_time(self.sky[k], self.sky[k_])
        # provide 0 time from None to any position, for the first move
        for k in self.visible:
            dSKY[None, k] = 0

        return dSKY


if __name__ == "__main__":
    from astro_tools import mk_obs_set

    sky, time_start, time_end = mk_obs_set()
    obs_data = ObsData(sky, time_start, time_end)
    print("observation time from {} to {}".format(obs_data.time_start, obs_data.time_end))
