"""
obs_data: Hold and preprocess information with experiment's details

classes:
    - ObsData(): convenience class for holding all information for one observation night
"""

from astropy.coordinates import SkyCoord
from constants import *
from motion_time import motion_time
# for getting moon coordinates:
from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_body
from astro_tools import calc_altaz, is_hidden_by_moon


class ObsData:
    def __init__(self, sky, time_start, time_end):
        """
        Convenience class holding information for an observation night

        Data to be stored in `self`:

        :ivar sky: list of `SkyCoord` positions to be observed
        :ivar K: list of indices in `sky`
        :ivar visible: indices for positions which become visible during the night
        :ivar total_obs_seconds: reference time for the end of the experiment (in seconds)
        :ivar obs_gap_seconds: minimum interval (in seconds) between two observations of the same position
        :ivar ivisible: subset of initially visible `sky` indices
        :ivar tvisible: dictionary, `tvisible[k]` -> time at which `k` appears in the sky (in seconds after starting time)
        :ivar thidden: dictionary, `thidden[k]` ->  time at which `k` disappears
        :ivar move_time: dictionary, `move_time[i,j]` ->  time to move from `self.sky[i]` to `self.sky[j]`

        :param list of SkyCoord sky: list of `SkyCoord` positions to be observed
        :param astropy.time.Time time_start: earliest time for an observation
        :param astropy.time.Time time_end: latest time for an observation
        """
        self.sky = sky.copy()
        self.K = range(len(sky))
        self.time_start = time_start
        self.time_end = time_end
        self.obs_gap_seconds = MIN_OBS_GAP.value * 3600
        self.total_obs_seconds = (time_end - time_start).value * 24 * 3600  # time (s) w.r.t. begin of experiment

        # ivisible: subset of initially visible `sky` indices
        # tvisible[k]: time at which k appears in the sky (in seconds after starting time)
        # thidden[k]: time at which k disappears
        self.ivisible, self.tvisible, self.thidden = self.preprocess()

        # remove locations hidden by the moon
        # moon location
        t = time_start + (time_end - time_start) / 2
        with solar_system_ephemeris.set('builtin'):
            moon = get_body('moon', t, telescope_pos)
        self.visible = [i for i in self.tvisible \
                        if not is_hidden_by_moon(moon, sky[i], telescope_pos, t, within_angle=ANGLE_FROM_MOON)]

        self.move_time = self.calc_times()

    def preprocess(self):
        """
        Prepare data for the scheduler

        :return: a tuple containing data to be stored in self
        """

        #
        t = self.time_start
        tvisible = {}  # tvisible[k] = time at which k appears in the sky
        thidden = {}  # thidden[k] = time at which k disappears
        alt, az = calc_altaz(t, telescope_pos, self.sky)
        visible = [k for k in self.K if alt[k] >= MIN_OBS_ANGLE]
        ivisible = visible.copy()
        new_visible = None
        while t <= self.time_end:
            if new_visible != None:
                changed = set(new_visible) ^ set(visible)  # set of locations for which visibility changed
                assert len(changed) > 0
                visible = new_visible
            for k in visible:
                if k not in tvisible:
                    tvisible[k] = (t - self.time_start).value * 24 * 3600
            for k in tvisible:
                if k not in visible and k not in thidden:
                    thidden[k] = (t - self.time_start).value * 24 * 3600

            # hunt search for exact moment when some point changes visibility
            step = 1
            while True:
                t_ub = t + step * u.second
                alt, az = calc_altaz(t_ub, telescope_pos, self.sky)
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
                        alt, az = calc_altaz(t_med, telescope_pos, self.sky)
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
                thidden[k] = (self.time_end - self.time_start).value * 24 * 3600
        return ivisible, tvisible, thidden

    def calc_times(self):
        """
        Returns a "distance" matrix: time for moving between each pair of positions in self.sky
        """
        dSKY = {}
        for k in self.visible:
            for k_ in self.visible:
                dSKY[k, k_] = motion_time(self.sky[k], self.sky[k_])
        # provide 0 time from None to any position, for the first move
        for k in self.visible:
            dSKY[None, k] = 0

        return dSKY


if __name__ == "__main__":
    from instance import mk_obs_time, mk_obs_set

    time_start, time_end = mk_obs_time()
    sky = mk_obs_set()
    obs_data = ObsData(sky, time_start, time_end)
    print("observation time from {} ({}) to {} ({})".format(time_start, obs_data.time_start,
                                                            time_end, obs_data.time_end))
    print("{} positions to be observed from a total of {}".format(len(obs_data.visible), len(obs_data.K)))
    for k in obs_data.visible:
        print(k, obs_data.tvisible[k], obs_data.thidden[k])
