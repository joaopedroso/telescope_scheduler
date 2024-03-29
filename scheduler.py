"""
scheduler: Use a nearest-neighbor variant to obtain an observation schedule

Function `scheduler` is the main tool in this software.
The aim is to obtain an observation schedule for a high-cadence telescope,
where the objective is to maximize the number of positions in the sky,
from a given list, which are observed at least three times.

See `simulator.py` or `main` part below for examples of usage.
"""

import random
from time import time
from datetime import timedelta
from constants import *
from solution import Solution

PRINTED_IDLING_MESSAGE = False  # to avoid cluttering the output with warning messages


def scheduler(target_obs, obs_data, time_now, last_obs, earliest, current_best, hidden_list, time_lim):
    """
    Obtain a schedule based on a variant of the nearest-neighbor heuristic

    Repeatedly construct observation schedules, until `time_lim` is reached.
    Each construction is based on a simulation, where at each point in time
    an observation is chosen among `obs_data.sky` positions available which
    are closest to the current position of the telescope.

    :rtype: object of class `Solution`, specifying observation sequence and its evaluation
    :param target_obs: number of desired observations (repetitions) for each position
    :param obs_data: data specifying the instance to be solved: sky positions and start/end times
    :param time_now: current time
    :param last_obs: sequence of observations made before
    :param earliest: dictionary, position `k` cannot be observed before `earliest[k]`
    :param current_best: previously best known solution (`None` if not available)
    :param hidden_list: indices in `obs_data.sky` that cannot be observed now (see `hidden.py`)
    :param time_lim: clock time allowed for obtaining a solution (in seconds) (use 0 for only one construction)
    :return: best sequence of observations found
    """
    cpu0 = time()
    d = obs_data  # typing convenience
    now = (time_now - d.time_start).value * 24 * 3600  # simulation time
    earl = {k: ((earliest[k] - d.time_start).value * 24 * 3600) for k in earliest}
    if LOG:
        print("{:6} {:>9}: {:6}".format("point", "sim.time", "#obs"))

    visible = set(k for k in d.visible if earl[k] <= now and d.tvisible[k] <= now <= d.thidden[k])
    cluttered = visible.intersection(set(hidden_list))

    if current_best is None or len(cluttered) > 0:  # if initial observation, or if previous best is invalid
        best_nXobs = -1
    else:  # previous solution may be used
        best = current_best
        best_nXobs = current_best.nXobs

    niter = 0  # number of iterations
    while True:
        niter += 1
        t = now
        start = earl.copy()  # start[k]: earliest time that 'k' may be observed
        obs = {k: 0 for k in d.K}  # number of observations
        prev = last_obs[-1]
        obs_seq = [prev]  # sequence of observations
        if LOG:
            print("{:6} {:>9}: {:6} {:6} {:^13} -> {:^14} | {!s:>8} {:>9} {:>10}".format(
                "point", "time", "#obs", "#cand", "sky", "telescope", "previous", "movement", "sim time"))
        while t <= d.total_obs_seconds:
            # check CPU limit was reached:
            if time_lim > EPS and time() - cpu0 >= time_lim:
                return best

            # add next observation
            visible = [k for k in d.visible if start[k] <= t and \
                       d.tvisible[k] <= t <= d.thidden[k] and k not in hidden_list]
            if len(visible) == 0:
                global PRINTED_IDLING_MESSAGE
                if not PRINTED_IDLING_MESSAGE:
                    PRINTED_IDLING_MESSAGE = True
                    print("WARNING: scheduler: NN has no available points for observing, advancing time {} seconds"\
                          .format(MAGICN4))
                t += MAGICN4
                continue

            # build candidate list for nearest neighbors
            mindist = INFINITY
            mindist2 = INFINITY  # for having more than one candidate, keep track of 2nd-NN
            minobs = min(obs[k] for k in visible)
            # has12obs = any(obs[k]==1 or obs[k]==2 for k in visible)
            hasobs = any(obs[k]>0 and obs[k]<target_obs for k in visible)
            cand = []
            for k in visible:
                kobs = obs[k]
                # # old heuristic:
                # if kobs > minobs + MAGICN3 or (kobs >= 3 and minobs <= 2):  # avoid 4-th time observations
                #     continue
                kdist = d.move_time[prev, k]
                #
                # new heuristic:
                if kobs >= target_obs and minobs <= target_obs-1:
                    continue
                if kobs == 0 and hasobs:
                    continue
                #
                if kdist < mindist - EPS:
                    if len(cand) > 0:
                        k2 = cand[0]
                        mindist2 = mindist
                    mindist = kdist
                    cand = [k]
                elif kdist < mindist + EPS:
                    mindist = kdist
                    cand.append(k)
                elif kdist < mindist2 - EPS:
                    mindist2 = kdist
                    k2 = k
            if len(cand) < 2 and mindist2 < MAGICN5 * mindist - EPS:  # avoid having only one candidate
                cand.append(k2)

            # choose a (random) candidate from the list
            curr = random.choice(cand)
            move = d.move_time[prev, curr]
            delay = move + EXPOSURE
            if LOG:
                from astro_tools import calc_altaz
                import astropy.units as u
                alt, az = calc_altaz(d.time_start + t * u.second, telescope_pos, d.sky)
                print("{:6} {:9.2f}: {:6} {:6} ({:5.1f},{:5.1f}) -> ({:5.1f}, {:5.1f}) | {!s:>7} {:9.2f} {:>11}".format(
                    curr, delay, len(visible), len(cand), d.sky[curr].ra.degree, d.sky[curr].dec.degree,
                    alt[curr], az[curr],
                    prev, move, str(timedelta(seconds=round(t, 0)))))

            if t + delay > d.total_obs_seconds:
                break
            t += delay
            start[curr] = t + delay + d.obs_gap_seconds
            obs[curr] += 1
            obs_seq.append(curr)
            prev = curr

        sol = Solution(target_obs, obs_seq, obs_data.K)
        assert len(obs_seq) - 1 == sum(sol.obs[k] for k in sol.obs)  # -1 => first "observation" is None
        #     if len(obs_seq) > best_n_obs:
        if sol.nXobs > best_nXobs:
            best_nXobs = sol.nXobs
            best = sol.copy()

        # values = [sum(1 for k in K if obs[k] == v) for v in range(n_max)]
        if LOG:
            print("{:6}\t{:9.2f}\t{:6}\t{:6}\t{}".format(niter, t, best.nXobs, len(best.seq), best.values))

        if time_lim <= 0:  # make only one construction
            return best


if __name__ == "__main__":
    import sys
    from astropy.time import Time
    import astropy.units as u
    from obs_data import ObsData
    from astro_tools import calc_altaz

    import importlib

    try:
        inst  = importlib.import_module(sys.argv[1])
        target_obs = int(sys.argv[2])
        time_lim = float(sys.argv[3])
    except:
        print("usage: python {} instance time\n"
              "where:\n"
              "  instance[.py] - file defining mk_obs_time and mk_obs_set\n"
              "  nobs - target number of observations per position\n"
              "  time - computing time allowed for finding a solution (s)\n"
              "e.g.:\n"
              "  python {} instance 3 60\n".format(sys.argv[0],sys.argv[0]))
        exit(-1)

    
    # prepare data
    time_start, time_end = inst.mk_obs_time()
    sky = inst.mk_obs_set()
    obs_data = ObsData(sky, time_start, time_end)

    t = time_start
    earliest = {k: t for k in obs_data.visible}
    PRINTED_IDLING_MESSAGE = True
    sol = scheduler(target_obs, obs_data, time_now=t, last_obs=[None], earliest=earliest,
                    current_best=None, hidden_list=[], time_lim=time_lim)


    prev = None
    t = 0
    for i,curr in enumerate(sol.seq):
        if i > 0:
            # if t < obs_data.tvisible[curr]:
            #     print("[delaying {} seconds, observing {} ({:5.1f},{:5.1f})] : ".format(
            #         MAGICN4, curr, obs_data.sky[curr].ra.degree, obs_data.sky[curr].dec.degree))
            move = obs_data.move_time[prev, curr]
            clock = obs_data.time_start + t * u.second

            alt, az = calc_altaz(clock, telescope_pos, obs_data.sky)
            # print("{:.7f}, {:.7f}".format(alt[curr], az[curr]))
            print("{}, {:.7f}, {:.7f}, {}".format(i, sky[curr].ra.degree, sky[curr].dec.degree, clock+9*u.hour))
            t += move + EXPOSURE
        prev = curr
