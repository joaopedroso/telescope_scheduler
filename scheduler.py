import random
from time import clock
from datetime import timedelta
from constants import *
from solution import Solution

PRINTED_IDLING_MESSAGE = False  # to avoid cluttering the output with warning messages


def scheduler(obs_data, time_now, last_obs, earliest, current_best, hidden_list, time_lim):
    cpu0 = clock()
    d = obs_data  # typing convenience
    now = (time_now - d.time_start).value * 24 * 3600  # simulation time
    earl = {k: ((earliest[k] - d.time_start).value * 24 * 3600) for k in earliest}

    if LOG:
        print("{:6} {:>9}: {:6}".format("point", "sim.time", "#obs"))

    visible = set(k for k in d.K if earl[k] <= now and d.tvisible[k] <= now <= d.thidden[k])
    cluttered = visible.intersection(set(hidden_list))

    if current_best is None or len(cluttered) > 0:  # initial observation, or previous best invalid
        best_n3obs = -1
    else:   # previous solution may be used
        best = current_best
        best_n3obs = current_best.n3obs
    # best_n3obs = -1
    niter = 0  # number of iterations
    while True:
        niter += 1
        t = now
        start = earl.copy()  # start[k]: earliest time that 'k' may be observed
        obs = {k: 0 for k in d.K}  # number of observations
        prev = last_obs[-1]
        obs_seq = [prev]  # sequence of observations
        obs_times = [None]
        if LOG:
            print("{:6} {:>9}: {:6} {:6} {:^13} -> {:^14} | {!s:>8} {:>9} {:>10}".format(
                "point", "time", "#obs", "#cand", "sky", "kiso", "previous", "movement", "sim time"))
        while t <= d.total_obs_seconds:
            # check CPU limit was reached:
            if clock() - cpu0 >= time_lim:
                return best

            # add next observation
            visible = [k for k in d.K if start[k] <= t and d.tvisible[k] <= t <= d.thidden[k] and k not in hidden_list]
            if len(visible) == 0:
                global PRINTED_IDLING_MESSAGE
                if not PRINTED_IDLING_MESSAGE:
                    PRINTED_IDLING_MESSAGE = True
                    print("WARNING: no available points for observing, waiting {} seconds".format(MAGICN4))
                t += MAGICN4
                continue

            mindist = INFINITY
            minobs = min(obs[k] for k in visible)
            for k in visible:
                if obs[k] > minobs + MAGICN3:
                    continue
                kdist = d.move_time[prev, k]
                if kdist < mindist - EPS:
                    mindist = kdist
                    cand = [k]
                elif kdist < mindist + EPS:
                    mindist = kdist
                    cand.append(k)

            curr = random.choice(cand)
            move = mindist
            delay = move + EXPOSURE
            if LOG:
                from astro_tools import calc_altaz
                alt, az = calc_altaz(d.time_start + t * u.seconds, kiso, d.sky)
                print("{:6} {:9.2f}: {:6} {:6} ({:5.1f},{:5.1f}) -> ({:5.1f}, {:5.1f}) | {!s:>7} {:9.2f} {:>11}".format(
                    curr, delay, len(visible), len(cand), d.sky[curr].ra.degree, d.sky[curr].dec.degree,
                    alt[curr], az[curr],
                    prev, move, str(timedelta(seconds=round(t, 0)))))

            if t + delay > d.total_obs_seconds:
                break
            obs_times.append(t)
            t += delay
            start[curr] = t + delay + d.delay_seconds
            obs[curr] += 1
            obs_seq.append(curr)
            prev = curr

        sol = Solution(obs_seq, obs_data.K)
        assert len(obs_seq) - 1 == sum(sol.obs[k] for k in sol.obs)  # -1 => first "observation" is None
        #     if len(obs_seq) > best_n_obs:
        if sol.n3obs > best_n3obs:
            best_n3obs = sol.n3obs
            best = sol.copy()

        # values = [sum(1 for k in K if obs[k] == v) for v in range(n_max)]
        if LOG:
            print("{:6}\t{:9.2f}\t{:6}\t{:6}\t{}".format(niter, t, best.n3obs, len(best.seq), best.values))
