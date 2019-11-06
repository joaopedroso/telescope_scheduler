"""
simulator: Example of usage of the telescope scheduler

This program illustrates how to use the telescope scheduler implemented
in module `scheduler.py`.

A dummy function `move_telescope` simply prints the position where
the real-world telescope should go.

The simulation starts at the begin of the observation period specified in
`instance.py`, making virtual observations of sky positions specified there
until reaching the end of the observation period.

For adapting to a different setting:
    - Rewrite file `instance.py`, or otherwise specify functions `mk_obs_time`
      and `mk_obs_set` for specifying a different instance.
    - Use function `hidden` specified in file `hidden.py` for updating in
      real-time the positions that are not available for observation (eg, due
      to being cluttered by clouds).
    - Define the telescope's movement speed in file `motion_time.py`.
"""
if __name__ == '__main__':
    from hidden import hidden
    from constants import *
    from obs_data import ObsData
    from solution import Solution
    from scheduler import scheduler

    TIME_ADVANCE_IF_NO_OBSERVABLE = 15  # time to wait if there are no observable positions (eg, bad weather)
    TIME_LIM = 15  # solving time limit (in seconds)
    print("starting simulation, using {} s CPU time limit for scheduling".format(TIME_LIM))


    def move_telescope(sky):
        print("telescope is now in position ({:5.1f}, {:5.1f})".format(sky.ra.degree, sky.dec.degree))


    # prepare/load data
    from instance import mk_obs_time, mk_obs_set

    time_start, time_end = mk_obs_time()
    sky = mk_obs_set()
    obs_data = ObsData(sky, time_start, time_end)

    # # if solving multiple times the same instance, pickling preprocessed data may be convenient:
    # # comment the previous lines, and uncomment the following for re-loading data quickly:
    # import pickle
    # filename = "telescope_sched_state.pck"
    # from os.path import exists
    # state = False
    # if exists(filename):
    #     with open(filename, 'rb') as f:
    #         sky, time_start, time_end, obs_data = pickle.load(f)
    #     state = True
    # time_start_, time_end_ = mk_obs_time()
    # if not state or time_start != time_start_ or time_end != time_end_:  # pickled data outdated
    #     time_start, time_end = mk_obs_time()
    #     sky = mk_obs_set()
    #     obs_data = ObsData(sky, time_start, time_end)
    #     with open(filename, 'wb') as f:
    #         # Pickle the 'data' dictionary using the highest protocol available.
    #         pickle.dump((sky, time_start, time_end, obs_data), f, pickle.HIGHEST_PROTOCOL)

    # start simulation of an observation night
    t = time_start
    earliest = {k: t for k in obs_data.visible}
    # construct initial (preliminary) solution
    sol = scheduler(obs_data, time_now=t, last_obs=[None], earliest=earliest,
                    current_best=None, hidden_list=[], time_lim=TIME_LIM)

    # time_end = time_start + 120 * u.second    # use for quick tests


    curr = None
    tpos = sol.seq[1]  # telescope position for the first observation
    print()
    print("PRELIMINARY SOLUTION:")
    print("{}\nn.obs: {}\tpos: {}\t3obs/total: {}/{}\t{}".format(t, 0, tpos, sol.n3obs, len(sol.seq), sol.values))
    print("current solution: {} ...".format(sol.seq[:10]))
    move_telescope(sky[tpos])
    print()
    seq = [tpos]
    nobs = 0
    best = Solution(seq, obs_data.K)
    while t <= time_end:
        nobs += 1
        print("Observation number {}, current time: {}".format(nobs, t))
        covered = hidden(obs_data.sky)
        print("Cluttered positions:", covered)
        print("Not observable:", [k for k in obs_data.visible if earliest[k] > t])
        sol = scheduler(obs_data, time_now=t, last_obs=seq, earliest=earliest,
                        current_best=best, hidden_list=covered, time_lim=TIME_LIM)
        if len(sol.seq) < 2:
            print("no observable positions, waiting {} seconds".format(TIME_ADVANCE_IF_NO_OBSERVABLE))
            t += TIME_ADVANCE_IF_NO_OBSERVABLE * u.second
            continue
        curr = sol.seq[1]
        print("\tcurrent solution: {} ...".format(sol.seq[:10]))
        print("\tnext: {}\t3obs/total: {}/{}\tobs.card: {}".format(curr, sol.n3obs, len(sol.seq), sol.values))
        print("Simulating telescope movement to position {}: ".format(curr), end="\t")
        move_telescope(sky[curr])

        # simulate time advance
        move = obs_data.move_time[tpos, curr]
        delay = move + EXPOSURE
        t += delay * u.second
        earliest[curr] = t + MIN_OBS_GAP

        # update variables
        best = Solution(sol.seq[1:], obs_data.K)  # starting solution for next iteration
        tpos = curr
        seq.append(curr)
        print("Previously observed positions:", seq[1:])
        print()

    print("SUMMARY OF OBSERVATIONS:")
    final = Solution(seq, obs_data.K)
    print("{}\nn.obs: {}\tpos: {}\t3obs/total: {}/{}\t{}".format(t, nobs, curr, final.n3obs, len(seq), final.values))
