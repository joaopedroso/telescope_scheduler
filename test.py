"""
test: Evaluate the nearest neighbor heuristic on some instances

Evaluates the telescope scheduler implemented in module `scheduler.py`.

"""
if __name__ == '__main__':
    from hidden import hidden
    from constants import *
    from obs_data import ObsData
    from solution import Solution
    from scheduler import scheduler

    # prepare/load data
    from instance import mk_obs_time, mk_obs_set

    # if solving multiple times the same instance, pickling preprocessed data may be convenient:
    # comment the previous lines, and uncomment the following for re-loading data quickly:
    import pickle
    filename = "telescope_sched_state.pck"
    from os.path import exists
    state = False
    if exists(filename):
        with open(filename, 'rb') as f:
            sky, time_start, time_end, obs_data = pickle.load(f)
        state = True
    time_start_, time_end_ = mk_obs_time()
    if not state or time_start != time_start_ or time_end != time_end_:  # pickled data outdated
        time_start, time_end = mk_obs_time()
        sky = mk_obs_set()
        obs_data = ObsData(sky, time_start, time_end)
        with open(filename, 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump((sky, time_start, time_end, obs_data), f, pickle.HIGHEST_PROTOCOL)

    # repeat a number of constructions
    earliest = {k: time_start for k in obs_data.visible}

    NITER = 1000
    n3obs_list = []
    best_n3obs = -1
    for i in range(NITER):
        # construct initial (preliminary) solution
        sol = scheduler(obs_data, time_now=time_start, last_obs=[None], earliest=earliest,
                        current_best=None, hidden_list=[], time_lim=0)

        n3obs_list.append(sol.n3obs)
        if sol.n3obs > best_n3obs:
            best_n3obs = sol.n3obs
            best = sol.copy()

        print("t3obs/total: {}/{}\t{}\t\t{}/{}".\
              format(sol.n3obs, len(sol.seq), sol.values, sol.n3obs, best.n3obs))

    print(n3obs_list)
    import matplotlib.pyplot as plt
    import numpy as np
    rmin = min(n3obs_list)
    rmax = max(n3obs_list)
    bins = np.arange(rmin, rmax + 1.5) - 0.5
    fig, ax = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(12, 7))
    ax.hist(n3obs_list, label="Number of 3-observations", bins=bins, alpha=0.5)
    ax.legend(loc='upper left')
    plt.grid(zorder=0)
    plt.savefig("histogram.pdf")