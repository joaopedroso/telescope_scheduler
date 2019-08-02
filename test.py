"""
test: Evaluate the nearest neighbor heuristic on some instances

Evaluates the telescope scheduler implemented in module `scheduler.py`.

"""

from astropy.time import Time
def mk_obs_times():
    """
    Define time interval for the observations
    Using Astronomical Twilight in Nagano
    https://www.timeanddate.com/sun/japan/nagano

    :returns: list of tuples with starting and ending time
    """
    utcoffset = +9 * u.hour  # JST
    times = []

    time_start = \
        [Time('2018-01-21T18:30') - utcoffset,
         Time('2018-02-21T18:59') - utcoffset,
         Time('2018-03-21T19:25') - utcoffset,
         Time('2018-04-21T19:59') - utcoffset,
         Time('2018-05-21T20:36') - utcoffset,
         Time('2018-06-21T21:00') - utcoffset,
         Time('2018-07-21T20:47') - utcoffset,
         Time('2018-08-21T20:04') - utcoffset,
         Time('2018-09-21T19:13') - utcoffset,
         Time('2018-10-21T18:30') - utcoffset,
         Time('2018-11-21T18:05') - utcoffset,
         Time('2018-12-21T18:07') - utcoffset]
    time_end = \
        [Time('2018-01-22T05:26') - utcoffset,
         Time('2018-02-22T05:03') - utcoffset,
         Time('2018-03-22T04:24') - utcoffset,
         Time('2018-04-22T03:33') - utcoffset,
         Time('2018-05-22T02:52') - utcoffset,
         Time('2018-06-22T02:37') - utcoffset,
         Time('2018-07-22T02:58') - utcoffset,
         Time('2018-08-22T03:35') - utcoffset,
         Time('2018-09-22T04:06') - utcoffset,
         Time('2018-10-22T04:33') - utcoffset,
         Time('2018-11-22T05:00') - utcoffset,
         Time('2018-12-22T05:22') - utcoffset]

    return (time_start, time_end) 



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
    filename = "telescope_sched_states.pck"
    from os.path import exists
    state = False
    if exists(filename):
        with open(filename, 'rb') as f:
            sky, time_start, time_end, obs_data = pickle.load(f)
        state = True
    time_start_, time_end_ = mk_obs_times()
    if not state or time_start != time_start_ or time_end != time_end_:  # pickled data outdated
        time_start, time_end = mk_obs_times()
        sky, obs_data = [], []
        for inst in range(len(time_start)):
            sky.append(mk_obs_set())
            obs_data.append(ObsData(sky[inst], time_start[inst], time_end[inst]))
        with open(filename, 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump((sky, time_start, time_end, obs_data), f, pickle.HIGHEST_PROTOCOL)

    for inst in range(len(time_start)):
        print("Observation limits: {}--{}".format(time_start[inst], time_end[inst]))
        earliest = {k: time_start[inst] for k in obs_data[inst].visible}
        # print(earliest[3])

        NITER = 1000
        n3obs_list = []
        best_n3obs = -1
        for i in range(NITER):
            # construct initial (preliminary) solution
            sol = scheduler(obs_data[inst], time_now=time_start[inst], last_obs=[None], earliest=earliest,
                            current_best=None, hidden_list=[], time_lim=0)
     
            n3obs_list.append(sol.n3obs)
            if sol.n3obs > best_n3obs:
                best_n3obs = sol.n3obs
                best = sol.copy()
     
            print("t3obs/total: {}/{}\t{}\t\t{}/{}".\
                  format(sol.n3obs, len(sol.seq), sol.values, sol.n3obs, best.n3obs))
     
        import matplotlib.pyplot as plt
        import numpy as np
        rmin = min(n3obs_list)
        rmax = max(n3obs_list)
        bins = np.arange(rmin, rmax + 1.5) - 0.5
        fig, ax = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(12, 7))
        ax.hist(n3obs_list, label="Number of 3-observations", bins=bins, alpha=0.5)
        ax.legend(loc='upper left')
        plt.grid(zorder=0)
        plt.title("{}--{}".format(time_start[inst], time_end[inst]))
        plt.savefig("histogram_{}.pdf".format(inst+1))
