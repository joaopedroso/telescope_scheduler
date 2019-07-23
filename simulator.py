"""This program illustrates how to use the telescope scheduler implemented
in file `scheduler.py`.
"""
from hidden import hidden
from constants import *
from obs_data import ObsData
from solution import Solution
from scheduler import scheduler

TIME_LIM = 1   # solving time limit (in seconds)
print("starting simulation, using {} s time limit for scheduling".format(TIME_LIM))

def move_telescope(sky):
    print("telescope is now in position ({:5.1f}, {:5.1f})".format(sky.ra.degree, sky.dec.degree))

# prepare/load data
from instance import mk_obs_time, mk_obs_set
import pickle
filename = "telescope_sched_state.pck"
from os.path import exists
state = False
if exists(filename):
    with open(filename, 'rb') as f:
        sky, time_start, time_end, obs_data = pickle.load(f)
    state = True

time_start_, time_end_ = mk_obs_time()
if not state or time_start != time_start_ or time_end != time_end_: # pickled data outdated
    sky, time_start, time_end = mk_obs_set()
    obs_data = ObsData(sky, time_start, time_end)
    with open(filename, 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump((sky, time_start, time_end, obs_data), f, pickle.HIGHEST_PROTOCOL)


# start simulation of an observation night
t = time_start
earliest = {k: t for k in obs_data.K}
# construct initial (preliminary) solution
sol = scheduler(obs_data, time_now=t, last_obs=[None], earliest=earliest,
                current_best=None, hidden_list=[], time_lim=TIME_LIM)

# time_end = time_start + 120 * u.second    # use for quick tests


tpos = sol.seq[1]   # telescope position
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
    print("Not observable:", [k for k in obs_data.K if earliest[k]>t])
    sol = scheduler(obs_data, time_now=t, last_obs=seq, earliest=earliest,
                    current_best=best, hidden_list=covered, time_lim=TIME_LIM)
    if len(sol.seq) < 2:
        break
    curr = sol.seq[1]
    print("\tcurrent solution: {} ...".format(sol.seq[:10]))
    print("\tnext: {}\t3obs/total: {}/{}\tobs.card: {}".format(curr, sol.n3obs, len(sol.seq), sol.values))
    print("Simulating telescope movement to position {}: ".format(curr), end="\t")
    move_telescope(sky[curr])

    # simulate time advance
    move = obs_data.move_time[tpos,curr]
    delay = move + EXPOSURE
    t += delay * u.second
    earliest[curr] = t + MIN_OBS_GAP

    # update variables
    best = Solution(sol.seq[1:], obs_data.K)   # starting solution for next iteration
    tpos = curr
    seq.append(curr)
    print("Previously observed positions:", seq[1:])
    print()

print("SUMMARY OF OBSERVATIONS:")
final = Solution(seq, obs_data.K)
print("{}\nn.obs: {}\tpos: {}\t3obs/total: {}/{}\t{}".format(t, nobs, curr, final.n3obs, len(seq), final.values))