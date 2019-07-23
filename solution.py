"""
solution: Hold and preprocess a scheduler's solution

classes:
    - `Solution()`: convenience class for holding all information concerning an observation schedule
"""

from copy import deepcopy


class Solution:
    def __init__(self, seq, K):
        """
        Convenience class holding parameters for a solution's information

        Data to be stored in self:

        :ivar seq: list with the sequence of positions being observed (first element is previous telescope position)
        :ivar obs: dictionary, `obs[k]` -> number of observations of position `obs_data.sky[k]`
        :ivar n_max: maximum number of times some position is observed
        :ivar values: dictionary, `values[i]` -> number of positions observed `i` times
        :ivar n3obs: total number of positions observed at least 3 times (the main objective of the problem)

        :param seq: list with sequence of indices from `obs_data.sky` being observed
        :param K: length of `obs_data.sky` array
        """

        self.seq = seq.copy()
        self.obs = {k: 0 for k in K}
        for k in seq[1:]:
            self.obs[k] += 1
        self.n_max = max(self.obs[k] for k in K) + 1
        self.values = {v: sum(1 for k in K if self.obs[k] == v) for v in range(self.n_max)}
        self.n3obs = sum(self.values[i] for i in range(3, self.n_max))

    def copy(self):
        return deepcopy(self)
