"""
solution: module to hold and preprocess a solution

classes:
    - Solution(): convenience class for holding all information concerning a schedule
"""

from copy import deepcopy

class Solution:
    def __init__(self, seq, K):
        self.seq = seq.copy()
        self.obs = {k: 0 for k in K}
        for k in seq[1:]:
            self.obs[k] += 1
        self.n_max = max(self.obs[k] for k in K) + 1
        self.values = {v: sum(1 for k in K if self.obs[k] == v) for v in range(self.n_max)}
        self.n3obs = sum(self.values[i] for i in range(3, self.n_max))

    def copy(self):
        return deepcopy(self)