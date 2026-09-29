"""Single-profile, literal O(T K²) reference for FastJointSLMLibraryI.f.

Matrices use Python orientation: transitions[t, source, destination],
emissions/scores/predecessors[t, state]. Paths and predecessors remain one-based.
The Fortran PI literal rounds to binary32 before promotion to binary64.
Keep that value and the original log-sum evaluation order.
"""

import math

import numpy as np


def elnsum(x: float, y: float) -> float:
    return x + math.log(1 + math.exp(y - x)) if x > y else y + math.log(1 + math.exp(x - y))


def matrices(values, means, mi, smu, sepsilon, eta):
    k = len(means)
    g = -((means - mi) ** 2 / (2 * smu**2))
    norm = float(g[0])
    for value in g[1:]:
        norm = elnsum(norm, float(value))
    g = g - norm
    transitions = np.empty((len(eta), k, k))
    for t, change in enumerate(eta):
        log_change = math.log(change)
        log_stay = math.log(1 - change) if change < 1 else -math.inf
        for source in range(k):
            for destination in range(k):
                jump = log_change + g[destination]
                transitions[t, source, destination] = (
                    elnsum(log_stay, jump) if source == destination else jump
                )
    emissions = np.empty((len(values), k))
    constant = math.log(1 / (math.sqrt(2 * 3.1415927410125732421875) * sepsilon))
    for t, value in enumerate(values):
        for j, mean in enumerate(means):
            emissions[t, j] = constant + -0.5 * ((value - mean) / sepsilon) ** 2
    return transitions, emissions


def viterbi(initial, transitions, emissions):
    n, k = emissions.shape
    scores = np.empty((n, k))
    predecessors = np.zeros((n, k), dtype=np.int32)
    scores[0] = initial + emissions[0]
    for t in range(1, n):
        for j in range(k):
            best = scores[t - 1, 0] + transitions[t - 1, 0, j]
            index = 0
            for source in range(1, k):
                candidate = scores[t - 1, source] + transitions[t - 1, source, j]
                if candidate > best:
                    best, index = candidate, source
            predecessors[t, j] = index + 1
            scores[t, j] = best + emissions[t, j]
    path = np.empty(n, dtype=np.int32)
    path[-1] = np.argmax(scores[-1]) + 1
    for t in range(n - 2, -1, -1):
        path[t] = predecessors[t + 1, path[t + 1] - 1]
    return path, scores, predecessors
