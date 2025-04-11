import numpy as np
import pandas as pd
from collections import Counter

alpha = 0.05
batch = 1
delta = 0.1
n0 = 24
k = 10
c = 1
n = 1 / 2 * (((2 * alpha) / (k - 1)) ** (-2 / (n0 - 1)) - 1)
h2 = 2 * c * n * (n0 - 1)


np.random.seed(0)

indexes = []
presets = {}
FINAL = []
for i in range(k):
    indexes.append(i)
    presets[i] = (np.random.uniform(0.5, 2), np.random.uniform(1, 2))

actual_means, actual_variances = zip(*list(presets.values()))
true_max_index = actual_means.index(max(actual_means))
lower_bound = actual_means[true_max_index] - delta
acceptable_indexes = []

# this does not work as intended
for i in range(k):
    if actual_means[i] >= lower_bound:
        acceptable_indexes.append(i)

for p in range(100):
    n0 = 24
    indexes = []
    data = {}
    s2 = {}
    means = []
    n_cap = {}
    system_pairs = {}

    for i in range(k):
        indexes.append(i)
        s2[i] = []
        n_cap[i] = []
        data[i] = np.random.normal(presets[i][0], presets[i][1], n0)
        means.append(sum(data[i]) / n0)

    for l in range(k):
        for i in range(k):
            s2[i].append(0)
            n_cap[i].append(0)
            for j in range(n0):
                s2[i][l] += (data[i][j] - data[l][j] - (means[i] - means[l])) ** 2
            s2[i][l] *= 1 / (n0 - 1)
            n_cap[i][l] = np.floor((s2[i][l] * h2) / (delta**2))

    n_i = []
    if n_i == []:
        for i in range(k):
            n_i.append(max(n_cap[i]))
        n_max = max(n_i)

    if n0 < n_max:
        r = n0
        for i in range(k):
            for l in range(k):
                if l == i or i not in indexes:
                    continue
                w = max(0, delta / (2 * c * r) * (h2 * s2[i][l] / delta**2 - r))
                if means[i] < means[l] - w:
                    indexes.remove(i)

        while len(indexes) > 1 and r <= n_max + 1:
            means = []
            r += 1
            for i in range(k):
                data[i] = np.append(
                    data[i], np.random.normal(presets[i][0], presets[i][1], 1)
                )
                means.append(sum(data[i]) / r)

            old_indexes = indexes[:]
            for i in range(k):
                for l in range(k):
                    if l == i or i not in indexes:
                        continue
                    w = max(0, delta / (2 * c * r) * (h2 * s2[i][l] / delta**2 - r))
                    if means[i] < means[l] - w:
                        indexes.remove(i)

        if len(indexes) == 1:
            FINAL.append(indexes[0])
            print(
                "Max mean is: {} and it is system/s: {} with mean/std: {} and took r rounds: {}".format(
                    round(means[indexes[0]], 2),
                    indexes[0],
                    (
                        round(presets[indexes[0]][0], 2),
                        round(presets[indexes[0]][1], 2),
                    ),
                    r,
                )
            )
        elif len(indexes) == 0:
            indexes = old_indexes
            print("0! Go back to old indexes!")
            print(indexes)
            means_candidates = [means[i] for i in indexes]
            A = means_candidates.index(max(means_candidates))
            FINAL.append(indexes[A])
            print(
                "Max mean is: {} and it is system/s: {} with mean/std: {} and took r rounds: {}".format(
                    round(means_candidates[A], 2),
                    indexes[A],
                    (
                        round(presets[indexes[A]][0], 2),
                        round(presets[indexes[A]][1], 2),
                    ),
                    r,
                )
            )
        else:
            print("More than 1!")
            means_candidates = [means[i] for i in indexes]
            A = means_candidates.index(max(means_candidates))
            FINAL.append(indexes[A])
            print(
                "Max mean is: {} and it is system/s: {} with mean/std: {} and took r rounds: {}".format(
                    round(means_candidates[A], 2),
                    indexes[A],
                    (
                        round(presets[indexes[A]][0], 2),
                        round(presets[indexes[A]][1], 2),
                    ),
                    r,
                )
            )

    else:

        A = means.index(max(means))
        print(
            "Max mean is: {} and it is system/s: {} with mean/std: {}".format(
                round(means[A], 2), A, presets[A]
            )
        )

print(Counter(FINAL))
print(
    "The true maximum mean was {} given by system {}".format(
        actual_means[true_max_index], true_max_index
    )
)
print(acceptable_indexes)
