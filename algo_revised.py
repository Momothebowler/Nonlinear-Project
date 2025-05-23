import numpy as np
import pandas as pd
from collections import Counter

runs = 10
alpha = 0.05
batch = 1
delta = 0.5 #delta has to be smaller if we want to actually calculate the largest mean; This will cause it run a lot longer...
n0 = 24
k = 10
c = 1
n = 1 / 2 * (((2 * alpha) / (k - 1)) ** (-2 / (n0 - 1)) - 1)
h2 = 2 * c * n * (n0 - 1)

#fix seed for testing
np.random.seed(0)

#generate the means/variances out here so that the for loop will handle the same systems repeatedly
indexes = []
presets = {}
FINAL = []
for i in range(k):
    indexes.append(i)
    presets[i] = (np.random.uniform(0.5, 2), np.random.uniform(12, 15)) #may have to increase range on the uniform if we use a small delta

actual_means, actual_variances = zip(*list(presets.values()))
true_max_index = actual_means.index(max(actual_means))
lower_bound = actual_means[true_max_index] - delta
acceptable_indexes = []
for i in range(k):
    if actual_means[i] >= lower_bound:
        acceptable_indexes.append(i)

for p in range(runs):
    #np.random.seed(p)
    #np.set_printoptions(legacy="1.25")

    n0 = 24
    indexes = []
    data = {}
    s2 = {}
    means = []
    n_cap = {}
    #presets = {}
    system_pairs = {}

    for i in range(k):
        indexes.append(i)
        s2[i] = []
        n_cap[i] = []
        #presets[i] = (np.random.uniform(0.5, 2), np.random.uniform(1, 2))
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

    # df = pd.DataFrame(n_cap)
    # print(n_max)

    if n0 < n_max:
        r = n0
        for i in range(k):
            for l in range(k):
                if l == i or i not in indexes:
                    continue
                w = max(0, delta / (2 * c * r) * (h2 * s2[i][l] / delta**2 - r))
                if means[i] < means[l] - w:
                    indexes.remove(i)

        # print(indexes)

        while len(indexes) > 1 and r <= n_max + 1: #used `or` originally
            means = []
            r += 1 #not batch
            for i in range(k):
                data[i] = np.append(data[i], np.random.normal(presets[i][0], presets[i][1], 1))#[0])
                means.append(sum(data[i]) / r)

            old_indexes = indexes[:] #you have to do this in case the upcoming for loop removes everything!
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
                    round(means[indexes[0]], 2), #you originally had max(means), but what if the maximum isn't the one?
                    indexes[0],
                    (round(presets[indexes[0]][0], 2), round(presets[indexes[0]][1], 2)),
                    r,
                )
            )
        elif len(indexes) == 0: #in case the screening removes all the indeces
            indexes = old_indexes
            print("0! Go back to old indexes!")
            print(indexes)
            means_candidates = [means[i] for i in indexes]
            A = means_candidates.index(max(means_candidates))
            FINAL.append(indexes[A])
            print(
                "Max mean is: {} and it is system/s: {} with mean/std: {} and took r rounds: {}".format(
                    round(means_candidates[A], 2), #you originally had max(means), but what if the maximum isn't the one?
                    indexes[A],
                    (round(presets[indexes[A]][0], 2), round(presets[indexes[A]][1], 2)),
                    r,
                )
            )
        else: #This handles the case when there is more than 1 index to choose from
            print("More than 1!")
            means_candidates = [means[i] for i in indexes]
            A = means_candidates.index(max(means_candidates))
            FINAL.append(indexes[A])
            print(
                "Max mean is: {} and it is system/s: {} with mean/std: {} and took r rounds: {}".format(
                    round(means_candidates[A], 2), #you originally had max(means), but what if the maximum isn't the one?
                    indexes[A],
                    (round(presets[indexes[A]][0], 2), round(presets[indexes[A]][1], 2)),
                    r,
                )
            )

    else: #used 0 before
        
        A = means.index(max(means))
        print(
            "Max mean is: {} and it is system/s: {} with mean/std: {}".format(
                round(means[A], 2), A, presets[A]
            )
        )

print(Counter(FINAL))
print("The true maximum mean was {} given by system {}".format(actual_means[true_max_index], true_max_index))
print(acceptable_indexes)
print(FINAL)
