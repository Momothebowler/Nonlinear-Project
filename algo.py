import numpy as np
import copy

alpha = 0.05
batch = 1000
delta = 0.1
n0 = 24
k = 10
c = 1
n = 1 / 2 * (((2 * alpha) / (k - 1)) ** (-2 / (n0 - 1)) - 1)
h2 = 2 * c * n * (n0 - 1)

for p in range(50):
    np.random.seed(p)
    np.set_printoptions(legacy="1.25")

    indexes = []
    data = {}
    s2 = {}
    means = []
    n_cap = {}
    presets = {}
    system_pairs = {}

    for i in range(k):
        indexes.append(i)
        s2[i] = []
        n_cap[i] = []
        presets[i] = (np.random.uniform(0.5, 1.5), np.random.uniform(1, 2))
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

        # print(indexes)

        while len(indexes) > 1 or r + 1 <= n_max + 1:
            if len(indexes) >= 1:
                index_old = copy.deepcopy(indexes)
            means = []
            r += batch
            for i in range(k):
                data[i] = np.append(
                    data[i], np.random.normal(presets[i][0], presets[i][1], batch)
                )
                means.append(sum(data[i]) / r)

            for i in range(k):
                for l in range(k):
                    if l == i or i not in indexes:
                        continue
                    w = max(0, delta / (2 * c * r) * (h2 * s2[i][l] / delta**2 - r))
                    if means[i] < means[l] - w:
                        indexes.remove(i)

        if indexes:
            try_max = 0
            try_max_index = -1
            for i in indexes:
                if try_max < means[i]:
                    try_max = max(try_max, means[i])
                    try_max_index = i
            print(
                "Max mean is: {} and it is system/s: {} with mean/std: {} and took r rounds: {} and seed: {}".format(
                    round(try_max, 2),
                    try_max_index,
                    (
                        round(presets[try_max_index][0], 2),
                        round(presets[try_max_index][1], 2),
                    ),
                    r,
                    p,
                )
            )
        else:
            try_max = 0
            try_max_index = -1
            for i in index_old:
                if try_max < means[i]:
                    try_max = max(try_max, means[i])
                    try_max_index = i
            print(
                "Max mean is: {} and it is system/s: {} with mean/std: {} and took r rounds: {} and seed: {}".format(
                    round(try_max, 2),
                    try_max_index,
                    (
                        round(presets[try_max_index][0], 2),
                        round(presets[try_max_index][1], 2),
                    ),
                    r,
                    p,
                )
            )

    else:
        try_max = 0
        try_max_index = -1
        for i in indexes:
            if try_max < means[i]:
                try_max = max(try_max, means[i])
                try_max_index = i
        print(
            "Max mean is: {} and it is system/s: {} with mean/std: {} and took r rounds: {} and seed: {}".format(
                round(try_max, 2),
                try_max_index,
                (
                    round(presets[try_max_index][0], 2),
                    round(presets[try_max_index][1], 2),
                ),
                r,
                p,
            )
        )
        print("hi")
