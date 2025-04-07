import numpy as np
import pandas as pd

np.random.seed(0)
np.set_printoptions(legacy="1.25")

alpha = 0.05
delta = 0.1
n0 = 24
k = 10
c = 1
n = 1 / 2 * (((2 * alpha) / (k - 1)) ** (-2 / (n0 - 1)) - 1)

data = {}
h2 = 2 * c * n * (n0 - 1)
s2 = {}
means = []
n_cap = {}

system_pairs = {}
for i in range(k):
    s2[i] = []
    n_cap[i] = []
    data[i] = np.random.normal(np.random.uniform(0.5, 2), np.random.uniform(1, 2), n0)
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
print(n_max)
