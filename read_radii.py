import matplotlib.pyplot as plt
import numpy as np

f = open('radii.txt', 'r')
for l in f:
    exec(l)
f.close()

def sizes__radii__radii_all(n2_min, n2_max):
    dots = 100

    x = [0, 100]
    X = [0, 100]

    for i in range(n2_min, n2_max+1):
        x.extend(radii[i])
        X.extend(radii_all[i])

    n, bins, patches = plt.hist(x, dots, facecolor='g', alpha = 0.75)
    N, Bins, patches = plt.hist(X, dots, facecolor='b', alpha = 0.75)

    return bins, n, N

sizes, N_will, N_all = sizes__radii__radii_all(2, 2)

for i in range(1, 99):
    if N_will[i] != 0:
        print(np.log(sizes[i]), np.log(N_will[i] / N_all[i]), sep = ', ')