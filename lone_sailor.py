import numpy as np
import random

steps = 10

bubble = np.zeros(3)

output = []

for t in range(100 * steps + 1):
    phi = random.random() * 2 * np.pi
    cos_theta = (random.random() - 0.5) * 2
    sin_theta = np.sqrt(1 - cos_theta**2)

    bubble[0] += sin_theta * np.cos(phi)
    bubble[1] += sin_theta * np.sin(phi)
    bubble[2] += cos_theta

    if t%steps == 0:
        output.append([t] + list(bubble))

print(output, file = open('lone_sailor_trajectory.txt', 'w'))