import numpy as np

def dr2(r1, r2):
    dr = np.array(r2) - np.array(r1)
    return np.dot(dr, dr)

def plot_data_for_a_trajectory(trajectory, plot_data):
    for i in range(1, len(trajectory)):
        if len(plot_data) <= i:
            plot_data.append([])

        for j in range(len(trajectory) - i):
            plot_data[i].append(dr2(trajectory[j][1:], trajectory[j+i][1:]))

plot_data = [[0]]

f = open('lone_sailor_trajectory.txt', 'r')
trajectory = eval(f.readline())
plot_data_for_a_trajectory(trajectory, plot_data)

for i in range(1, len(plot_data)):
    # total_data += len(plot_data[i])
    # dr2_over_dt += np.mean(plot_data[i])/i * len(plot_data[i])
    print(i, np.mean(plot_data[i]), sep = ', ')